//% @file MathFunctions.scala
//
// Math functions module
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import chisel3._
import chisel3.util._

import rial.approx._
import rial.arith.{RealSpec, RoundSpec}
import rial.util.PipelineStageConfig

import scala.collection.mutable.ArrayBuffer
import scala.reflect.ClassTag

/** A module that calculates math functions
 *
 * @constructor Constructs a new MathFunctions module.
 * @param cfg       Configuration that determines which functions are supported.
 * @param spec      RealSpec for input/output floating point numbers.
 * @param stage     Pipeline-stage configuration (each phase configurable).
 */
class MathFunctions(
  val cfg:  MathFuncConfig,
  val spec: RealSpec,
  val stage: MathFuncPipelineConfig
) extends Module {

  import FuncKind._

  assert(!spec.disableSign, "assuming the float type allows negative values")

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val bitW   = cfg.bitW
  val coeffs = cfg.coeffs

  val io = IO(new Bundle {
    val sel = Input (UInt(cfg.signalW.W))
    val x   = Input (UInt(spec.W.W))
    val z   = Output(UInt(spec.W.W))
  })

  val selIn   = ShiftRegister(io.sel,                 1, 0.U, true.B)
  val selTab  = ShiftRegister(selIn,   stage.preStage  , 0.U, true.B)
  val selCalc = ShiftRegister(selTab,  stage.tableStage, 0.U, true.B)
  val selPost = ShiftRegister(selCalc, stage.calcStage , 0.U, true.B)
  val selOut  = ShiftRegister(selPost, stage.postStage , 0.U, true.B)

  // ==========================================================================
  // Output selection
  // ==========================================================================

  val zs = cfg.funcs.map(fn => {
    fn -> WireDefault(0.U.asTypeOf(Valid(UInt(spec.W.W))))
  }).toMap

  io.z := MuxCase(0.U, cfg.funcs.map(fn => {
    (selOut === cfg.signal(fn) && zs(fn).valid) -> zs(fn).bits
  }))

  // ==========================================================================
  // Decompose input and create delayed versions
  // ==========================================================================

  val xdecomp = Module(new DecomposeReal(spec))
  xdecomp.io.real.valid := cfg.isValid(io.sel)
  xdecomp.io.real.bits  := io.x

  val xIn   = Pipe(xdecomp.io.decomp, 1)
  val xTab  = Pipe(xIn,   stage.preStage  )
  val xCalc = Pipe(xTab,  stage.tableStage)
  val xPost = Pipe(xCalc, stage.calcStage )
  val xOut  = Pipe(xPost, stage.postStage )

  // ==========================================================================
  // Shared pre/post-process multipliers (only built if some configured
  // function needs them; sized to the minimum width that satisfies all
  // such functions).
  // ==========================================================================

  val preMul: Option[PreProcMultiplier] =
    if (PreProcMultiplier.isNeeded(cfg.funcs)) {
      Some(Module(new PreProcMultiplier(
        lhsW  = PreProcMultiplier.lhsW(cfg.funcs, spec),
        manW  = manW,
        stage = stage.preMulStage,
      )))
    } else { None }

  // Pre-mult driver wires: each function that wants the multiplier overrides
  // these in its own when-block below. Default = 0 (no spurious switching).
  val (preMulValid, preMulLhs, preMulRhs) = if (PreProcMultiplier.isNeeded(cfg.funcs)) {
    (
      Some(WireDefault(false.B)),
      Some(WireDefault(0.U(PreProcMultiplier.lhsW(cfg.funcs, spec).W))),
      Some(WireDefault(0.U((1+manW).W)))
    )
  } else {(None, None, None)}

  if (preMul.isDefined) {
    preMul.get.io.in.valid    := preMulValid.get
    preMul.get.io.in.bits.lhs := preMulLhs.get
    preMul.get.io.in.bits.rhs := preMulRhs.get
  }

  // --------------------------------------------------------------------------

  val postMul: Option[PostProcMultiplier] =
    if (PostProcMultiplier.isNeeded(cfg.funcs)) {
      Some(Module(new PostProcMultiplier(
        fracW     = bitW.outBitW,
        manW      = manW,
        roundSpec = RoundSpec.roundToEven,
        stage     = stage.postMulStage,
      )))
    } else { None }

  val (postMulValid, postMulLhs, postMulRhs) = if (PostProcMultiplier.isNeeded(cfg.funcs)) {
    (
      Some(WireDefault(false.B)),
      Some(WireDefault(0.U((1 + bitW.outBitW).W))),
      Some(WireDefault(0.U((1+manW).W)))
    )
  } else {(None, None, None)}

  if (postMul.isDefined) {
    postMul.get.io.in.valid     := postMulValid.get
    postMul.get.io.in.bits.lhs  := postMulLhs.get
    postMul.get.io.in.bits.rhs  := postMulRhs.get
  }

  // ==========================================================================
  // Polynomial Evaluator (shared)
  // ==========================================================================

  val polynomialEval = Module(new PolynomialEval(bitW, stage.calcStage))

  polynomialEval.io.in.valid := false.B

  val polynomialCoefW = bitW.coefBits.reduce(_+_)
  val polynomialCoefs = cfg.funcs.map(fn => {
    fn -> WireDefault(0.U.asTypeOf(Valid(UInt(polynomialCoefW.W))))
  }).toMap

  val polynomialDxs = if (bitW.nOrder == 0) { None } else { Some(
    cfg.funcs.map(fn => {
      fn -> WireDefault(0.U.asTypeOf(Valid(UInt(bitW.inBitW.W))))
    }).toMap
  )}

  val polynomialCoef = MuxCase(0.U.asTypeOf(Valid(UInt(polynomialCoefW.W))),
    cfg.funcs.map(fn => {
      (selCalc === cfg.signal(fn)) -> polynomialCoefs(fn)
    })
  )
  val polynomialDx = if (bitW.nOrder == 0) { None } else { Some(
    MuxCase(0.U.asTypeOf(Valid(UInt(bitW.inBitW.W))),
      cfg.funcs.map(fn => {
        (selCalc === cfg.signal(fn)) -> polynomialDxs.get(fn)
      })
    )
  )}

  if(bitW.nOrder == 0) { // w/o dx, table only
    polynomialEval.io.in.valid       := polynomialCoef.valid
    polynomialEval.io.in.bits.coeffs := polynomialCoef.bits.asTypeOf(
      MixedVec(bitW.coefBits.map(w => SInt(w.W)))
    )
  } else { // w/ dx
    polynomialEval.io.in.valid       := polynomialCoef.valid && polynomialDx.get.valid
    polynomialEval.io.in.bits.coeffs := polynomialCoef.bits.asTypeOf(
      MixedVec(bitW.coefBits.map(w => SInt(w.W)))
    )
    polynomialEval.io.in.bits.x.get  := polynomialDx.get.bits.asSInt
  }

  // ==========================================================================
  // Per-function pre/tab/nontab/post modules
  //
  // `.collect` so that Sin/Cos (which share their HW) are skipped here and
  // handled by the dedicated block below.
  // ==========================================================================

  val pres: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt       => Sqrt       -> Module(new SqrtPreProcess      (spec, bitW, stage.preStage))
    case InvSqrt    => InvSqrt    -> Module(new InvSqrtPreProcess   (spec, bitW, stage.preStage))
    case Reciprocal => Reciprocal -> Module(new ReciprocalPreProcess(spec, bitW, stage.preStage))
  }.toMap[FuncKind, Module]

  val tabs: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt       => Sqrt       -> Module(new SqrtTableCoeff      (spec, coeffs(Sqrt      ), stage.tableStage))
    case InvSqrt    => InvSqrt    -> Module(new InvSqrtTableCoeff   (spec, coeffs(InvSqrt   ), stage.tableStage))
    case Reciprocal => Reciprocal -> Module(new ReciprocalTableCoeff(spec, coeffs(Reciprocal), stage.tableStage))
  }.toMap

  val nontabs: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt       => Sqrt       -> Module(new SqrtNonTable      (spec, bitW, stage.nonTableStage))
    case InvSqrt    => InvSqrt    -> Module(new InvSqrtNonTable   (spec, bitW, stage.nonTableStage))
    case Reciprocal => Reciprocal -> Module(new ReciprocalNonTable(spec, bitW, stage.nonTableStage))
  }.toMap

  val posts: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt       => Sqrt       -> Module(new SqrtPostProcess      (spec, bitW, stage.postStage))
    case InvSqrt    => InvSqrt    -> Module(new InvSqrtPostProcess   (spec, bitW, stage.postStage))
    case Reciprocal => Reciprocal -> Module(new ReciprocalPostProcess(spec, bitW, stage.postStage))
  }.toMap

  private def toModule[T <: Module](x: Module)(implicit ct: ClassTag[T]): T = {
    require(ct.runtimeClass.isInstance(x),
      s"expected ${ct.runtimeClass.getSimpleName} but got ${x.getClass.getSimpleName}")
    x.asInstanceOf[T]
  }

  // --------------------------------------------------------------------------

  if(cfg.has(Sqrt)) {
    val pre  = toModule[SqrtPreProcess ](pres   (Sqrt))
    val tab  = toModule[SqrtTableCoeff ](tabs   (Sqrt))
    val ntab = toModule[SqrtNonTable   ](nontabs(Sqrt))
    val post = toModule[SqrtPostProcess](posts  (Sqrt))
    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(Sqrt) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Sqrt) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    ntab.io.x := xTab

    post.io.in.valid     := ntab.io.z.valid && polynomialEval.io.out.valid
    post.io.in.bits.info := ntab.io.z.bits
    post.io.in.bits.zres := polynomialEval.io.out.bits

    zs(Sqrt) := post.io.z
  }

  // --------------------------------------------------------------------------

  if(cfg.has(InvSqrt)) {
    val pre  = toModule[InvSqrtPreProcess ](pres   (InvSqrt))
    val tab  = toModule[InvSqrtTableCoeff ](tabs   (InvSqrt))
    val ntab = toModule[InvSqrtNonTable   ](nontabs(InvSqrt))
    val post = toModule[InvSqrtPostProcess](posts  (InvSqrt))
    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(InvSqrt) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(InvSqrt) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    ntab.io.x := xTab

    post.io.in.valid     := ntab.io.z.valid && polynomialEval.io.out.valid
    post.io.in.bits.info := ntab.io.z.bits
    post.io.in.bits.zres := polynomialEval.io.out.bits

    zs(InvSqrt) := post.io.z
  }

  // --------------------------------------------------------------------------

  if(cfg.has(Reciprocal)) {
    val pre  = toModule[ReciprocalPreProcess ](pres   (Reciprocal))
    val tab  = toModule[ReciprocalTableCoeff ](tabs   (Reciprocal))
    val ntab = toModule[ReciprocalNonTable   ](nontabs(Reciprocal))
    val post = toModule[ReciprocalPostProcess](posts  (Reciprocal))
    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(Reciprocal) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Reciprocal) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    ntab.io.x := xTab

    post.io.in.valid     := ntab.io.z.valid && polynomialEval.io.out.valid
    post.io.in.bits.info := ntab.io.z.bits
    post.io.in.bits.zres := polynomialEval.io.out.bits

    zs(Reciprocal) := post.io.z
  }

  // --------------------------------------------------------------------------
  // Sin / Cos -- one shared datapath, distinguished by `isSin` at runtime.
  //
  // Pipeline alignment:
  //   preMul       : 0 .. preMulStage         (consumes xIn at t=0)
  //   sincosPre    : preMulStage .. preStage  (consumes preMul output)
  //   sincosTab    : preStage    .. preStage+tableStage
  //   polynomialEval (shared) and sincosPre.info pipelined to here
  //   postMul      : preStage+tableStage+calcStage .. + postMulStage
  //   sincosPost   :   ...                    .. postStage (final)

  if(cfg.has(Sin) || cfg.has(Cos)) {
    assert(preMul.isDefined,  "Sin/Cos requires the pre-process multiplier")
    assert(postMul.isDefined, "Sin/Cos requires the post-process multiplier")

    // a single SinCos pre/tab/post; sel signal isSin chooses between them
    val sincosPreStage  = stage.preStage  - stage.preMulStage
    val sincosPostStage = stage.postStage - stage.postMulStage

    val sincosPre  = Module(new SinCosPreProcess (spec, bitW, sincosPreStage))
    val sincosTab  = Module(new SinCosTableCoeff (spec, coeffs(Sin), stage.tableStage))
    val sincosPost = Module(new SinCosPostProcess(spec, bitW, sincosPostStage))

    val isSinSelIn    = selIn === cfg.signal(Sin)
    val isCosSelIn    = selIn === cfg.signal(Cos)
    val sincosSelIn   = isSinSelIn || isCosSelIn

    // ----- pre-mult driver: 1/π * (1 + x.man) -----
    when (sincosSelIn) {
      preMulLhs.get   := SinCosPreMulArgs.lhs(spec, preMul.get.lhsW)
      preMulRhs.get   := Cat(1.U(1.W), xIn.bits.man)
      preMulValid.get := xIn.valid
    }

    // ----- pre stage: x delayed by preMulStage to align with preMul.out -----
    val xAtSinCosPre  = Pipe(xIn, stage.preMulStage)
    val isSinAtPre    = ShiftRegister(isSinSelIn, stage.preMulStage, false.B, true.B)

    sincosPre.io.isSin := isSinAtPre
    sincosPre.io.x     := xAtSinCosPre
    sincosPre.io.prod  := preMul.get.io.out

    // ----- table -----
    sincosTab.io.adr := sincosPre.io.adr

    polynomialCoefs(Sin) := sincosTab.io.cs
    polynomialCoefs(Cos) := sincosTab.io.cs
    if(bitW.nOrder != 0) {
      val dxAtCalc = Pipe(sincosPre.io.dx.get, stage.tableStage)
      polynomialDxs.get(Sin) := dxAtCalc
      polynomialDxs.get(Cos) := dxAtCalc
    }

    // ----- post-mult driver: polynomial * (1 + yman) -----
    //
    // sincosPre.io.info is at xTab timing; the polynomial result emerges at
    // xPost timing (= xTab + tableStage + calcStage). To feed both into the
    // post-mult at the same cycle, pipe `info` along.
    val sincosInfoAtPost = Pipe(sincosPre.io.info, stage.tableStage + stage.calcStage)
    val sincosSelAtPost  = ShiftRegister(sincosSelIn,
      stage.preStage + stage.tableStage + stage.calcStage, false.B, true.B)

    when (sincosSelAtPost) {
      postMulLhs.get := Cat(1.U(1.W), polynomialEval.io.out.bits)
      postMulRhs.get := Cat(1.U(1.W), sincosInfoAtPost.bits.yman)
      postMulValid.get := polynomialEval.io.out.valid && sincosInfoAtPost.valid
    }

    // ----- post -----
    //
    // sincosPost consumes the pre's `info` (yex etc.) at xPost+postMulStage
    // timing -- same as postMul.io.out.
    sincosPost.io.pre := Pipe(sincosInfoAtPost, stage.postMulStage)
    sincosPost.io.mul := postMul.get.io.out

    zs(Sin) := sincosPost.io.z
    zs(Cos) := sincosPost.io.z
  }
}
