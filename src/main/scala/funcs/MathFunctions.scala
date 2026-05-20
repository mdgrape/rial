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
    val y   = if (cfg.has(FuncKind.ATan2Phase1)) Some(Input(UInt(spec.W.W))) else None
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

  // ATan2 also needs the y operand decomposed/pipelined.
  val ydecomp: Option[DecomposeReal] = if (cfg.has(ATan2Phase1)) {
    val m = Module(new DecomposeReal(spec))
    m.io.real.valid := cfg.isValid(io.sel)
    m.io.real.bits  := io.y.get
    Some(m)
  } else { None }

  val yIn: Option[Valid[DecomposedReal]] = if (cfg.has(ATan2Phase1)) {
    Some(Pipe(ydecomp.get.io.decomp, 1))
  } else { None }

  // ==========================================================================
  // Shared pre/post-process multipliers (only built if some configured
  // function needs them; sized to the minimum width that satisfies all
  // such functions).
  // ==========================================================================

  val needPreMul = PreProcMultiplier.isNeeded(cfg.funcs)

  val preMul: Option[PreProcMultiplier] = if (needPreMul) { Some(
    Module(new PreProcMultiplier(
      lhsW  = PreProcMultiplier.lhsW(cfg.funcs, spec),
      manW  = manW,
      stage = stage.preMulStage,
    ))
  )} else { None }

  // Pre-mult driver wires: each function that wants the multiplier overrides
  // these in its own when-block below. Default = 0 (no spurious switching).
  val (preMulValid, preMulLhs, preMulRhs) = if (needPreMul) { (
    Some(WireDefault(false.B)),
    Some(WireDefault(0.U(PreProcMultiplier.lhsW(cfg.funcs, spec).W))),
    Some(WireDefault(0.U((1+manW).W)))
  )} else {(None, None, None)}

  if (preMul.isDefined) {
    preMul.get.io.in.valid    := preMulValid.get
    preMul.get.io.in.bits.lhs := preMulLhs.get
    preMul.get.io.in.bits.rhs := preMulRhs.get
  }

  // --------------------------------------------------------------------------

  val needPostMul = PostProcMultiplier.isNeeded(cfg.funcs)

  val postMul: Option[PostProcMultiplier] = if (needPostMul) { Some(
    Module(new PostProcMultiplier(
      fracW     = bitW.outBitW,
      manW      = manW,
      roundSpec = RoundSpec.roundToEven,
      stage     = stage.postMulStage,
    ))
  )} else { None }

  val (postMulValid, postMulLhs, postMulRhs) = if (needPostMul) { (
    Some(WireDefault(false.B)),
    Some(WireDefault(0.U((1 + bitW.outBitW).W))),
    Some(WireDefault(0.U((1+manW).W)))
  )} else {(None, None, None)}

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
  // `.collect` so that some funcs does not have the corresponding modules
  // ==========================================================================

  val pres: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtPreProcess       (spec, bitW, stage.preStage))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtPreProcess    (spec, bitW, stage.preStage))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalPreProcess (spec, bitW, stage.preStage))
    case Sin         => Sin         -> Module(new SinCosPreProcess     (spec, bitW, stage.preStage-stage.preMulStage))
    case Exp         => Exp         -> Module(new ExpPreProcess        (spec, bitW, stage.preStage-stage.preMulStage))
    case Log         => Log         -> Module(new LogPreProcess        (spec, bitW, stage.preStage))
    case ACosPhase1  => ACosPhase1  -> Module(new ACosPhase1PreProcess (spec, bitW, stage.preStage))
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2PreProcess (spec, bitW, stage.preStage))
    case ATan2Phase1 => ATan2Phase1 -> Module(new ATan2Phase1PreProcess(spec, bitW, stage.preStage))
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2PreProcess(spec, bitW, stage.preStage))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidPreProcess    (spec, bitW, stage.preStage))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusPreProcess   (spec, bitW, stage.preStage))
    case ScaleMixtureGaussian => {
      val (sA, sB) = cfg.scaleMixtureGaussianSigma.get
      ScaleMixtureGaussian -> Module(new ScaleMixtureGaussianPreProcess(
        sA, sB, spec, bitW, stage.preStage - stage.preMulStage))
    }
  }.toMap[FuncKind, Module]

  val tabs: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtTableCoeff       (spec, coeffs(Sqrt       ), stage.tableStage))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtTableCoeff    (spec, coeffs(InvSqrt    ), stage.tableStage))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalTableCoeff (spec, coeffs(Reciprocal ), stage.tableStage))
    case Sin         => Sin         -> Module(new SinCosTableCoeff     (spec, coeffs(Sin        ), stage.tableStage))
    case Log         => Log         -> Module(new LogTableCoeff        (spec, coeffs(Log        ), stage.tableStage))
    case Exp         => Exp         -> Module(new ExpTableCoeff        (spec, coeffs(Exp        ), stage.tableStage))
    // ACosPhase1 shares Sqrt's TableCoeff
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2TableCoeff (spec, coeffs(ACosPhase2 ), stage.tableStage))
    // ATan2Phase1 shares Reciprocal's TableCoeff
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2TableCoeff(spec, coeffs(ATan2Phase2), stage.tableStage))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidTableCoeff    (spec, coeffs(Sigmoid    ), stage.tableStage))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusTableCoeff   (spec, coeffs(SoftPlus   ), stage.tableStage))
    case ScaleMixtureGaussian => ScaleMixtureGaussian ->
      Module(new ScaleMixtureGaussianTableCoeff(spec, coeffs(ScaleMixtureGaussian), stage.tableStage))
  }.toMap

  val nontabs: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt       => Sqrt       -> Module(new SqrtNonTable      (spec, bitW, stage.nonTableStage))
    case InvSqrt    => InvSqrt    -> Module(new InvSqrtNonTable   (spec, bitW, stage.nonTableStage))
    case Reciprocal => Reciprocal -> Module(new ReciprocalNonTable(spec, bitW, stage.nonTableStage))
    // sincos / acos do not have a dedicated non-table path
    case Exp        => Exp        -> Module(new ExpNonTable       (spec, bitW, stage.nonTableStage))
    case Log        => Log        -> Module(new LogNonTable       (spec, bitW, stage.nonTableStage))
    // ACosPhase1 reuses SqrtNonTable
    // ACosPhase2 doesn't need one.
  }.toMap

  val posts: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtPostProcess       (spec, bitW, stage.postStage))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtPostProcess    (spec, bitW, stage.postStage))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalPostProcess (spec, bitW, stage.postStage))
    case Sin         => Sin         -> Module(new SinCosPostProcess     (spec, bitW, stage.postStage-stage.postMulStage))
    case Exp         => Exp         -> Module(new ExpPostProcess        (spec, bitW, stage.postStage))
    case Log         => Log         -> Module(new LogPostProcess        (spec, bitW, stage.postStage-stage.postMulStage))
    case ACosPhase1  => ACosPhase1  -> Module(new ACosPhase1PostProcess (spec, bitW, stage.postStage))
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2PostProcess (spec, bitW, stage.postStage-stage.postMulStage))
    case ATan2Phase1 => ATan2Phase1 -> Module(new ATan2Phase1PostProcess(spec, bitW, stage.postStage-stage.postMulStage))
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2PostProcess(spec, bitW, stage.postStage-stage.postMulStage))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidPostProcess    (spec, bitW, stage.postStage))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusPostProcess   (spec, bitW, stage.postStage))
    case ScaleMixtureGaussian => ScaleMixtureGaussian ->
      Module(new ScaleMixtureGaussianPostProcess(spec, bitW, stage.postStage - stage.postMulStage))
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
  // Sin / Cos

  if(cfg.has(Sin) || cfg.has(Cos)) {
    assert(preMul .isDefined, "Sin/Cos requires the pre-process multiplier")
    assert(postMul.isDefined, "Sin/Cos requires the post-process multiplier")

    val sincosPre  = toModule[SinCosPreProcess ](pres (Sin))
    val sincosTab  = toModule[SinCosTableCoeff ](tabs (Sin))
    val sincosPost = toModule[SinCosPostProcess](posts(Sin))

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
    val sincosInfoAtPost = Pipe(sincosPre.io.info, stage.tableStage + stage.calcStage)
    val sincosSelAtPost  = ShiftRegister(sincosSelIn,
      stage.preStage + stage.tableStage + stage.calcStage, false.B, true.B)

    when (sincosSelAtPost) {
      postMulLhs.get := Cat(1.U(1.W), polynomialEval.io.out.bits)
      postMulRhs.get := Cat(1.U(1.W), sincosInfoAtPost.bits.yman)
      postMulValid.get := polynomialEval.io.out.valid && sincosInfoAtPost.valid
    }

    // ----- post -----
    sincosPost.io.pre := Pipe(sincosInfoAtPost, stage.postMulStage)
    sincosPost.io.mul := postMul.get.io.out

    zs(Sin) := sincosPost.io.z
    zs(Cos) := sincosPost.io.z
  }

  // --------------------------------------------------------------------------
  // Exp

  if(cfg.has(Exp)) {
    assert(preMul.isDefined, "Exp requires the pre-process multiplier")

    val pre  = toModule[ExpPreProcess ](pres   (Exp))
    val tab  = toModule[ExpTableCoeff ](tabs   (Exp))
    val ntab = toModule[ExpNonTable   ](nontabs(Exp))
    val post = toModule[ExpPostProcess](posts  (Exp))

    val expSelIn = selIn === cfg.signal(Exp)

    // ----- pre-mult driver: log2(e) * (1 + x.man) -----
    when (expSelIn) {
      preMulLhs.get   := ExpPreMulArgs.lhs(spec, preMul.get.lhsW)
      preMulRhs.get   := Cat(1.U(1.W), xIn.bits.man)
      preMulValid.get := xIn.valid
    }

    // ----- pre stage -----
    val xAtExpPre  = Pipe(xIn, stage.preMulStage)

    pre.io.x    := xAtExpPre
    pre.io.prod.valid := preMul.get.io.out.valid
    pre.io.prod.bits  := ExpPreMulArgs.prod(spec, preMul.get.io.out.bits)

    // ----- table -----
    tab.io.adr := pre.io.adr

    polynomialCoefs(Exp) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Exp) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    // ----- non-table ------
    ntab.io.x    := xTab
    ntab.io.info := pre.io.info

    // ----- post -----
    post.io.in.valid     := ntab.io.z.valid && polynomialEval.io.out.valid
    post.io.in.bits.info := ntab.io.z.bits
    post.io.in.bits.zres := polynomialEval.io.out.bits

    zs(Exp) := post.io.z
  }

  // --------------------------------------------------------------------------
  // Log -- uses the post-process multiplier (poly * const), no pre-mult.

  if(cfg.has(Log)) {
    assert(postMul.isDefined, "Log requires the post-process multiplier")

    val pre  = toModule[LogPreProcess ](pres   (Log))
    val tab  = toModule[LogTableCoeff ](tabs   (Log))
    val ntab = toModule[LogNonTable   ](nontabs(Log))
    val post = toModule[LogPostProcess](posts  (Log))

    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(Log) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Log) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    ntab.io.x := xTab

    // ----- post-mult arg builder (combinational, at xPost timing) -----
    val mulArgs = Module(new LogPostMultArgs(spec, bitW))
    mulArgs.io.info := ntab.io.z                  // at xPost timing
    mulArgs.io.zres := polynomialEval.io.out.bits // at xPost timing

    val logSelAtPost = ShiftRegister(selIn === cfg.signal(Log),
      stage.preStage + stage.tableStage + stage.calcStage, false.B, true.B)

    when (logSelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid && polynomialEval.io.out.valid
    }

    // ----- post -----
    post.io.info := Pipe(ntab.io.z, stage.postMulStage)
    post.io.mul  := postMul.get.io.out

    zs(Log) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ACosPhase1 -- reuses Sqrt's TableCoeff and NonTable.

  if(cfg.has(ACosPhase1)) {
    assert(cfg.has(Sqrt),
      "ACosPhase1 requires Sqrt (cfg.normalize should pull it in automatically)")

    val pre  = toModule[ACosPhase1PreProcess ](pres   (ACosPhase1))
    val post = toModule[ACosPhase1PostProcess](posts  (ACosPhase1))
    val sqrtTab  = toModule[SqrtTableCoeff](tabs   (Sqrt))
    val sqrtNtab = toModule[SqrtNonTable  ](nontabs(Sqrt))

    pre.io.x := xIn

    when (selTab === cfg.signal(ACosPhase1)) {
      sqrtTab.io.adr := pre.io.adr
      sqrtNtab.io.x  := pre.io.y
    }

    polynomialCoefs(ACosPhase1) := sqrtTab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ACosPhase1) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    // pre.io.flags is at xTab timing; pipeline to xPost = xTab + tableStage + calcStage.
    val flagsAtPost = Pipe(pre.io.flags, stage.tableStage + stage.calcStage)

    post.io.in.valid       := sqrtNtab.io.z.valid && polynomialEval.io.out.valid && flagsAtPost.valid
    post.io.in.bits.info   := sqrtNtab.io.z.bits
    post.io.in.bits.zres   := polynomialEval.io.out.bits
    post.io.in.bits.flags  := flagsAtPost.bits

    zs(ACosPhase1) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ACosPhase2 -- own table, uses the post-process multiplier (poly * x).

  if(cfg.has(ACosPhase2)) {
    assert(postMul.isDefined, "ACosPhase2 requires the post-process multiplier")

    val pre  = toModule[ACosPhase2PreProcess ](pres (ACosPhase2))
    val tab  = toModule[ACosPhase2TableCoeff ](tabs (ACosPhase2))
    val post = toModule[ACosPhase2PostProcess](posts(ACosPhase2))

    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(ACosPhase2) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ACosPhase2) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    // post-mult args at xPost timing: lhs = polynomial w/ hidden 1, rhs = x.manW1.
    val xManW1AtPost = Pipe(pre.io.xManW1, stage.tableStage + stage.calcStage)
    val mulArgs = Module(new ACosPhase2PostMultArgs(spec, bitW))
    mulArgs.io.zres   := polynomialEval.io.out.bits
    mulArgs.io.xManW1 := xManW1AtPost

    val acos2SelAtPost = ShiftRegister(selIn === cfg.signal(ACosPhase2),
      stage.preStage + stage.tableStage + stage.calcStage, false.B, true.B)

    when (acos2SelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid && polynomialEval.io.out.valid
    }

    // flags & x.ex pipelined to the post stage (after postMul)
    val flagsAtAfterMul = Pipe(pre.io.flags,
      stage.tableStage + stage.calcStage + stage.postMulStage)
    val xexAtAfterMul   = Pipe(pre.io.xex,
      stage.tableStage + stage.calcStage + stage.postMulStage)

    post.io.flags := flagsAtAfterMul
    post.io.xex   := xexAtAfterMul
    post.io.mul   := postMul.get.io.out

    zs(ACosPhase2) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ATan2Phase1 -- reuses Reciprocal's TableCoeff.
  //
  // Uses the post-process multiplier to assemble (1+min) * (2/(1+max)).

  if(cfg.has(ATan2Phase1)) {
    assert(cfg.has(Reciprocal),
      "ATan2Phase1 requires Reciprocal (cfg.normalize should pull it in automatically)")
    assert(postMul.isDefined, "ATan2Phase1 requires the post-process multiplier")
    assert(yIn.isDefined,     "ATan2Phase1 requires y as an additional input")

    val pre    = toModule[ATan2Phase1PreProcess ](pres (ATan2Phase1))
    val post   = toModule[ATan2Phase1PostProcess](posts(ATan2Phase1))
    val recTab = toModule[ReciprocalTableCoeff  ](tabs (Reciprocal))

    pre.io.x := xIn
    pre.io.y := yIn.get

    when (selTab === cfg.signal(ATan2Phase1)) {
      recTab.io.adr := pre.io.adr
    }

    polynomialCoefs(ATan2Phase1) := recTab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ATan2Phase1) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    // ---- post-mult args at xPost timing: lhs (poly+hidden 1), rhs (1+minMan) ----
    val infoAtPost = Pipe(pre.io.info, stage.tableStage + stage.calcStage)
    val mulArgs = Module(new ATan2Phase1PostMultArgs(spec, bitW))
    mulArgs.io.zres      := polynomialEval.io.out.bits
    mulArgs.io.maxXYMan0 := infoAtPost.bits.maxXYMan0
    mulArgs.io.minMan    := infoAtPost.bits.minMan
    mulArgs.io.valid     := infoAtPost.valid && polynomialEval.io.out.valid

    val atan2P1SelAtPost = ShiftRegister(selIn === cfg.signal(ATan2Phase1),
      stage.preStage + stage.tableStage + stage.calcStage, false.B, true.B)

    when (atan2P1SelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid
    }

    // info & flags pipelined to the post stage (after postMul)
    val infoAtAfterMul  = Pipe(pre.io.info,  stage.tableStage + stage.calcStage + stage.postMulStage)
    val flagsAtAfterMul = Pipe(pre.io.flags, stage.tableStage + stage.calcStage + stage.postMulStage)

    post.io.info  := infoAtAfterMul
    post.io.flags := flagsAtAfterMul
    post.io.mul   := postMul.get.io.out

    zs(ATan2Phase1) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ATan2Phase2 -- own table; uses the post-process multiplier (atan(t)/t * t).

  if(cfg.has(ATan2Phase2)) {
    assert(postMul.isDefined, "ATan2Phase2 requires the post-process multiplier")

    val pre  = toModule[ATan2Phase2PreProcess ](pres (ATan2Phase2))
    val tab  = toModule[ATan2Phase2TableCoeff ](tabs (ATan2Phase2))
    val post = toModule[ATan2Phase2PostProcess](posts(ATan2Phase2))

    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(ATan2Phase2) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ATan2Phase2) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    // ---- post-mult args at xPost timing ----
    val xManW1AtPost   = Pipe(pre.io.xManW1,   stage.tableStage + stage.calcStage)
    val isLinearAtPost = Pipe(pre.io.isLinear, stage.tableStage + stage.calcStage)
    val mulArgs = Module(new ATan2Phase2PostMultArgs(spec, bitW))
    mulArgs.io.zres     := polynomialEval.io.out.bits
    mulArgs.io.isLinear := isLinearAtPost.bits
    mulArgs.io.xManW1   := xManW1AtPost

    val atan2P2SelAtPost = ShiftRegister(selIn === cfg.signal(ATan2Phase2),
      stage.preStage + stage.tableStage + stage.calcStage, false.B, true.B)

    when (atan2P2SelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid && polynomialEval.io.out.valid
    }

    // flags, xex, isLinear pipelined to the post stage (after postMul)
    val flagsAtAfterMul    = Pipe(pre.io.flags,    stage.tableStage + stage.calcStage + stage.postMulStage)
    val xexAtAfterMul      = Pipe(pre.io.xex,      stage.tableStage + stage.calcStage + stage.postMulStage)
    val isLinearAtAfterMul = Pipe(pre.io.isLinear, stage.tableStage + stage.calcStage + stage.postMulStage)

    post.io.flags    := flagsAtAfterMul
    post.io.xex      := xexAtAfterMul
    post.io.isLinear := isLinearAtAfterMul
    post.io.mul      := postMul.get.io.out

    zs(ATan2Phase2) := post.io.z
  }

  // --------------------------------------------------------------------------
  // Sigmoid -- own table, no pre/post-mult (the post-process fold is a shift).

  if(cfg.has(Sigmoid)) {
    val pre  = toModule[SigmoidPreProcess ](pres (Sigmoid))
    val tab  = toModule[SigmoidTableCoeff ](tabs (Sigmoid))
    val post = toModule[SigmoidPostProcess](posts(Sigmoid))

    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(Sigmoid) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(Sigmoid) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    // info pipelined from xTab to xPost (across tableStage + calcStage)
    val infoAtPost = Pipe(pre.io.info, stage.tableStage + stage.calcStage)
    post.io.in.valid     := infoAtPost.valid && polynomialEval.io.out.valid
    post.io.in.bits.info := infoAtPost.bits
    post.io.in.bits.zres := polynomialEval.io.out.bits

    zs(Sigmoid) := post.io.z
  }

  // --------------------------------------------------------------------------
  // SoftPlus -- own table (concat pos+neg, Mux'd by xsgn), no pre/post-mult.

  if(cfg.has(SoftPlus)) {
    val pre  = toModule[SoftPlusPreProcess ](pres (SoftPlus))
    val tab  = toModule[SoftPlusTableCoeff ](tabs (SoftPlus))
    val post = toModule[SoftPlusPostProcess](posts(SoftPlus))

    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(SoftPlus) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(SoftPlus) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    val infoAtPost = Pipe(pre.io.info, stage.tableStage + stage.calcStage)
    post.io.in.valid     := infoAtPost.valid && polynomialEval.io.out.valid
    post.io.in.bits.info := infoAtPost.bits
    post.io.in.bits.zres := polynomialEval.io.out.bits

    zs(SoftPlus) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ScaleMixtureGaussian -- uses BOTH the pre-process and post-process
  // multipliers.
  //   pre-mult: (1/sgmA^2).manW1 * x.manW1   (= x/sgmA^2)
  //   post-mult: z0.manW1 * xsgmA2.manW1    (where z0 = polynomial*2^maxBitDigit + 1)

  if(cfg.has(ScaleMixtureGaussian)) {
    assert(preMul .isDefined, "ScaleMixtureGaussian requires the pre-process multiplier")
    assert(postMul.isDefined, "ScaleMixtureGaussian requires the post-process multiplier")
    val (sgmA, sgmB) = cfg.scaleMixtureGaussianSigma.get

    val pre  = toModule[ScaleMixtureGaussianPreProcess ](pres (ScaleMixtureGaussian))
    val tab  = toModule[ScaleMixtureGaussianTableCoeff ](tabs (ScaleMixtureGaussian))
    val post = toModule[ScaleMixtureGaussianPostProcess](posts(ScaleMixtureGaussian))

    val smgSelIn = selIn === cfg.signal(ScaleMixtureGaussian)

    // ----- pre-mult driver: (1/sgmA^2) * (1 + x.man) -----
    when (smgSelIn) {
      preMulLhs.get   := ScaleMixtureGaussianPreMulArgs.lhs(sgmA, spec, preMul.get.lhsW)
      preMulRhs.get   := Cat(1.U(1.W), xIn.bits.man)
      preMulValid.get := xIn.valid
    }

    // ----- pre stage: x delayed by preMulStage to align with preMul.out -----
    val xAtSmgPre = Pipe(xIn, stage.preMulStage)
    pre.io.x          := xAtSmgPre
    pre.io.prod.valid := preMul.get.io.out.valid
    pre.io.prod.bits  := ScaleMixtureGaussianPreMulArgs.prod(spec, preMul.get.io.out.bits)

    // ----- table -----
    tab.io.adr := pre.io.adr
    polynomialCoefs(ScaleMixtureGaussian) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ScaleMixtureGaussian) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    // ----- post-mult args at xPost timing -----
    val infoAtPost = Pipe(pre.io.info, stage.tableStage + stage.calcStage)
    val mulArgs    = Module(new ScaleMixtureGaussianPostMulArgs(sgmA, sgmB, spec, bitW))
    mulArgs.io.useTable  := infoAtPost.bits.useTable
    mulArgs.io.xsgmA2Man := infoAtPost.bits.xsgmA2Man
    mulArgs.io.zres      := polynomialEval.io.out.bits
    mulArgs.io.valid     := infoAtPost.valid && polynomialEval.io.out.valid

    val smgSelAtPost = ShiftRegister(smgSelIn,
      stage.preStage + stage.tableStage + stage.calcStage, false.B, true.B)
    when (smgSelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid
    }

    // ----- post stage: pipe info & z0Ex to post-postMul timing -----
    val infoAtAfterMul = Pipe(infoAtPost, stage.postMulStage)
    val z0exAtAfterMul = Pipe(mulArgs.io.mul.valid, mulArgs.io.mul.bits.z0ex,
                              stage.postMulStage)

    post.io.info := infoAtAfterMul
    post.io.z0ex := z0exAtAfterMul
    post.io.mul  := postMul.get.io.out
    zs(ScaleMixtureGaussian) := post.io.z
  }
}
