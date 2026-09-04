//% @file MathFunctions.scala
//
// Math functions module
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{RealSpec, RoundSpec}
import rial.util._

import scala.collection.mutable.ArrayBuffer
import scala.reflect.ClassTag

/** Evaluates one of a configured set of math functions on a shared datapath.
 *
 * The function to evaluate is selected every cycle by `io.in.bits.sel`, a
 * signal built with [[rial.math.MathFuncConfig.signal]]; `io.in.bits.x` is its
 * argument and `io.in.bits.y`, present only when
 * [[rial.math.FuncKind.ATan2Phase1]] is configured, its second one. The result
 * of the selected function appears on `io.out`.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`        the module input
 *  - `"preMul"`    after the shared pre-process product
 *  - `"pre"`       after the argument reduction
 *  - `"table"`     after the coefficient table lookup
 *  - `"poly.*"`    the boundaries of the [[rial.table.PolynomialEval]]
 *                  submodule that evaluates the polynomial of the interval,
 *                  see its doc for the names; the group name `"poly"` stands
 *                  for `"poly.out"`
 *  - `"postMul.*"` the boundaries of the shared post-process multiplier,
 *                  `"postMul.prod"` after its product and `"postMul.out"`
 *                  after its rounding; `"postMul"` stands for `"postMul.out"`
 *  - `"post"`      after the result assembly
 *  - `"out"`       the module output
 *
 * `"preMul"` and `"postMul.*"` exist only when one of the configured functions
 * needs the corresponding multiplier; [[MathFunctions.stagePoints]] gives the
 * points of a given configuration.
 *
 * [[getStage]] gives the resulting latency, and `io.out.valid` is
 * `io.in.valid` with a `sel` that names a configured function, delayed by that
 * many cycles.
 *
 * {{{
 * // six registers, placed by the module
 * new MathFunctions(cfg, spec, 6)
 * // hand-picked positions, one of them inside the polynomial
 * new MathFunctions(cfg, spec,
 *   PipelineStageConfig.specify("pre" -> 1, "poly.prod1" -> 1, "post" -> 1))
 * // combinational
 * new MathFunctions(cfg, spec, PipelineStageConfig.none)
 * }}}
 *
 * @constructor Constructs a new MathFunctions module.
 * @param cfg   Configuration that determines which functions are supported.
 * @param spec  RealSpec for input/output floating point numbers.
 * @param stage pipeline registers of this module.
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class MathFunctions(
  val cfg:   MathFuncConfig,
  val spec:  RealSpec,
  val stage: PipelineStageConfig
) extends Module {

  import FuncKind._

  assert(!spec.disableSign, "assuming the float type allows negative values")

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val bitW   = cfg.bitW
  val coeffs = cfg.coeffs

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = MathFunctions.stagePoints(cfg)

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new MathFunctions.Input(cfg, spec)))
    val out = Valid(UInt(spec.W.W))
  })

  /** pipes a `Valid` through one cut of the datapath. */
  private def cutValid[T <: Data](cut: SharedPipe, x: Valid[T]): Valid[T] = {
    val w = Wire(Valid(chiselTypeOf(x.bits)))
    w.valid := cut.valid
    w.bits  := cut(x.bits)
    w
  }

  // ==========================================================================
  // stage "in": the input, decomposed into its floating point fields
  // ==========================================================================

  private val cutIn = new SharedPipe("in",
    io.in.valid && cfg.isValid(io.in.bits.sel), stages.at("in"))

  val selIn = cutIn(io.in.bits.sel)

  val xdecomp = Module(new DecomposeReal(spec))
  xdecomp.io.real.valid := cutIn.valid
  xdecomp.io.real.bits  := cutIn(io.in.bits.x)

  val xIn: Valid[DecomposedReal] = xdecomp.io.decomp

  // ATan2 also needs the y operand decomposed.
  val ydecomp: Option[DecomposeReal] = if (cfg.has(ATan2Phase1)) {
    val m = Module(new DecomposeReal(spec))
    m.io.real.valid := cutIn.valid
    m.io.real.bits  := cutIn(io.in.bits.y.get)
    Some(m)
  } else {
    None
  }

  val yIn: Option[Valid[DecomposedReal]] = ydecomp.map(m => m.io.decomp)

  // ==========================================================================
  // Shared pre/post-process multipliers (only built if some configured
  // function needs them; sized to the minimum width that satisfies all
  // such functions).
  // ==========================================================================

  val needPreMul = PreProcMultiplier.isNeeded(cfg.funcs)

  val nPreMul = if (needPreMul) {
    stages.at("preMul")
  } else {
    0
  }

  val preMul: Option[PreProcMultiplier] = if (needPreMul) { Some(
    Module(new PreProcMultiplier(
      lhsW  = PreProcMultiplier.lhsW(cfg.funcs, spec),
      manW  = manW,
      stage = nPreMul,
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

  // ==========================================================================
  // stage "preMul": everything that bypasses the pre-process multiplier
  // ==========================================================================

  private val cutPreMul = new SharedPipe("preMul", cutIn.valid, nPreMul)

  val selPre = cutPreMul(selIn)
  val xPre   = cutValid(cutPreMul, xIn)
  val yPre   = yIn.map(y => cutValid(cutPreMul, y))

  // --------------------------------------------------------------------------

  val needPostMul = PostProcMultiplier.isNeeded(cfg.funcs)

  val nPostMul = if (needPostMul) {
    stages.sub("postMul").total
  } else {
    0
  }

  val postMul: Option[PostProcMultiplier] = if (needPostMul) { Some(
    Module(new PostProcMultiplier(
      fracW     = bitW.outBitW,
      manW      = manW,
      roundSpec = RoundSpec.roundToEven,
      stage     = stages.sub("postMul"),
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
  // stage "pre": the argument reduction of every function
  // ==========================================================================

  val nPre = stages.at("pre")

  private val cutPre = new SharedPipe("pre", cutPreMul.valid, nPre)

  val selTab = cutPre(selPre)
  val xTab   = cutValid(cutPre, xPre)

  // ==========================================================================
  // stage "table": the coefficient table lookup of every function
  // ==========================================================================

  val nTable = stages.at("table")

  private val cutTable = new SharedPipe("table", cutPre.valid, nTable)

  val selCalc = cutTable(selTab)

  // ==========================================================================
  // Polynomial Evaluator (shared), stages "poly.*"
  // ==========================================================================

  val polynomialEval = Module(new PolynomialEval(bitW, stages.sub("poly")))

  val nPoly = stages.sub("poly").total

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
  // stage "poly": everything that bypasses the polynomial
  // ==========================================================================

  private val cutPoly = new SharedPipe("poly", cutTable.valid, nPoly)

  val selPost = cutPoly(selCalc)

  // ==========================================================================
  // stage "postMul": everything that bypasses the post-process multiplier
  // ==========================================================================

  private val cutPostMul = new SharedPipe("postMul", cutPoly.valid, nPostMul)

  val selPostMul = cutPostMul(selPost)

  // the polynomial result as the functions that do not use the post-process
  // multiplier see it, so that every post-process sits at the same boundary
  private val cutPostMulRes = new SharedPipe("postMulRes",
    polynomialEval.io.out.valid, nPostMul)

  val zresPostMul = cutValid(cutPostMulRes, polynomialEval.io.out)

  // ==========================================================================
  // stage "post": the result assembly of every function
  // ==========================================================================

  val nPost = stages.at("post")

  private val cutPost = new SharedPipe("post", cutPostMul.valid, nPost)

  val selOut = cutPost(selPostMul)

  // ==========================================================================
  // Output selection, stage "out"
  // ==========================================================================

  val zs = cfg.funcs.map(fn => {
    fn -> WireDefault(0.U.asTypeOf(Valid(UInt(spec.W.W))))
  }).toMap

  val zMux = MuxCase(0.U, cfg.funcs.map(fn => {
    (selOut === cfg.signal(fn) && zs(fn).valid) -> zs(fn).bits
  }))

  io.out := Pipe(cutPost.valid, zMux, stages.at("out"))

  // ==========================================================================
  // Per-function pre/tab/nontab/post modules
  //
  // `.collect` so that some funcs does not have the corresponding modules
  // ==========================================================================

  val pres: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtPreProcess       (spec, bitW, nPre))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtPreProcess    (spec, bitW, nPre))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalPreProcess (spec, bitW, nPre))
    case Sin         => Sin         -> Module(new SinCosPreProcess     (spec, bitW, nPre))
    case Exp         => Exp         -> Module(new ExpPreProcess        (spec, bitW, nPre))
    case Log         => Log         -> Module(new LogPreProcess        (spec, bitW, nPre))
    case ACosPhase1  => ACosPhase1  -> Module(new ACosPhase1PreProcess (spec, bitW, nPre))
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2PreProcess (spec, bitW, nPre))
    case ATan2Phase1 => ATan2Phase1 -> Module(new ATan2Phase1PreProcess(spec, bitW, nPre))
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2PreProcess(spec, bitW, nPre))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidPreProcess    (spec, bitW, nPre))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusPreProcess   (spec, bitW, nPre))
    case ScaleMixtureGaussian => {
      val (sA, sB) = cfg.scaleMixtureGaussianSigma.get
      ScaleMixtureGaussian -> Module(new ScaleMixtureGaussianPreProcess(
        sA, sB, spec, bitW, nPre))
    }
  }.toMap[FuncKind, Module]

  val tabs: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtTableCoeff       (spec, coeffs(Sqrt       ), nTable))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtTableCoeff    (spec, coeffs(InvSqrt    ), nTable))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalTableCoeff (spec, coeffs(Reciprocal ), nTable))
    case Sin         => Sin         -> Module(new SinCosTableCoeff     (spec, coeffs(Sin        ), nTable))
    case Log         => Log         -> Module(new LogTableCoeff        (spec, coeffs(Log        ), nTable))
    case Exp         => Exp         -> Module(new ExpTableCoeff        (spec, coeffs(Exp        ), nTable))
    // ACosPhase1 shares Sqrt's TableCoeff
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2TableCoeff (spec, coeffs(ACosPhase2 ), nTable))
    // ATan2Phase1 shares Reciprocal's TableCoeff
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2TableCoeff(spec, coeffs(ATan2Phase2), nTable))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidTableCoeff    (spec, coeffs(Sigmoid    ), nTable))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusTableCoeff   (spec, coeffs(SoftPlus   ), nTable))
    case ScaleMixtureGaussian => ScaleMixtureGaussian ->
      Module(new ScaleMixtureGaussianTableCoeff(spec, coeffs(ScaleMixtureGaussian), nTable))
  }.toMap

  val nontabs: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt       => Sqrt       -> Module(new SqrtNonTable      (spec, bitW, nTable + nPoly))
    case InvSqrt    => InvSqrt    -> Module(new InvSqrtNonTable   (spec, bitW, nTable + nPoly))
    case Reciprocal => Reciprocal -> Module(new ReciprocalNonTable(spec, bitW, nTable + nPoly))
    // sincos / acos do not have a dedicated non-table path
    case Exp        => Exp        -> Module(new ExpNonTable       (spec, bitW, nTable + nPoly))
    case Log        => Log        -> Module(new LogNonTable       (spec, bitW, nTable + nPoly))
    // ACosPhase1 reuses SqrtNonTable
    // ACosPhase2 doesn't need one.
  }.toMap

  val posts: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtPostProcess       (spec, bitW, nPost))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtPostProcess    (spec, bitW, nPost))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalPostProcess (spec, bitW, nPost))
    case Sin         => Sin         -> Module(new SinCosPostProcess     (spec, bitW, nPost))
    case Exp         => Exp         -> Module(new ExpPostProcess        (spec, bitW, nPost))
    case Log         => Log         -> Module(new LogPostProcess        (spec, bitW, nPost))
    case ACosPhase1  => ACosPhase1  -> Module(new ACosPhase1PostProcess (spec, bitW, nPost))
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2PostProcess (spec, bitW, nPost))
    case ATan2Phase1 => ATan2Phase1 -> Module(new ATan2Phase1PostProcess(spec, bitW, nPost))
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2PostProcess(spec, bitW, nPost))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidPostProcess    (spec, bitW, nPost))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusPostProcess   (spec, bitW, nPost))
    case ScaleMixtureGaussian => ScaleMixtureGaussian ->
      Module(new ScaleMixtureGaussianPostProcess(spec, bitW, nPost))
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
    pre.io.x.valid := xPre.valid && selPre === cfg.signal(Sqrt)
    pre.io.x.bits  := xPre.bits
    tab.io.adr := pre.io.adr

    polynomialCoefs(Sqrt) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Sqrt) := Pipe(pre.io.dx.get, nTable)
    }

    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(Sqrt)
    ntab.io.x.bits  := xTab.bits

    val ntabAtPost = Pipe(ntab.io.z, nPostMul)

    post.io.in.valid     := ntabAtPost.valid && zresPostMul.valid
    post.io.in.bits.info := ntabAtPost.bits
    post.io.in.bits.zres := zresPostMul.bits

    zs(Sqrt) := post.io.z
  }

  // --------------------------------------------------------------------------

  if(cfg.has(InvSqrt)) {
    val pre  = toModule[InvSqrtPreProcess ](pres   (InvSqrt))
    val tab  = toModule[InvSqrtTableCoeff ](tabs   (InvSqrt))
    val ntab = toModule[InvSqrtNonTable   ](nontabs(InvSqrt))
    val post = toModule[InvSqrtPostProcess](posts  (InvSqrt))
    pre.io.x.valid := xPre.valid && selPre === cfg.signal(InvSqrt)
    pre.io.x.bits  := xPre.bits
    tab.io.adr := pre.io.adr

    polynomialCoefs(InvSqrt) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(InvSqrt) := Pipe(pre.io.dx.get, nTable)
    }

    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(InvSqrt)
    ntab.io.x.bits  := xTab.bits

    val ntabAtPost = Pipe(ntab.io.z, nPostMul)

    post.io.in.valid     := ntabAtPost.valid && zresPostMul.valid
    post.io.in.bits.info := ntabAtPost.bits
    post.io.in.bits.zres := zresPostMul.bits

    zs(InvSqrt) := post.io.z
  }

  // --------------------------------------------------------------------------

  if(cfg.has(Reciprocal)) {
    val pre  = toModule[ReciprocalPreProcess ](pres   (Reciprocal))
    val tab  = toModule[ReciprocalTableCoeff ](tabs   (Reciprocal))
    val ntab = toModule[ReciprocalNonTable   ](nontabs(Reciprocal))
    val post = toModule[ReciprocalPostProcess](posts  (Reciprocal))
    pre.io.x.valid := xPre.valid && selPre === cfg.signal(Reciprocal)
    pre.io.x.bits  := xPre.bits
    tab.io.adr := pre.io.adr

    polynomialCoefs(Reciprocal) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Reciprocal) := Pipe(pre.io.dx.get, nTable)
    }

    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(Reciprocal)
    ntab.io.x.bits  := xTab.bits

    val ntabAtPost = Pipe(ntab.io.z, nPostMul)

    post.io.in.valid     := ntabAtPost.valid && zresPostMul.valid
    post.io.in.bits.info := ntabAtPost.bits
    post.io.in.bits.zres := zresPostMul.bits

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

    // ----- pre stage: x and sel on the far side of the pre-mult -----
    val isSinAtPre     = selPre === cfg.signal(Sin)
    val sincosSelAtPre = isSinAtPre || (selPre === cfg.signal(Cos))

    sincosPre.io.isSin := isSinAtPre
    sincosPre.io.x.valid := xPre.valid && sincosSelAtPre
    sincosPre.io.x.bits  := xPre.bits

    sincosPre.io.prod.valid := preMul.get.io.out.valid && sincosSelAtPre
    sincosPre.io.prod.bits  := SinCosPreMulArgs.prod(spec, preMul.get.io.out.bits)

    // ----- table -----
    sincosTab.io.adr := sincosPre.io.adr

    polynomialCoefs(Sin) := sincosTab.io.cs
    polynomialCoefs(Cos) := sincosTab.io.cs
    if(bitW.nOrder != 0) {
      val dxAtCalc = Pipe(sincosPre.io.dx.get, nTable)
      polynomialDxs.get(Sin) := dxAtCalc
      polynomialDxs.get(Cos) := dxAtCalc
    }

    // ----- post-mult driver: polynomial * (1 + yman) -----
    val sincosInfoAtPost = Pipe(sincosPre.io.info, nTable + nPoly)
    val sincosSelAtPost  = (selPost === cfg.signal(Sin)) ||
      (selPost === cfg.signal(Cos))

    when (sincosSelAtPost) {
      postMulLhs.get := Cat(1.U(1.W), polynomialEval.io.out.bits)
      postMulRhs.get := Cat(1.U(1.W), sincosInfoAtPost.bits.yman)
      postMulValid.get := polynomialEval.io.out.valid && sincosInfoAtPost.valid
    }

    // ----- post -----
    sincosPost.io.pre := Pipe(sincosInfoAtPost, nPostMul)
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
    val expSelAtPre = selPre === cfg.signal(Exp)

    pre.io.x.valid := xPre.valid && expSelAtPre
    pre.io.x.bits  := xPre.bits
    pre.io.prod.valid := preMul.get.io.out.valid && expSelAtPre
    pre.io.prod.bits  := ExpPreMulArgs.prod(spec, preMul.get.io.out.bits)

    // ----- table -----
    tab.io.adr := pre.io.adr

    polynomialCoefs(Exp) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Exp) := Pipe(pre.io.dx.get, nTable)
    }

    // ----- non-table ------
    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(Exp)
    ntab.io.x.bits  := xTab.bits
    ntab.io.info := pre.io.info

    // ----- post -----
    val ntabAtPost = Pipe(ntab.io.z, nPostMul)

    post.io.in.valid     := ntabAtPost.valid && zresPostMul.valid
    post.io.in.bits.info := ntabAtPost.bits
    post.io.in.bits.zres := zresPostMul.bits

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

    pre.io.x.valid := xPre.valid && selPre === cfg.signal(Log)
    pre.io.x.bits  := xPre.bits
    tab.io.adr := pre.io.adr

    polynomialCoefs(Log) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Log) := Pipe(pre.io.dx.get, nTable)
    }

    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(Log)
    ntab.io.x.bits  := xTab.bits

    // ----- post-mult arg builder (combinational, at the polynomial output) -----
    val mulArgs = Module(new LogPostMultArgs(spec, bitW))
    mulArgs.io.info := ntab.io.z
    mulArgs.io.zres := polynomialEval.io.out.bits

    val logSelAtPost = selPost === cfg.signal(Log)

    when (logSelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid && polynomialEval.io.out.valid
    }

    // ----- post -----
    post.io.info := Pipe(ntab.io.z, nPostMul)
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

    pre.io.x.valid := xPre.valid && selPre === cfg.signal(ACosPhase1)
    pre.io.x.bits  := xPre.bits

    when (selTab === cfg.signal(ACosPhase1)) {
      sqrtTab.io.adr := pre.io.adr
      sqrtNtab.io.x  := pre.io.y
    }

    polynomialCoefs(ACosPhase1) := sqrtTab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ACosPhase1) := Pipe(pre.io.dx.get, nTable)
    }

    // pre.io.flags leaves the pre-process at the table input; pipeline it to
    // the input of the post-process
    val flagsAtPost = Pipe(pre.io.flags, nTable + nPoly + nPostMul)
    val ntabAtPost  = Pipe(sqrtNtab.io.z, nPostMul)

    post.io.in.valid       := ntabAtPost.valid && zresPostMul.valid && flagsAtPost.valid
    post.io.in.bits.info   := ntabAtPost.bits
    post.io.in.bits.zres   := zresPostMul.bits
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

    pre.io.x.valid := xPre.valid && selPre === cfg.signal(ACosPhase2)
    pre.io.x.bits  := xPre.bits
    tab.io.adr := pre.io.adr

    polynomialCoefs(ACosPhase2) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ACosPhase2) := Pipe(pre.io.dx.get, nTable)
    }

    // post-mult args at the polynomial output: lhs = polynomial w/ hidden 1,
    // rhs = x.manW1.
    val xManW1AtPost = Pipe(pre.io.xManW1, nTable + nPoly)
    val mulArgs = Module(new ACosPhase2PostMultArgs(spec, bitW))
    mulArgs.io.zres   := polynomialEval.io.out.bits
    mulArgs.io.xManW1 := xManW1AtPost

    val acos2SelAtPost = selPost === cfg.signal(ACosPhase2)

    when (acos2SelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid && polynomialEval.io.out.valid
    }

    // flags & x.ex pipelined to the input of the post-process
    val flagsAtAfterMul = Pipe(pre.io.flags, nTable + nPoly + nPostMul)
    val xexAtAfterMul   = Pipe(pre.io.xex,   nTable + nPoly + nPostMul)

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

    pre.io.x.valid := xPre.valid && selPre === cfg.signal(ATan2Phase1)
    pre.io.x.bits  := xPre.bits
    pre.io.y.valid := yPre.get.valid && selPre === cfg.signal(ATan2Phase1)
    pre.io.y.bits  := yPre.get.bits

    when (selTab === cfg.signal(ATan2Phase1)) {
      recTab.io.adr := pre.io.adr
    }

    polynomialCoefs(ATan2Phase1) := recTab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ATan2Phase1) := Pipe(pre.io.dx.get, nTable)
    }

    // ---- post-mult args at the polynomial output: lhs (poly+hidden 1), rhs (1+minMan) ----
    val infoAtPost = Pipe(pre.io.info, nTable + nPoly)
    val mulArgs = Module(new ATan2Phase1PostMultArgs(spec, bitW))
    mulArgs.io.zres      := polynomialEval.io.out.bits
    mulArgs.io.maxXYMan0 := infoAtPost.bits.maxXYMan0
    mulArgs.io.minMan    := infoAtPost.bits.minMan
    mulArgs.io.valid     := infoAtPost.valid && polynomialEval.io.out.valid

    val atan2P1SelAtPost = selPost === cfg.signal(ATan2Phase1)

    when (atan2P1SelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid
    }

    // info & flags pipelined to the input of the post-process
    val infoAtAfterMul  = Pipe(pre.io.info,  nTable + nPoly + nPostMul)
    val flagsAtAfterMul = Pipe(pre.io.flags, nTable + nPoly + nPostMul)

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

    pre.io.x.valid := xPre.valid && selPre === cfg.signal(ATan2Phase2)
    pre.io.x.bits  := xPre.bits
    tab.io.adr := pre.io.adr

    polynomialCoefs(ATan2Phase2) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ATan2Phase2) := Pipe(pre.io.dx.get, nTable)
    }

    // ---- post-mult args at the polynomial output ----
    val xManW1AtPost   = Pipe(pre.io.xManW1,   nTable + nPoly)
    val isLinearAtPost = Pipe(pre.io.isLinear, nTable + nPoly)
    val mulArgs = Module(new ATan2Phase2PostMultArgs(spec, bitW))
    mulArgs.io.zres     := polynomialEval.io.out.bits
    mulArgs.io.isLinear := isLinearAtPost.bits
    mulArgs.io.xManW1   := xManW1AtPost

    val atan2P2SelAtPost = selPost === cfg.signal(ATan2Phase2)

    when (atan2P2SelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid && polynomialEval.io.out.valid
    }

    // flags, xex, isLinear pipelined to the input of the post-process
    val flagsAtAfterMul    = Pipe(pre.io.flags,    nTable + nPoly + nPostMul)
    val xexAtAfterMul      = Pipe(pre.io.xex,      nTable + nPoly + nPostMul)
    val isLinearAtAfterMul = Pipe(pre.io.isLinear, nTable + nPoly + nPostMul)

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

    pre.io.x.valid := xPre.valid && selPre === cfg.signal(Sigmoid)
    pre.io.x.bits  := xPre.bits
    tab.io.adr := pre.io.adr

    polynomialCoefs(Sigmoid) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(Sigmoid) := Pipe(pre.io.dx.get, nTable)
    }

    // info pipelined from the table input to the input of the post-process
    val infoAtPost = Pipe(pre.io.info, nTable + nPoly + nPostMul)
    post.io.in.valid     := infoAtPost.valid && zresPostMul.valid
    post.io.in.bits.info := infoAtPost.bits
    post.io.in.bits.zres := zresPostMul.bits

    zs(Sigmoid) := post.io.z
  }

  // --------------------------------------------------------------------------
  // SoftPlus -- own table (concat pos+neg, Mux'd by xsgn), no pre/post-mult.

  if(cfg.has(SoftPlus)) {
    val pre  = toModule[SoftPlusPreProcess ](pres (SoftPlus))
    val tab  = toModule[SoftPlusTableCoeff ](tabs (SoftPlus))
    val post = toModule[SoftPlusPostProcess](posts(SoftPlus))

    pre.io.x.valid := xPre.valid && selPre === cfg.signal(SoftPlus)
    pre.io.x.bits  := xPre.bits
    tab.io.adr := pre.io.adr

    polynomialCoefs(SoftPlus) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(SoftPlus) := Pipe(pre.io.dx.get, nTable)
    }

    val infoAtPost = Pipe(pre.io.info, nTable + nPoly + nPostMul)
    post.io.in.valid     := infoAtPost.valid && zresPostMul.valid
    post.io.in.bits.info := infoAtPost.bits
    post.io.in.bits.zres := zresPostMul.bits

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

    // ----- pre stage: x and sel on the far side of the pre-mult -----
    val smgSelAtPre = selPre === cfg.signal(ScaleMixtureGaussian)
    pre.io.x.valid := xPre.valid && smgSelAtPre
    pre.io.x.bits  := xPre.bits
    pre.io.prod.valid := preMul.get.io.out.valid && smgSelAtPre
    pre.io.prod.bits  := ScaleMixtureGaussianPreMulArgs.prod(spec, preMul.get.io.out.bits)

    // ----- table -----
    tab.io.adr := pre.io.adr
    polynomialCoefs(ScaleMixtureGaussian) := tab.io.cs
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ScaleMixtureGaussian) := Pipe(pre.io.dx.get, nTable)
    }

    // ----- post-mult args at the polynomial output -----
    val infoAtPost = Pipe(pre.io.info, nTable + nPoly)
    val mulArgs    = Module(new ScaleMixtureGaussianPostMulArgs(sgmA, sgmB, spec, bitW))
    mulArgs.io.useTable  := infoAtPost.bits.useTable
    mulArgs.io.xsgmA2Man := infoAtPost.bits.xsgmA2Man
    mulArgs.io.zres      := polynomialEval.io.out.bits
    mulArgs.io.valid     := infoAtPost.valid && polynomialEval.io.out.valid

    val smgSelAtPost = selPost === cfg.signal(ScaleMixtureGaussian)
    when (smgSelAtPost) {
      postMulLhs.get   := mulArgs.io.mul.bits.lhs
      postMulRhs.get   := mulArgs.io.mul.bits.rhs
      postMulValid.get := mulArgs.io.mul.valid
    }

    // ----- post stage: pipe info & z0Ex to the input of the post-process -----
    val infoAtAfterMul = Pipe(infoAtPost, nPostMul)
    val z0exAtAfterMul = Pipe(mulArgs.io.mul.valid, mulArgs.io.mul.bits.z0ex,
                              nPostMul)

    post.io.info := infoAtAfterMul
    post.io.z0ex := z0exAtAfterMul
    post.io.mul  := postMul.get.io.out
    zs(ScaleMixtureGaussian) := post.io.z
  }
}

/** Input bundle and pipeline stage points of [[rial.math.MathFunctions]].
 *
 * @see [[rial.util.StagePoint]]
 * @see [[rial.util.StageGroup]]
 */
object MathFunctions {

  /** Input of [[rial.math.MathFunctions]]: the selected function and its
   * argument(s).
   *
   * `sel` is a signal built with [[rial.math.MathFuncConfig.signal]]; `y` is
   * present only when [[rial.math.FuncKind.ATan2Phase1]] is configured.
   *
   * @param cfg  configuration of the module this input is for
   * @param spec floating point spec of the arguments
   */
  class Input(val cfg: MathFuncConfig, val spec: RealSpec) extends Bundle {
    val sel = UInt(cfg.signalW.W)
    val x   = UInt(spec.W.W)
    val y   = if (cfg.has(FuncKind.ATan2Phase1)) {
      Some(UInt(spec.W.W))
    } else {
      None
    }
  }

  /** Boundaries at which pipeline registers can be placed; see
   * [[rial.math.MathFunctions]]. The boundaries of the polynomial and of the
   * post-process multiplier appear as the groups `"poly"` and `"postMul"`.
   *
   * `"preMul"` and the group `"postMul"` are declared only when one of the
   * functions of `cfg` needs the corresponding multiplier.
   *
   * @param cfg configuration of the module these boundaries are for
   */
  def stagePoints(cfg: MathFuncConfig): Seq[StageElem] = {
    val preMul: Seq[StageElem] = if (PreProcMultiplier.isNeeded(cfg.funcs)) {
      Seq(StagePoint("preMul", 3.0))
    } else {
      Seq.empty[StageElem]
    }
    val postMul: Seq[StageElem] = if (PostProcMultiplier.isNeeded(cfg.funcs)) {
      Seq(StageGroup("postMul", PostProcMultiplier.stagePoints))
    } else {
      Seq.empty[StageElem]
    }
    preMul ++
    Seq(StagePoint("pre",   4.0),
        StagePoint("table", 2.0),
        StageGroup("poly",  PolynomialEval.stagePoints(cfg.bitW.nOrder))) ++
    postMul ++
    Seq(StagePoint("post", 3.0),
        StagePoint("out",  1.0))
  }
}
