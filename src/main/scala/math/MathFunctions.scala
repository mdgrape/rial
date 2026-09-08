//% @file MathFunctions.scala
//
// Math functions module
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{MulImpl, Multiplier, RealSpec, RoundSpec}
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
 *  - `"preMul.*"`  the boundaries of the shared pre-process multiplier, see
 *                  [[rial.arith.Multiplier]] for the names
 *  - `"pre.*"`     the argument reduction of each function, one path per
 *                  function, see below
 *  - `"table"`     after the coefficient table lookup
 *  - `"coef"`      after the coefficient selection
 *  - `"poly.*"`    the boundaries of the [[rial.table.PolynomialEval]]
 *                  submodule that evaluates the polynomial of the interval,
 *                  see its doc for the names; the group name `"poly"` stands
 *                  for `"poly.out"`
 *  - `"mulArgs"`   after the post-process multiplier operand selection
 *  - `"postMul.*"` the boundaries of the shared post-process multiplier,
 *                  see [[rial.math.PostProcMultiplier]] for the names
 *  - `"post.*"`    the result assembly of each function, one path per
 *                  function, see below
 *  - `"out"`       the module output
 *
 * `"pre"` and `"post"` are sets of parallel paths, one per configured
 * function, named after it with [[rial.math.FuncKind.pointName]]: the
 * boundaries of the argument reduction of `Sqrt` are `"pre.sqrt.*"` and those
 * of its result assembly `"post.sqrt.*"`, see the doc of the corresponding
 * submodule for the names. `Sin` and `Cos` share the path `"sinCos"`. Every
 * path of a set has the same latency, so `"pre" -> n` (or `"post" -> n`) adds
 * `n` registers to every one of them.
 *
 * `"preMul.*"`, `"mulArgs"` and `"postMul.*"` exist only when one of the
 * configured functions needs the corresponding multiplier;
 * [[MathFunctions.stagePoints]] gives the points of a given configuration.
 *
 * The modules that compute the exponent and the special values of a function
 * beside the table and the polynomial place their own registers, so they have
 * no user-visible point.
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
 * // one register inside the argument reduction of sin/cos
 * new MathFunctions(cfg, spec, PipelineStageConfig.specify("pre.sinCos.fold" -> 1))
 * // combinational
 * new MathFunctions(cfg, spec, PipelineStageConfig.none)
 * }}}
 *
 * @constructor Constructs a new MathFunctions module.
 * @param cfg         Configuration that determines which functions are
 *                    supported.
 * @param spec        RealSpec for input/output floating point numbers.
 * @param stage       pipeline registers of this module.
 * @param preMulImpl  structure of the product of the pre-process multiplier.
 * @param postMulImpl structure of the product of the post-process multiplier.
 *                    The two multipliers are chosen separately because they
 *                    are not the same shape: the pre-process one is a wide
 *                    operand times a mantissa, the post-process one two
 *                    mantissas followed by the rounding. The polynomial keeps
 *                    a plain product; its operands fit one hard multiplier.
 *
 * @see [[rial.arith.MulImpl]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class MathFunctions(
  val cfg:         MathFuncConfig,
  val spec:        RealSpec,
  val stage:       PipelineStageConfig,
  val preMulImpl:  MulImpl = MulImpl.plain,
  val postMulImpl: MulImpl = MulImpl.plain
) extends Module {

  import FuncKind._

  assert(!spec.disableSign, "assuming the float type allows negative values")

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val bitW   = cfg.bitW
  val coeffs = cfg.coeffs

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = MathFunctions.stagePoints(cfg, preMulImpl, postMulImpl)

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
    stages.latency("preMul")
  } else {
    0
  }

  /** width of the left hand operand of the pre-process multiplier. */
  val preMulLhsW = PreProcMultiplier.lhsW(cfg.funcs, spec)

  val preMul: Option[Multiplier] = if (needPreMul) { Some(
    Module(new Multiplier(preMulLhsW, 1 + manW, preMulImpl, stages.sub("preMul")))
  )} else { None }

  /** the operand set one function offers to the pre-process multiplier. */
  private def preMulOperands = new MathFunctions.MulOperands(preMulLhsW, 1 + manW)

  // Pre-mult operands, one set per function: each function that wants the
  // multiplier drives its own set in its block below and offers 0 otherwise.
  // Only the selected function's set is ever asked for, so the multiplier
  // input is one flat selection on `selIn`.
  val preMulArgs: Option[Map[FuncKind, Valid[MathFunctions.MulOperands]]] =
    if (needPreMul) {
      Some(cfg.funcs.map(fn => {
        fn -> WireDefault(0.U.asTypeOf(Valid(preMulOperands)))
      }).toMap)
    } else {
      None
    }

  if (preMul.isDefined) {
    val args = MuxLookup(selIn, 0.U.asTypeOf(Valid(preMulOperands)))(
      cfg.funcs.map(fn => cfg.signal(fn) -> preMulArgs.get(fn)))
    preMul.get.io.in.valid  := args.valid
    preMul.get.io.in.bits.x := args.bits.lhs
    preMul.get.io.in.bits.y := args.bits.rhs
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
    stages.latency("postMul")
  } else {
    0
  }

  val nMulArgs = if (needPostMul) {
    stages.at("mulArgs")
  } else {
    0
  }

  /** latency from the polynomial output to the post-process input. */
  val nMulToPost = nMulArgs + nPostMul

  val postMul: Option[PostProcMultiplier] = if (needPostMul) { Some(
    Module(new PostProcMultiplier(
      fracW     = bitW.outBitW,
      manW      = manW,
      roundSpec = RoundSpec.roundToEven,
      stage     = stages.sub("postMul"),
      mulImpl   = postMulImpl,
    ))
  )} else { None }

  /** the operand set one function offers to the post-process multiplier. */
  private def postMulOperands = new MathFunctions.MulOperands(
    1 + bitW.outBitW, 1 + manW)

  // Post-mult operands, one set per function, driven the same way as the
  // pre-mult ones; the selection itself sits at the "mulArgs" boundary below,
  // where the function signal has reached the polynomial output.
  val postMulArgs: Option[Map[FuncKind, Valid[MathFunctions.MulOperands]]] =
    if (needPostMul) {
      Some(cfg.funcs.map(fn => {
        fn -> WireDefault(0.U.asTypeOf(Valid(postMulOperands)))
      }).toMap)
    } else {
      None
    }

  // ==========================================================================
  // stage "pre": the argument reduction of every function
  // ==========================================================================

  val nPre = stages.latency("pre")

  private val cutPre = new SharedPipe("pre", cutPreMul.valid, nPre)

  val selTab = cutPre(selPre)
  val xTab   = cutValid(cutPre, xPre)

  // ==========================================================================
  // stage "table": the coefficient table lookup of every function
  // ==========================================================================

  val nTable = stages.at("table")

  private val cutTable = new SharedPipe("table", cutPre.valid, nTable)

  val selCalc = cutTable(selTab)

  val nCoef = stages.at("coef")

  // ==========================================================================
  // Polynomial Evaluator (shared), stages "poly.*"
  // ==========================================================================

  val polynomialEval = Module(new PolynomialEval(bitW, stages.sub("poly")))

  val nPoly = stages.latency("poly")

  /** latency from the pre-process output to the polynomial output. */
  val nTabToPoly = nTable + nCoef + nPoly

  /** latency from the pre-process output to the post-process input. */
  val nTabToPost = nTabToPoly + nMulToPost

  val polynomialCoefW = bitW.coefBits.reduce(_+_)
  val polynomialCoefs = cfg.funcs.map(fn => {
    fn -> WireDefault(0.U.asTypeOf(Valid(UInt(polynomialCoefW.W))))
  }).toMap

  val polynomialDxs = if (bitW.nOrder == 0) { None } else { Some(
    cfg.funcs.map(fn => {
      fn -> WireDefault(0.U.asTypeOf(Valid(UInt(bitW.inBitW.W))))
    }).toMap
  )}

  // one flat selection on the function signal, at the two boundaries the
  // coefficients and the interval offset reach it
  val polynomialCoef =
    MuxLookup(selCalc, 0.U.asTypeOf(Valid(UInt(polynomialCoefW.W))))(
      cfg.funcs.map(fn => cfg.signal(fn) -> polynomialCoefs(fn)))

  val polynomialDx = if (bitW.nOrder == 0) { None } else { Some(
    Pipe(MuxLookup(selTab, 0.U.asTypeOf(Valid(UInt(bitW.inBitW.W))))(
      cfg.funcs.map(fn => cfg.signal(fn) -> polynomialDxs.get(fn))
    ), nTable)
  )}

  // ==========================================================================
  // stage "coef": the selected coefficients, on their way to the polynomial
  // ==========================================================================

  private val cutCoef = new SharedPipe("coef", cutTable.valid, nCoef)

  val selPoly = cutCoef(selCalc)

  private val polyInValid = if (bitW.nOrder == 0) {
    polynomialCoef.valid
  } else {
    polynomialCoef.valid && polynomialDx.get.valid
  }

  private val cutCoefIn = new SharedPipe("coefIn", polyInValid, nCoef)

  polynomialEval.io.in.valid       := cutCoefIn.valid
  polynomialEval.io.in.bits.coeffs := cutCoefIn(polynomialCoef.bits).asTypeOf(
    MixedVec(bitW.coefBits.map(w => SInt(w.W)))
  )
  if(bitW.nOrder != 0) {
    polynomialEval.io.in.bits.x.get := cutCoefIn(polynomialDx.get.bits).asSInt
  }

  // ==========================================================================
  // stage "poly": everything that bypasses the polynomial
  // ==========================================================================

  private val cutPoly = new SharedPipe("poly", cutCoef.valid, nPoly)

  val selPost = cutPoly(selPoly)

  // ==========================================================================
  // stages "mulArgs" and "postMul": everything that bypasses the post-process
  // multiplier
  // ==========================================================================

  // the operand set of the function `selPost` names, on its way to the
  // post-process multiplier
  private val postMulIn = postMulArgs.map(args => {
    MuxLookup(selPost, 0.U.asTypeOf(Valid(postMulOperands)))(
      cfg.funcs.map(fn => cfg.signal(fn) -> args(fn)))
  })

  // stage "mulArgs": the selected operands
  private val cutMulArgsIn = postMulIn.map(in => {
    new SharedPipe("mulArgsIn", in.valid, nMulArgs)
  })

  if (postMul.isDefined) {
    postMul.get.io.in.valid    := cutMulArgsIn.get.valid
    postMul.get.io.in.bits.lhs := cutMulArgsIn.get(postMulIn.get.bits.lhs)
    postMul.get.io.in.bits.rhs := cutMulArgsIn.get(postMulIn.get.bits.rhs)
  }

  private val cutMulArgs = new SharedPipe("mulArgs", cutPoly.valid, nMulArgs)

  val selMul = cutMulArgs(selPost)

  private val cutPostMul = new SharedPipe("postMul", cutMulArgs.valid, nPostMul)

  val selPostMul = cutPostMul(selMul)

  // the polynomial result as the functions that do not use the post-process
  // multiplier see it, so that every post-process sits at the same boundary
  private val cutPostMulRes = new SharedPipe("postMulRes",
    polynomialEval.io.out.valid, nMulToPost)

  val zresPostMul = cutValid(cutPostMulRes, polynomialEval.io.out)

  // ==========================================================================
  // stage "post": the result assembly of every function
  // ==========================================================================

  val nPost = stages.latency("post")

  private val cutPost = new SharedPipe("post", cutPostMul.valid, nPost)

  val selOut = cutPost(selPostMul)

  // ==========================================================================
  // Output selection, stage "out"
  // ==========================================================================

  val zs = cfg.funcs.map(fn => {
    fn -> WireDefault(0.U.asTypeOf(Valid(UInt(spec.W.W))))
  }).toMap

  // `io.out.valid` is `cutPost.valid`, which is low in every cycle in which
  // the selected function's `zs` is, so the result is selected on `selOut`
  // alone.
  val zMux = MuxLookup(selOut, 0.U(spec.W.W))(
    cfg.funcs.map(fn => cfg.signal(fn) -> zs(fn).bits))

  io.out := Pipe(cutPost.valid, zMux, stages.at("out"))

  // ==========================================================================
  // Per-function pre/tab/nontab/post modules
  //
  // `.collect` so that some funcs does not have the corresponding modules
  // ==========================================================================

  /** the boundaries of the pre-process of `fn`, as the module sees them. */
  private def preStage(fn: FuncKind): PipelineStageConfig = {
    stages.sub("pre." + FuncKind.pointName(fn))
  }

  /** the boundaries of the post-process of `fn`, as the module sees them. */
  private def postStage(fn: FuncKind): PipelineStageConfig = {
    stages.sub("post." + FuncKind.pointName(fn))
  }

  val pres: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtPreProcess       (spec, bitW, preStage(Sqrt       )))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtPreProcess    (spec, bitW, preStage(InvSqrt    )))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalPreProcess (spec, bitW, preStage(Reciprocal )))
    case Sin         => Sin         -> Module(new SinCosPreProcess     (spec, bitW, preStage(Sin        )))
    case Exp         => Exp         -> Module(new ExpPreProcess        (spec, bitW, preStage(Exp        )))
    case Log         => Log         -> Module(new LogPreProcess        (spec, bitW, preStage(Log        )))
    case ACosPhase1  => ACosPhase1  -> Module(new ACosPhase1PreProcess (spec, bitW, preStage(ACosPhase1 )))
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2PreProcess (spec, bitW, preStage(ACosPhase2 )))
    case ATan2Phase1 => ATan2Phase1 -> Module(new ATan2Phase1PreProcess(spec, bitW, preStage(ATan2Phase1)))
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2PreProcess(spec, bitW, preStage(ATan2Phase2)))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidPreProcess    (spec, bitW, preStage(Sigmoid    )))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusPreProcess   (spec, bitW, preStage(SoftPlus   )))
    case ScaleMixtureGaussian => {
      val (sA, sB) = cfg.scaleMixtureGaussianSigma.get
      ScaleMixtureGaussian -> Module(new ScaleMixtureGaussianPreProcess(
        sA, sB, spec, bitW, preStage(ScaleMixtureGaussian)))
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

  // the non-table modules run beside the table, the coefficient selection and
  // the polynomial, and place their registers themselves
  private val nonTableStage = PipelineStageConfig.default(nTabToPoly)

  val nontabs: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt       => Sqrt       -> Module(new SqrtNonTable      (spec, bitW, nonTableStage))
    case InvSqrt    => InvSqrt    -> Module(new InvSqrtNonTable   (spec, bitW, nonTableStage))
    case Reciprocal => Reciprocal -> Module(new ReciprocalNonTable(spec, bitW, nonTableStage))
    // sincos / acos do not have a dedicated non-table path
    case Exp        => Exp        -> Module(new ExpNonTable       (spec, bitW, nonTableStage))
    case Log        => Log        -> Module(new LogNonTable       (spec, bitW, nonTableStage))
    // ACosPhase1 reuses SqrtNonTable
    // ACosPhase2 doesn't need one.
  }.toMap

  val posts: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtPostProcess       (spec, bitW, postStage(Sqrt       )))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtPostProcess    (spec, bitW, postStage(InvSqrt    )))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalPostProcess (spec, bitW, postStage(Reciprocal )))
    case Sin         => Sin         -> Module(new SinCosPostProcess     (spec, bitW, postStage(Sin        )))
    case Exp         => Exp         -> Module(new ExpPostProcess        (spec, bitW, postStage(Exp        )))
    case Log         => Log         -> Module(new LogPostProcess        (spec, bitW, postStage(Log        )))
    case ACosPhase1  => ACosPhase1  -> Module(new ACosPhase1PostProcess (spec, bitW, postStage(ACosPhase1 )))
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2PostProcess (spec, bitW, postStage(ACosPhase2 )))
    case ATan2Phase1 => ATan2Phase1 -> Module(new ATan2Phase1PostProcess(spec, bitW, postStage(ATan2Phase1)))
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2PostProcess(spec, bitW, postStage(ATan2Phase2)))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidPostProcess    (spec, bitW, postStage(Sigmoid    )))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusPostProcess   (spec, bitW, postStage(SoftPlus   )))
    case ScaleMixtureGaussian => ScaleMixtureGaussian ->
      Module(new ScaleMixtureGaussianPostProcess(spec, bitW, postStage(ScaleMixtureGaussian)))
  }.toMap

  private def toModule[T <: Module](x: Module)(implicit ct: ClassTag[T]): T = {
    require(ct.runtimeClass.isInstance(x),
      s"expected ${ct.runtimeClass.getSimpleName} but got ${x.getClass.getSimpleName}")
    x.asInstanceOf[T]
  }

  // --------------------------------------------------------------------------

  /** [[SqrtNonTable]]'s output at the post-process input.  Sqrt and ACosPhase1
   *  read the same net, so it is delayed once. */
  private val sqrtNtabAtPost = if (cfg.has(Sqrt)) {
    Some(Pipe(toModule[SqrtNonTable](nontabs(Sqrt)).io.z, nMulToPost))
  } else {
    None
  }

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
      polynomialDxs.get(Sqrt) := pre.io.dx.get
    }

    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(Sqrt)
    ntab.io.x.bits  := xTab.bits

    val ntabAtPost = sqrtNtabAtPost.get

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
      polynomialDxs.get(InvSqrt) := pre.io.dx.get
    }

    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(InvSqrt)
    ntab.io.x.bits  := xTab.bits

    val ntabAtPost = Pipe(ntab.io.z, nMulToPost)

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
      polynomialDxs.get(Reciprocal) := pre.io.dx.get
    }

    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(Reciprocal)
    ntab.io.x.bits  := xTab.bits

    val ntabAtPost = Pipe(ntab.io.z, nMulToPost)

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

    // ----- pre-mult operands: 1/π * (1 + x.man), the same set for both -----
    Seq(Sin, Cos).foreach(fn => {
      preMulArgs.get(fn).valid    := xIn.valid
      preMulArgs.get(fn).bits.lhs := SinCosPreMulArgs.lhs(spec, preMulLhsW)
      preMulArgs.get(fn).bits.rhs := Cat(1.U(1.W), xIn.bits.man)
    })

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
      polynomialDxs.get(Sin) := sincosPre.io.dx.get
      polynomialDxs.get(Cos) := sincosPre.io.dx.get
    }

    // ----- post-mult operands: polynomial * (1 + yman) -----
    val sincosInfoAtPost = Pipe(sincosPre.io.info, nTabToPoly)

    Seq(Sin, Cos).foreach(fn => {
      postMulArgs.get(fn).valid    :=
        polynomialEval.io.out.valid && sincosInfoAtPost.valid
      postMulArgs.get(fn).bits.lhs := Cat(1.U(1.W), polynomialEval.io.out.bits)
      postMulArgs.get(fn).bits.rhs := Cat(1.U(1.W), sincosInfoAtPost.bits.yman)
    })

    // ----- post -----
    sincosPost.io.pre := Pipe(sincosInfoAtPost, nMulToPost)
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

    // ----- pre-mult operands: log2(e) * (1 + x.man) -----
    preMulArgs.get(Exp).valid    := xIn.valid
    preMulArgs.get(Exp).bits.lhs := ExpPreMulArgs.lhs(spec, preMulLhsW)
    preMulArgs.get(Exp).bits.rhs := Cat(1.U(1.W), xIn.bits.man)

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
      polynomialDxs.get(Exp) := pre.io.dx.get
    }

    // ----- non-table ------
    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(Exp)
    ntab.io.x.bits  := xTab.bits
    ntab.io.info := pre.io.info

    // ----- post -----
    val ntabAtPost = Pipe(ntab.io.z, nMulToPost)

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
      polynomialDxs.get(Log) := pre.io.dx.get
    }

    ntab.io.x.valid := xTab.valid && selTab === cfg.signal(Log)
    ntab.io.x.bits  := xTab.bits

    // ----- post-mult arg builder (combinational, at the polynomial output) -----
    val mulArgs = Module(new LogPostMultArgs(spec, bitW))
    mulArgs.io.info := ntab.io.z
    mulArgs.io.zres := polynomialEval.io.out.bits

    postMulArgs.get(Log).valid    :=
      mulArgs.io.mul.valid && polynomialEval.io.out.valid
    postMulArgs.get(Log).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(Log).bits.rhs := mulArgs.io.mul.bits.rhs

    // ----- post -----
    post.io.info := Pipe(ntab.io.z, nMulToPost)
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
      polynomialDxs.get(ACosPhase1) := pre.io.dx.get
    }

    // pre.io.flags leaves the pre-process at the table input; pipeline it to
    // the input of the post-process
    val flagsAtPost = Pipe(pre.io.flags, nTabToPost)
    val ntabAtPost  = sqrtNtabAtPost.get

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
      polynomialDxs.get(ACosPhase2) := pre.io.dx.get
    }

    // post-mult args at the polynomial output: lhs = polynomial w/ hidden 1,
    // rhs = x.manW1.
    val xManW1AtPost = Pipe(pre.io.xManW1, nTabToPoly)
    val mulArgs = Module(new ACosPhase2PostMultArgs(spec, bitW))
    mulArgs.io.zres   := polynomialEval.io.out.bits
    mulArgs.io.xManW1 := xManW1AtPost

    postMulArgs.get(ACosPhase2).valid    :=
      mulArgs.io.mul.valid && polynomialEval.io.out.valid
    postMulArgs.get(ACosPhase2).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(ACosPhase2).bits.rhs := mulArgs.io.mul.bits.rhs

    // flags & x.ex pipelined to the input of the post-process
    val flagsAtAfterMul = Pipe(pre.io.flags, nTabToPost)
    val xexAtAfterMul   = Pipe(pre.io.xex,   nTabToPost)

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
      polynomialDxs.get(ATan2Phase1) := pre.io.dx.get
    }

    // ---- post-mult args at the polynomial output: lhs (poly+hidden 1), rhs (1+minMan) ----
    val infoAtPost = Pipe(pre.io.info, nTabToPoly)
    val mulArgs = Module(new ATan2Phase1PostMultArgs(spec, bitW))
    mulArgs.io.zres      := polynomialEval.io.out.bits
    mulArgs.io.maxXYMan0 := infoAtPost.bits.maxXYMan0
    mulArgs.io.minMan    := infoAtPost.bits.minMan
    mulArgs.io.valid     := infoAtPost.valid && polynomialEval.io.out.valid

    postMulArgs.get(ATan2Phase1).valid    := mulArgs.io.mul.valid
    postMulArgs.get(ATan2Phase1).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(ATan2Phase1).bits.rhs := mulArgs.io.mul.bits.rhs

    // info & flags pipelined to the input of the post-process; the
    // post-process reads only zex and xySameMan
    val postInfo = Wire(Valid(new ATan2Phase1PostInfo(spec)))
    postInfo.valid          := infoAtPost.valid
    postInfo.bits.zex       := infoAtPost.bits.zex
    postInfo.bits.xySameMan := infoAtPost.bits.xySameMan

    val infoAtAfterMul  = Pipe(postInfo,     nMulToPost)
    val flagsAtAfterMul = Pipe(pre.io.flags, nTabToPost)

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
      polynomialDxs.get(ATan2Phase2) := pre.io.dx.get
    }

    // ---- post-mult args at the polynomial output ----
    val xManW1AtPost   = Pipe(pre.io.xManW1,   nTabToPoly)
    val isLinearAtPost = Pipe(pre.io.isLinear, nTabToPoly)
    val mulArgs = Module(new ATan2Phase2PostMultArgs(spec, bitW))
    mulArgs.io.zres     := polynomialEval.io.out.bits
    mulArgs.io.isLinear := isLinearAtPost.bits
    mulArgs.io.xManW1   := xManW1AtPost

    postMulArgs.get(ATan2Phase2).valid    :=
      mulArgs.io.mul.valid && polynomialEval.io.out.valid
    postMulArgs.get(ATan2Phase2).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(ATan2Phase2).bits.rhs := mulArgs.io.mul.bits.rhs

    // flags, xex, isLinear pipelined to the input of the post-process
    val flagsAtAfterMul    = Pipe(pre.io.flags,    nTabToPost)
    val xexAtAfterMul      = Pipe(pre.io.xex,      nTabToPost)
    val isLinearAtAfterMul = Pipe(isLinearAtPost,  nMulToPost)

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
      polynomialDxs.get(Sigmoid) := pre.io.dx.get
    }

    // info pipelined from the table input to the input of the post-process
    val infoAtPost = Pipe(pre.io.info, nTabToPost)
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
      polynomialDxs.get(SoftPlus) := pre.io.dx.get
    }

    val infoAtPost = Pipe(pre.io.info, nTabToPost)
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

    // ----- pre-mult operands: (1/sgmA^2) * (1 + x.man) -----
    preMulArgs.get(ScaleMixtureGaussian).valid    := xIn.valid
    preMulArgs.get(ScaleMixtureGaussian).bits.lhs :=
      ScaleMixtureGaussianPreMulArgs.lhs(sgmA, spec, preMulLhsW)
    preMulArgs.get(ScaleMixtureGaussian).bits.rhs := Cat(1.U(1.W), xIn.bits.man)

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
      polynomialDxs.get(ScaleMixtureGaussian) := pre.io.dx.get
    }

    // ----- post-mult args at the polynomial output -----
    val infoAtPost = Pipe(pre.io.info, nTabToPoly)
    val mulArgs    = Module(new ScaleMixtureGaussianPostMulArgs(sgmA, sgmB, spec, bitW))
    mulArgs.io.useTable  := infoAtPost.bits.useTable
    mulArgs.io.xsgmA2Man := infoAtPost.bits.xsgmA2Man
    mulArgs.io.zres      := polynomialEval.io.out.bits
    mulArgs.io.valid     := infoAtPost.valid && polynomialEval.io.out.valid

    postMulArgs.get(ScaleMixtureGaussian).valid    := mulArgs.io.mul.valid
    postMulArgs.get(ScaleMixtureGaussian).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(ScaleMixtureGaussian).bits.rhs := mulArgs.io.mul.bits.rhs

    // ----- post stage: pipe info & z0Ex to the input of the post-process -----
    val infoAtAfterMul = Pipe(infoAtPost, nMulToPost)
    val z0exAtAfterMul = Pipe(mulArgs.io.mul.valid, mulArgs.io.mul.bits.z0ex,
                              nMulToPost)

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

  /** The operands one function offers to a shared multiplier.
   *
   * Every configured function drives its own set and the set of the selected
   * function reaches the multiplier, so the operands are picked by one flat
   * mux on the function signal.
   *
   * @param lhsW width of the left operand
   * @param rhsW width of the right operand
   */
  class MulOperands(val lhsW: Int, val rhsW: Int) extends Bundle {
    val lhs = UInt(lhsW.W)
    val rhs = UInt(rhsW.W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   * [[rial.math.MathFunctions]]. The argument reduction and the result
   * assembly of the configured functions appear as the parallels `"pre"` and
   * `"post"`, the polynomial and the multipliers as the groups `"poly"`,
   * `"preMul"` and `"postMul"`.
   *
   * `"preMul"`, `"mulArgs"` and `"postMul"` are declared only when one of the
   * functions of `cfg` needs the corresponding multiplier.
   *
   * @param cfg         configuration of the module these boundaries are for
   * @param preMulImpl  structure of the pre-process multiplier's product
   * @param postMulImpl structure of the post-process multiplier's product
   */
  def stagePoints(cfg: MathFuncConfig,
    preMulImpl:  MulImpl = MulImpl.plain,
    postMulImpl: MulImpl = MulImpl.plain): Seq[StageElem] = {
    val preMul: Seq[StageElem] = if (PreProcMultiplier.isNeeded(cfg.funcs)) {
      Seq(StageGroup("preMul", Multiplier.stagePoints(preMulImpl)))
    } else {
      Seq.empty[StageElem]
    }
    val postMul: Seq[StageElem] = if (PostProcMultiplier.isNeeded(cfg.funcs)) {
      Seq(StagePoint("mulArgs", 1.0),
          StageGroup("postMul", PostProcMultiplier.stagePoints(postMulImpl)))
    } else {
      Seq.empty[StageElem]
    }
    preMul ++
    Seq(StageParallel("pre", prePaths(cfg)),
        StagePoint("table", 2.0),
        StagePoint("coef",  1.0),
        StageGroup("poly",  PolynomialEval.stagePoints(cfg.bitW.nOrder))) ++
    postMul ++
    Seq(StageParallel("post", postPaths(cfg)),
        StagePoint("out",  1.0))
  }

  /** The argument reduction of each configured function, in `cfg.funcs` order.
   *
   * `Sin` and `Cos` share one pre-process, so they contribute one path.
   *
   * @param cfg configuration of the module these paths are for
   */
  def prePaths(cfg: MathFuncConfig): Seq[StageGroup] = {
    import FuncKind._
    cfg.funcs.collect {
      case Sqrt        => path(Sqrt,        SqrtPreProcess               .stagePoints)
      case InvSqrt     => path(InvSqrt,     InvSqrtPreProcess            .stagePoints)
      case Reciprocal  => path(Reciprocal,  ReciprocalPreProcess         .stagePoints)
      case Sin         => path(Sin,         SinCosPreProcess             .stagePoints)
      case ACosPhase1  => path(ACosPhase1,  ACosPhase1PreProcess         .stagePoints)
      case ACosPhase2  => path(ACosPhase2,  ACosPhase2PreProcess         .stagePoints)
      case ATan2Phase1 => path(ATan2Phase1, ATan2Phase1PreProcess        .stagePoints)
      case ATan2Phase2 => path(ATan2Phase2, ATan2Phase2PreProcess        .stagePoints)
      case Exp         => path(Exp,         ExpPreProcess                .stagePoints)
      case Log         => path(Log,         LogPreProcess                .stagePoints)
      case Sigmoid     => path(Sigmoid,     SigmoidPreProcess            .stagePoints)
      case SoftPlus    => path(SoftPlus,    SoftPlusPreProcess           .stagePoints)
      case ScaleMixtureGaussian =>
                          path(ScaleMixtureGaussian,
                               ScaleMixtureGaussianPreProcess.stagePoints)
    }
  }

  /** The result assembly of each configured function, in `cfg.funcs` order.
   *
   * `Sin` and `Cos` share one post-process, so they contribute one path.
   *
   * @param cfg configuration of the module these paths are for
   */
  def postPaths(cfg: MathFuncConfig): Seq[StageGroup] = {
    import FuncKind._
    cfg.funcs.collect {
      case Sqrt        => path(Sqrt,        SqrtPostProcess               .stagePoints)
      case InvSqrt     => path(InvSqrt,     InvSqrtPostProcess            .stagePoints)
      case Reciprocal  => path(Reciprocal,  ReciprocalPostProcess         .stagePoints)
      case Sin         => path(Sin,         SinCosPostProcess             .stagePoints)
      case ACosPhase1  => path(ACosPhase1,  ACosPhase1PostProcess         .stagePoints)
      case ACosPhase2  => path(ACosPhase2,  ACosPhase2PostProcess         .stagePoints)
      case ATan2Phase1 => path(ATan2Phase1, ATan2Phase1PostProcess        .stagePoints)
      case ATan2Phase2 => path(ATan2Phase2, ATan2Phase2PostProcess        .stagePoints)
      case Exp         => path(Exp,         ExpPostProcess                .stagePoints)
      case Log         => path(Log,         LogPostProcess                .stagePoints)
      case Sigmoid     => path(Sigmoid,     SigmoidPostProcess            .stagePoints)
      case SoftPlus    => path(SoftPlus,    SoftPlusPostProcess           .stagePoints)
      case ScaleMixtureGaussian =>
                          path(ScaleMixtureGaussian,
                               ScaleMixtureGaussianPostProcess.stagePoints)
    }
  }

  /** One path of the `"pre"` or `"post"` parallel, named after its function. */
  private def path(fn: FuncKind.FuncKind, points: Seq[StageElem]): StageGroup = {
    StageGroup(FuncKind.pointName(fn), points)
  }
}
