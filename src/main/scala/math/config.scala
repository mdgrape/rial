package rial.math

import chisel3._
import chisel3.util._

import rial.arith._
import rial.util.PipelineStageConfig
import rial.util.PipelineStageConfig._

import scala.math.exp
import scala.collection.mutable.{ArrayBuffer, HashMap}

/** Enumerator to represent which function is used in [[rial.math.MathFunctions]]. */
object FuncKind extends Enumeration {
  type FuncKind = Value

  /** sqrt(x) */
  val Sqrt        = Value
  /** 1 / sqrt(x) */
  val InvSqrt     = Value
  /** 1 / x */
  val Reciprocal  = Value
  /** sin(x) */
  val Sin         = Value
  /** cos(x) */
  val Cos         = Value
  /** acos(x) is calculated in 2 phases. it represents the first phase of acos */
  val ACosPhase1  = Value
  /** acos(x) is calculated in 2 phases. it represents the second phase of acos */
  val ACosPhase2  = Value
  /** atan2(y, x) is calculated in 2 phases. it represents the second phase of atan2 */
  val ATan2Phase1 = Value
  /** atan2(y, x) is calculated in 2 phases. it represents the second phase of atan2 */
  val ATan2Phase2 = Value
  /** exp(x) */
  val Exp         = Value
  /** log(x). base e. */
  val Log         = Value
  /** 1 / (1+exp(-x)) */
  val Sigmoid     = Value
  /** log(1 + exp(x)) */
  val SoftPlus    = Value
  val ScaleMixtureGaussian = Value

  /** Returns enum name (e.g. FuncKind.Sqrt -> "Sqrt").
   */
  def getString(fn: FuncKind): String = {
    fn match {
      case Sqrt => "Sqrt"
      case InvSqrt => "InvSqrt"
      case Reciprocal => "Reciprocal"
      case Sin => "Sin"
      case Cos => "Cos"
      case ACosPhase1  => "ACosPhase1"
      case ACosPhase2  => "ACosPhase2"
      case ATan2Phase1 => "ATan2Phase1"
      case ATan2Phase2 => "ATan2Phase2"
      case Exp => "Exp"
      case Log => "Log"
      case Sigmoid => "Sigmoid"
      case SoftPlus => "SoftPlus"
      case ScaleMixtureGaussian => "ScaleMixtureGaussian"
    }
  }

  /** (Internal use) Returns true if `fn` requires preprocess-multiplier.
   */
  def needPreMult(fn: FuncKind): Boolean = {
    fn == Exp || fn == Sin || fn == Cos || fn == ScaleMixtureGaussian
  }

  /** (Internal use) Returns true if `fn` requires postprocess-multiplier.
   */
  def needPostMult(fn: FuncKind): Boolean = {
    fn == ACosPhase2 || fn == ATan2Phase1 || fn == ATan2Phase2 ||
    fn == Log || fn == Sin || fn == Cos || fn == ScaleMixtureGaussian
  }

  /** Returns list of funcs required to support that function.
   *
   *  Some functions perform variable transformations to calculate more accurate
   *  approximations, utilizing the approximation tables of other functions.
   */
  def requiredFuncs(fn: FuncKind): Seq[FuncKind] = {
    fn match {
      case Sqrt                 => Seq(Sqrt)
      case InvSqrt              => Seq(InvSqrt)
      case Reciprocal           => Seq(Reciprocal)
      case Sin                  => Seq(Sin, Cos)
      case Cos                  => Seq(Sin, Cos)
      case ACosPhase1           => Seq(ACosPhase1,  ACosPhase2, Sqrt)
      case ACosPhase2           => Seq(ACosPhase1,  ACosPhase2, Sqrt)
      case ATan2Phase1          => Seq(ATan2Phase1, ATan2Phase2, Reciprocal)
      case ATan2Phase2          => Seq(ATan2Phase1, ATan2Phase2, Reciprocal)
      case Exp                  => Seq(Exp)
      case Log                  => Seq(Log)
      case Sigmoid              => Seq(Sigmoid)
      case SoftPlus             => Seq(SoftPlus)
      case ScaleMixtureGaussian => Seq(ScaleMixtureGaussian)
    }
  }

  /** Returns list of funcs that are required by specified functions.
   */
  def normalize(funcs: Seq[FuncKind]): Seq[FuncKind] = {
    funcs.map(f => requiredFuncs(f)).reduce(_++_).distinct.sorted
  }
}

/** A Config class for [[rial.math.MathFunctions]] Module.

 *  @constructor create a new MathFuncConfig.
 *  @param funcs the list of functions that should be supported.
 *  @param ScaleMixtureGaussianSigma sigmaA and sigmaB for ScaleMixtureGaussian.
 */
class MathFuncConfig(
  val funcs: Seq[FuncKind.FuncKind],
  val scaleMixtureGaussianSigma: Option[(Double, Double)] = None,
) {

  assert(funcs.length > 0, "At least one function should be supported")
  assert(funcs == FuncKind.normalize(funcs),
    """| some function depends on each other. like, ATan2 reuqires Reciprocal
       | function. do `normalize` before passing it.
       |""".stripMargin)

  import FuncKind._

  /** Checks if a function is supported.
   *
   *  @param fn An enumerator of the function.
   *  @return true if the function is supported. false if not.
   */
  def has(fn: FuncKind): Boolean = {
    funcs.exists(_==fn)
  }

  if(has(ScaleMixtureGaussian)) {
    assert(scaleMixtureGaussianSigma.isDefined)
  }

  // ---------------------------------------------------------------------------

  /** (Internal use) Determines bit widths of the polynomial coefficients of
   *  the function.
   *
   *  @return bit width of polynomial coefficients of the function per order.
   */
  def getCBits(fn: FuncKind, spec: RealSpec, polySpec: PolynomialSpec): Seq[Int] = {
    fn match {
      case Sqrt        => SqrtTableCoeff       .getCBits(spec, polySpec)
      case InvSqrt     => InvSqrtTableCoeff    .getCBits(spec, polySpec)
      case Reciprocal  => ReciprocalTableCoeff .getCBits(spec, polySpec)
      case Sin         => SinCosTableCoeff     .getCBits(spec, polySpec)
      case Cos         => SinCosTableCoeff     .getCBits(spec, polySpec)
      case ACosPhase1  => SqrtTableCoeff       .getCBits(spec, polySpec)
      case ACosPhase2  => ACosTableCoeff       .getCBits(spec, polySpec)
      case ATan2Phase1 => ReciprocalTableCoeff .getCBits(spec, polySpec)
      case ATan2Phase2 => ATan2Phase2TableCoeff.getCBits(spec, polySpec)
      case Exp         => ExpTableCoeff        .getCBits(spec, polySpec)
      case Log         => LogTableCoeff        .getCBits(spec, polySpec)
      case Sigmoid     => SigmoidTableCoeff    .getCBits(spec, polySpec)
      case SoftPlus    => SoftPlusTableCoeff   .getCBits(spec, polySpec)
      case ScaleMixtureGaussian => {
        val (sA, sB) = scaleMixtureGaussianSigma.get
        ScaleMixtureGaussianTableCoeff.getCBits(sA, sB, spec, polySpec)
      }
    }
  }

  /** (Internal use) Determines the bit width of the temporary used to calculate
   *  polynomial.
   *
   *  @return bit width of temporary used in the polynomial of the function.
   */
  def getCalcW(fn: FuncKind, spec: RealSpec, polySpec: PolynomialSpec): Seq[Int] = {
    fn match {
      case Sqrt        => SqrtTableCoeff       .getCalcW(spec, polySpec)
      case InvSqrt     => InvSqrtTableCoeff    .getCalcW(spec, polySpec)
      case Reciprocal  => ReciprocalTableCoeff .getCalcW(spec, polySpec)
      case Sin         => SinCosTableCoeff     .getCalcW(spec, polySpec)
      case Cos         => SinCosTableCoeff     .getCalcW(spec, polySpec)
      case ACosPhase1  => SqrtTableCoeff       .getCalcW(spec, polySpec)
      case ACosPhase2  => ACosTableCoeff       .getCalcW(spec, polySpec)
      case ATan2Phase1 => ReciprocalTableCoeff .getCalcW(spec, polySpec)
      case ATan2Phase2 => ATan2Phase2TableCoeff.getCalcW(spec, polySpec)
      case Exp         => ExpTableCoeff        .getCalcW(spec, polySpec)
      case Log         => LogTableCoeff        .getCalcW(spec, polySpec)
      case Sigmoid     => SigmoidTableCoeff    .getCalcW(spec, polySpec)
      case SoftPlus    => SoftPlusTableCoeff   .getCalcW(spec, polySpec)
      case ScaleMixtureGaussian => {
        val (sA, sB) = scaleMixtureGaussianSigma.get
        ScaleMixtureGaussianTableCoeff.getCalcW(sA, sB, spec, polySpec)
      }
    }
  }

  // --------------------------------------------------------------------------

  /** The width of UInt to represent function select signal.
   *
   *  signal contains: `None ++ [funcs..] ++ Invalid`.
   *  `io.sel` of each functions corresponds to non-zero index in `funcs`.
   *  If `None` is given, it calculates nothing.
   *  If `Invalid` is given, it calculates something invalid.
   *
   *  `Invalid` is introduced to make implementation easy. There is no sense in
   *  sending `Invalid` as `io.sel` signal.
   */
  val signalW = log2Up(2 + funcs.length)

  /** Returns function select signal to be passed to [[rial.math.MathFunction]].
   *
   * {{{
   * val mathfunc = Module(new MathFunction(...))
   *
   * when(doSqrt) {
   *   mathfunc.io.sel := MathFuncConfig.signal(FuncKind.Sqrt)
   * }.otherwise {
   *   mathfunc.io.sel := MathFuncConfig.signalNone
   * }
   * }}}
   *
   * If the passed function is not supported, it returns invalid value (`max+1`).
   *
   * @param fn An enumerator of the function
   * @return the signal that corresponds to the function
   */
  def signal(fn: FuncKind): UInt = {
    if (has(fn)) {
      (funcs.indexWhere(_==fn) + 1).U(signalW.W)
    } else {
      (funcs.length + 1).U(signalW.W)
    }
  }

  /** Returns function select signal that runs no function.
   *
   *  @return the signal that corresponds to "no function".
   */
  def signalNone(): UInt = {
    0.U(signalW.W)
  }

  def getString: String = {
    funcs.map(fn => FuncKind.getString(fn)).mkString("[", ", ", "]")
  }
}

/** Factory for [[rial.math.MathFuncConfig]].
 */
object MathFuncConfig {
  import FuncKind._

  /** Defines simple functions that does not require multiplier in preprocess / postprocess.
   *
   * includes:
   * <ul>
   *   <li> sqrt       </li>
   *   <li> invsqrt    </li>
   *   <li> reciprocal </li>
   * </ul>
   */
  val simple = new MathFuncConfig(Seq(
    Sqrt, InvSqrt, Reciprocal
  ))

  /** Defines standard math functions that are frequently used.
   *
   * includes:
   * <ul>
   *   <li> sqrt       </li>
   *   <li> invsqrt    </li>
   *   <li> reciprocal </li>
   *   <li> sin        </li>
   *   <li> cos        </li>
   *   <li> acos       </li>
   *   <li> atan2      </li>
   *   <li> exp        </li>
   *   <li> log        </li>
   * </ul>
   */
  val standard = new MathFuncConfig(Seq(
    Sqrt, InvSqrt, Reciprocal, Sin, Cos,
    ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2,
    Exp, Log
  ))

  /** Defines all the supported math functions including too task-specific ones.
   *
   * includes:
   * <ul>
   *   <li> everything. </li>
   * </ul>
   */
  val all = new MathFuncConfig(Seq(
    Sqrt, InvSqrt, Reciprocal, Sin, Cos,
    ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2, Exp, Log,
    Sigmoid, SoftPlus, ScaleMixtureGaussian),
    Some((exp(-1.0), exp(-6.0))
  ))
}

/** Config class to set pipeline stages in a [[rial.math.MathFunctions]] module.
 *
 * Overview:
 *
 * {{{
 * //             .--preStage    .--tableCalcGap .-- calcPostGap
 * //             |      .--preCalcGap   .--calcStage    .--postStage
 * //         ____|____  |       |  _____|_____  |  _____|_____
 * //        '         ' '       ' '           ' ' '           '
 * //       .----.-.----.-.-----.-.-----.-.-----.-.-----.-.-----.
 * //       |    |v|    |v|table|v|Calc1|v|Calc2|v|     |v|     |
 * // in -> |Pre1| |Pre2| :-----'-'-----'-'-----: |Post1| |Post2| -> out
 * //       |    | |    | |      non-table      | |     | |     |
 * //       '----'-'----'-'---------------------'-'-----'-'-----'
 * }}}
 *
 * @constructor create a new MathFuncConfig.
 * @param preStage     pipeline stages of preprocess.
 * @param calcStage    pipeline stages of table/polynomial and non-table path.
 * @param postStage    pipeline stages of postprocess.
 * @param postMulStage pipeline stages of multiplier in postprocess. should be smaller than postStage.
 * @param preCalcGap   if true, add register between preprocess and calculation stage
 * @param tableCalcGap if true, add register between table and calculation stage (+1 to calcStage for OtherPath)
 * @param calcPostGap  if true, add register between calculation and postprocess stage
 *
 */
class MathFuncPipelineConfig(
  val preStage:     PipelineStageConfig,
  val preMulStage:  PipelineStageConfig,
  val calcStage:    PipelineStageConfig,
  val postStage:    PipelineStageConfig,
  val postMulStage: PipelineStageConfig,
  val preCalcGap:   Boolean,
  val tableCalcGap: Boolean,
  val calcPostGap:  Boolean,
  ) {

  assert(preMulStage.total <= preStage.total,
    "preStage includes preMulStage, so should be larger than preMulStage.")
  assert(postMulStage.total <= postStage.total,
    "postStage includes postMulStage, so should be larger than postMulStage.")

  /** Calculates the total latency in clock cycles.
   */
  def total = {
    preStage.total +
    calcStage.total +
    postStage.total +
    (if(preCalcGap)   {1} else {0}) +
    (if(tableCalcGap) {1} else {0}) +
    (if(calcPostGap)  {1} else {0})
  }
  /** Generates a string that represents the current config.
   */
  def getString = {
    f"pre:${preStage.total}(mul=${preMulStage.total}), " +
    f"calc:${calcStage.total}, " +
    f"post:${postStage.total}(mul=${postMulStage.total})" +
    (if(preCalcGap){" pre/c"} else {""}) +
    (if(tableCalcGap){" t/c"} else {""}) +
    (if(calcPostGap){" c/post"} else {""})
  }
}

/** Factory for [[rial.math.MathFuncPipelineConfig]].
 */
object MathFuncPipelineConfig {

  /** constructs a config corresponding to the single cycle mathfunc.
   */
  def none = {
    new MathFuncPipelineConfig(
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      false,
      false,
      false)
  }
}

