package rial.funcs

import chisel3._
import chisel3.util._

import rial.arith.RealSpec
import rial.util.PipelineStageConfig
import rial.funcs._

import scala.collection.mutable.{ArrayBuffer, HashMap}
import scala.math.exp

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
  /** the 1st phase of acos(x).
   *
   * `acos` is calculated in 2 phases.
   * the second phase takes the result of the first phase.
   */
  val ACosPhase1  = Value
  /** the 2nd phase of acos(x).
   *
   * `acos` is calculated in 2 phases.
   * the second phase takes the result of the first phase.
   */
  val ACosPhase2  = Value
  /** the 1st phase of atan2(y,x).
   *
   * `atan2` is calculated in 2 phases.
   * the second phase takes the result of the first phase.
   */
  val ATan2Phase1 = Value
  /** the 2nd phase of atan2(y,x).
   *
   * `atan2` is calculated in 2 phases.
   * the second phase takes the result of the first phase.
   */
  val ATan2Phase2 = Value
  /** exp(x) e^x */
  val Exp         = Value
  /** log(x) base e. */
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
      case Sqrt        => "Sqrt"
      case InvSqrt     => "InvSqrt"
      case Reciprocal  => "Reciprocal"
      case Sin         => "Sin"
      case Cos         => "Cos"
      case ACosPhase1  => "ACosPhase1"
      case ACosPhase2  => "ACosPhase2"
      case ATan2Phase1 => "ATan2Phase1"
      case ATan2Phase2 => "ATan2Phase2"
      case Exp         => "Exp"
      case Log         => "Log"
      case Sigmoid     => "Sigmoid"
      case SoftPlus    => "SoftPlus"
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
   *
   * Some functions requires others. For example, sin and cos shares polynomial
   * table, so Sin requires Cos and Cos requires Sin. Also, the first phase of
   * acos and atan2 requires other math functions, so they requires others.
   * this function lists all the functions required to support the given list
   * of functions.
   */
  def normalize(funcs: Seq[FuncKind]): Seq[FuncKind] = {
    funcs.map(f => requiredFuncs(f)).reduce(_++_).distinct.sorted
  }
}

// ============================================================================

/** A Config class for [[rial.math.MathFunctions]] Module.

 *  @constructor create a new MathFuncConfig.
 *  @param funcs the list of functions that should be supported.
 *  @param ScaleMixtureGaussianSigma sigmaA and sigmaB for ScaleMixtureGaussian.
 */
class MathFuncConfig(
  val funcs:   Seq[FuncKind.FuncKind],
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
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
  val signalW = log2Ceil(2 + funcs.length)

  /** Returns function select signal to be passed to [[rial.math.MathFunctions]].
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

  def isValid(s: UInt): Bool = {
    s =/= 0.U && s < (funcs.length + 1).U
  }

  /** Returns function select signal that runs no function.
   *
   *  @return the signal that corresponds to "no function".
   */
  def signalNone(): UInt = {
    0.U(signalW.W)
  }

  // ---------------------------------------------------------------------------

  val tables: Map[FuncKind, TableBase] = funcs.map {
    case Sqrt        => (Sqrt       -> new SqrtTable       (nOrder, addrW, inBitW, outBitW))
    case InvSqrt     => (InvSqrt    -> new InvSqrtTable    (nOrder, addrW, inBitW, outBitW))
    case Reciprocal  => (Reciprocal -> new ReciprocalTable (nOrder, addrW, inBitW, outBitW))
    case Sin         => (Sin        -> new SinCosTable     (nOrder, addrW, inBitW, outBitW))
    case Cos         => (Cos        -> new SinCosTable     (nOrder, addrW, inBitW, outBitW))
    case Exp         => (Exp        -> new ExpTable        (nOrder, addrW, inBitW, outBitW))
    case Log         => (Log        -> new LogTable        (nOrder, addrW, inBitW, outBitW))
    // Phase 1 reuses Sqrt's polynomial path; we still register a SqrtTable
    // here so the merged bitW accounts for it.  At HW time, the Sqrt block
    // shares its TableCoeff/NonTable with Phase 1 (no extra ROM).
    case ACosPhase1  => (ACosPhase1 -> new SqrtTable       (nOrder, addrW, inBitW, outBitW))
    case ACosPhase2  => (ACosPhase2 -> new ACosPhase2Table (nOrder, addrW, inBitW, outBitW))
    // case ATan2Phase1 => new ReciprocalTable (nOrder, addrW, inBitW, outBitW)
    // case ATan2Phase2 => new ATan2Phase2Table(nOrder, addrW, inBitW, outBitW)
    // case Sigmoid     => new SigmoidTable    (nOrder, addrW, inBitW, outBitW)
    // case SoftPlus    => new SoftPlusTable   (nOrder, addrW, inBitW, outBitW)
    // case ScaleMixtureGaussian => {
    //   val (sA, sB) = scaleMixtureGaussianSigma.get
    //   new ScaleMixtureGaussianTable(sA, sB, nOrder, nIntervals, inBitW, outBitW)
    // }
  }.toMap

  val bitW   = tables.map(_._2.getBitWidth()).reduce((l, r) => { l.merge(r) } )
  val coeffs = tables.map{ case (k,v) => (k -> v.generateTable(bitW)) }.toMap

  // --------------------------------------------------------------------------

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
  val simple = new MathFuncConfig(
    funcs   = Seq(Sqrt, InvSqrt, Reciprocal),
    nOrder  = 2,
    addrW   = 8,
    inBitW  = RealSpec.Float32Spec.manW - 8,
    outBitW = RealSpec.Float32Spec.manW + 3,
  )

  // TODO: uncomment when sin/cos/exp/log are implemented in funcs/
  // /** Defines standard math functions that are frequently used.
  //  *
  //  * includes:
  //  * <ul>
  //  *   <li> sqrt       </li>
  //  *   <li> invsqrt    </li>
  //  *   <li> reciprocal </li>
  //  *   <li> sin        </li>
  //  *   <li> cos        </li>
  //  *   <li> acos       </li>
  //  *   <li> atan2      </li>
  //  *   <li> exp        </li>
  //  *   <li> log        </li>
  //  * </ul>
  //  */
  // val standard = new MathFuncConfig(
  //   funcs = Seq(Sqrt, InvSqrt, Reciprocal, Sin, Cos,
  //     ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2,
  //     Exp, Log),
  //   nOrder = 2,
  //   addrW = 8,
  //   inBitW = RealSpec.Float32Spec.manW - 8,
  //   outBitW = RealSpec.Float32Spec.manW + 3
  // )
  //
  // /** Defines all the supported math functions including too task-specific ones.
  //  *
  //  * includes:
  //  * <ul>
  //  *   <li> everything. </li>
  //  * </ul>
  //  */
  // val all = new MathFuncConfig(
  //   funcs = Seq(Sqrt, InvSqrt, Reciprocal, Sin, Cos,
  //     ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2, Exp, Log,
  //     Sigmoid, SoftPlus, ScaleMixtureGaussian),
  //   nOrder = 2,
  //   addrW = 8,
  //   inBitW = RealSpec.Float32Spec.manW - 8,
  //   outBitW = RealSpec.Float32Spec.manW + 3,
  //   scaleMixtureGaussianSigma = Some((exp(-1.0), exp(-6.0)))
  // )
}

/** Config class to set pipeline stages in a [[rial.math.MathFunctions]] module.
 *
 * @constructor create a new MathFuncConfig.
 * @param preStage     pipeline stages of preprocess.
 * @param calcStage    pipeline stages of table/polynomial and non-table path.
 * @param postStage    pipeline stages of postprocess.
 *
 */
class MathFuncPipelineConfig(
  val preStage:     Int,
  val tableStage:   Int,
  val calcStage:    Int,
  val postStage:    Int,
  val preMulStage:  Int = 0,
  val postMulStage: Int = 0,
) {

  assert(preMulStage <= preStage,
    "preStage includes preMulStage, so should be larger than preMulStage.")
  assert(postMulStage <= postStage,
    "postStage includes postMulStage, so should be larger than postMulStage.")

  /** Calculates the total latency in clock cycles.
   */
  def total = {
             1 +
    preStage   +
    tableStage +
    calcStage  +
    postStage
  }

  val nonTableStage = tableStage + calcStage

  /** Generates a string that represents the current config.
   */
  def getString = {
    f"pre:${preStage}(mul=${preMulStage}), " +
    f"table:${tableStage}, " +
    f"calc:${calcStage}, " +
    f"post:${postStage}(mul=${postMulStage})"
  }
}

/** Factory for [[rial.math.MathFuncPipelineConfig]].
 */
object MathFuncPipelineConfig {
  /** constructs a config corresponding to the single cycle mathfunc.
   */
  def none = {
    new MathFuncPipelineConfig(0,0,0,0,0,0)
  }
}

