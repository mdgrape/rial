package rial.math

import chisel3._
import chisel3.util._

import rial.arith.RealSpec
import rial.table.{PolynomialBitWidth, PolynomialCoefficients, PolynomialGenerator}

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

  /** 2^x, for an argument the caller has already scaled.
   *
   * `pow2(x * log2(e))` is exp(x), up to the rounding of that product.
   */
  val Pow2        = Value
  /** sin(π x), for an argument the caller has already scaled.
   *
   * `sinPi(x / π)` is sin(x), up to the rounding of that quotient.
   */
  val SinPi       = Value
  /** cos(π x), for an argument the caller has already scaled.
   *
   * `cosPi(x / π)` is cos(x), up to the rounding of that quotient.
   */
  val CosPi       = Value
  /** sin(π x), for an argument the caller has already scaled, without either
   *  shared multiplier.
   *
   * Its table holds sin(π y) of the reduced argument itself, one region per
   * binade of y, each scaled by a power of two, so its precision is relative:
   * the result stays within a few units of its last place wherever it is
   * read, the arguments at which the sine vanishes included.
   */
  val SinPiFixed  = Value
  /** cos(π x), for an argument the caller has already scaled, without either
   *  shared multiplier.
   *
   * Its table holds sin(π y) of the reduced argument itself, one region per
   * binade of y, each scaled by a power of two, so its precision is relative:
   * the result stays within a few units of its last place wherever it is
   * read, the arguments at which the cosine vanishes included.
   */
  val CosPiFixed  = Value
  /** log2(x), without either shared multiplier.
   *
   * Its table holds |log2(x)| itself, one region per binade of the distance
   * of x from one on either side of it, each scaled by a power of two, so its
   * precision is relative: the result stays within a few units of its last
   * place wherever it is read, the arguments that approach one included.
   */
  val Log2Fixed   = Value
  /** atan(t) for |t| <= 1, without either shared multiplier.
   *
   * Its table approximates atan(t) itself, so its precision is absolute: the
   * relative error grows as the result approaches zero, until t is small
   * enough for atan(t) to be t, which it returns exactly. It takes the ratio
   * itself, so the caller does the min/max, the reciprocal and the quadrant
   * of an atan2; for |t| > 1 the result is finite but meaningless.
   */
  val ATanFixed   = Value
  /** acos(x) for |x| <= 1, in one pass, without either shared multiplier.
   *
   * Its table holds the angle itself, one region per binade of 1 - |x|, each
   * scaled by a power of two, so its precision is relative: the result stays
   * within a few units of its last place over the whole of [-1, 1], the end
   * where the angle vanishes included. A finite argument beyond an end of
   * that range gets the value of that end, so an argument that overshoots it
   * by rounding needs no clamp; an infinity or a NaN gives NaN.
   */
  val ACosFixed   = Value

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
      case Pow2        => "Pow2"
      case SinPi       => "SinPi"
      case CosPi       => "CosPi"
      case SinPiFixed  => "SinPiFixed"
      case CosPiFixed  => "CosPiFixed"
      case Log2Fixed   => "Log2Fixed"
      case ATanFixed   => "ATanFixed"
      case ACosFixed   => "ACosFixed"
    }
  }

  /** Returns the name of `fn` in the pipeline stage points of
   *  [[rial.math.MathFunctions]] (e.g. FuncKind.InvSqrt -> "invSqrt").
   *
   * The pre- and post-process of a function are one path of the parallels
   * `"pre"` and `"post"`, named after the function: the boundaries of the
   * pre-process of `fn` are `"pre.<pointName(fn)>.*"`. `Sin` and `Cos` share
   * one pre-process and one post-process, so both give `"sinCos"`, and
   * `SinPi` and `CosPi` likewise both give `"sinPi"`, `SinPiFixed` and
   * `CosPiFixed` both `"sinPiFixed"`.
   *
   * @param fn An enumerator of the function.
   */
  def pointName(fn: FuncKind): String = {
    fn match {
      case Sqrt        => "sqrt"
      case InvSqrt     => "invSqrt"
      case Reciprocal  => "reciprocal"
      case Sin         => "sinCos"
      case Cos         => "sinCos"
      case ACosPhase1  => "acosPhase1"
      case ACosPhase2  => "acosPhase2"
      case ATan2Phase1 => "atan2Phase1"
      case ATan2Phase2 => "atan2Phase2"
      case Exp         => "exp"
      case Log         => "log"
      case Sigmoid     => "sigmoid"
      case SoftPlus    => "softPlus"
      case ScaleMixtureGaussian => "scaleMixtureGaussian"
      case Pow2        => "pow2"
      case SinPi       => "sinPi"
      case CosPi       => "sinPi"
      case SinPiFixed  => "sinPiFixed"
      case CosPiFixed  => "sinPiFixed"
      case Log2Fixed   => "log2Fixed"
      case ATanFixed   => "atanFixed"
      case ACosFixed   => "acosFixed"
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
    fn == Log || fn == Sin || fn == Cos || fn == ScaleMixtureGaussian ||
    fn == SinPi || fn == CosPi
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
      case ACosPhase1           => Seq(ACosPhase1,  Sqrt)
      case ACosPhase2           => Seq(ACosPhase1,  ACosPhase2, Sqrt)
      case ATan2Phase1          => Seq(ATan2Phase1, ATan2Phase2, Reciprocal)
      case ATan2Phase2          => Seq(ATan2Phase1, ATan2Phase2, Reciprocal)
      case Exp                  => Seq(Exp)
      case Log                  => Seq(Log)
      case Sigmoid              => Seq(Sigmoid)
      case SoftPlus             => Seq(SoftPlus)
      case ScaleMixtureGaussian => Seq(ScaleMixtureGaussian)
      case Pow2                 => Seq(Pow2)
      case SinPi                => Seq(SinPi, CosPi)
      case CosPi                => Seq(SinPi, CosPi)
      case SinPiFixed           => Seq(SinPiFixed, CosPiFixed)
      case CosPiFixed           => Seq(SinPiFixed, CosPiFixed)
      case Log2Fixed            => Seq(Log2Fixed)
      case ATanFixed            => Seq(ATanFixed)
      case ACosFixed            => Seq(ACosFixed)
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

/** The width of the table address of each function of a [[MathFuncConfig]].
 *
 * A table whose address is `n` bits wide holds `2^n` sub-intervals, so the
 * address width is what one function costs in ROM. The functions of a
 * configuration rarely need the same one, so a configuration takes a width
 * per function rather than a single number.
 *
 * A function that reads the table of another one reads it as it is, so a
 * width belongs to the function that owns the table; [[TableAddrWidth.owner]]
 * says which function that is.
 *
 * An `Int` converts to [[UniformAddrWidth]], so one width for every table is
 * written as that number alone.
 *
 * {{{
 * import rial.math.FuncKind._
 *
 * new MathFuncConfig(
 *   funcs   = MathFuncConfig.standard,
 *   spec    = RealSpec.Float32Spec,
 *   nOrder  = 2,
 *   addrW   = PerFuncAddrWidth(8, Map(Sqrt -> 5, InvSqrt -> 6, Exp -> 6)),
 *   inBitW  = RealSpec.Float32Spec.manW - 8,
 *   outBitW = RealSpec.Float32Spec.manW + 3,
 * )
 * }}}
 */
sealed trait TableAddrWidth {
  /** The width of the address of the table `fn` reads, in bits.
   *
   *  @param fn An enumerator of the function.
   */
  def of(fn: FuncKind.FuncKind): Int
}

/** One address width for every table of a configuration.
 *
 *  @param addrW the width of every table address, in bits
 */
case class UniformAddrWidth(addrW: Int) extends TableAddrWidth {
  def of(fn: FuncKind.FuncKind): Int = { addrW }
}

/** An address width of its own for the tables that ask for one.
 *
 * A function `widths` does not name reads its table at `default`. Every key
 * has to own its table, as [[TableAddrWidth.owner]] says: a function that
 * reads another one's table cannot be given a width of its own.
 *
 *  @param default the width of the address of every table `widths` leaves out
 *  @param widths  the width of the address of the table of each function named
 */
case class PerFuncAddrWidth(
  default: Int,
  widths:  Map[FuncKind.FuncKind, Int],
) extends TableAddrWidth {

  widths.keys.foreach(fn => {
    val own = TableAddrWidth.owner(fn)
    require(own == fn,
      f"${FuncKind.getString(fn)} reads the table of " +
      f"${FuncKind.getString(own)}, so it has no address width of its own: " +
      f"give the width to ${FuncKind.getString(own)} instead")
  })

  def of(fn: FuncKind.FuncKind): Int = {
    widths.getOrElse(TableAddrWidth.owner(fn), default)
  }
}

/** Factory for [[TableAddrWidth]].
 */
object TableAddrWidth {
  import FuncKind._
  import scala.language.implicitConversions

  /** Returns the function that owns the table `fn` reads.
   *
   * `Cos` reads the table of `Sin`, `CosPi` that of `SinPi`, `CosPiFixed`
   * that of `SinPiFixed`, the first phase of `acos` that of `Sqrt` and the
   * first phase of `atan2` that of `Reciprocal`. Every other function owns
   * its table and is its own owner.
   *
   * @param fn An enumerator of the function.
   */
  def owner(fn: FuncKind): FuncKind = {
    fn match {
      case Cos         => Sin
      case CosPi       => SinPi
      case CosPiFixed  => SinPiFixed
      case ACosPhase1  => Sqrt
      case ATan2Phase1 => Reciprocal
      case other       => other
    }
  }

  /** One width, for every table of the configuration alike. */
  implicit def uniform(addrW: Int): TableAddrWidth = { UniformAddrWidth(addrW) }
}

// ============================================================================

/** Where the table of each function of a [[MathFuncConfig]] sits in the ROM
 *  they share.
 *
 * [[rial.math.MathFunctions]] holds the coefficients of every function it is
 * configured for in one ROM, one table after another. The address a function
 * reads it at is the address of its own table under the row that table starts
 * at, which [[offset]] gives; [[depth]] is the rows of that table, [[rows]]
 * the rows of the whole ROM and [[packed]] what it holds, row by row, each
 * packed as [[rial.table.PolynomialGenerator.packTable]] packs it.
 *
 * A function that reads the table of another one, as [[TableAddrWidth.owner]]
 * says, reads it where its owner put it, so the two give the same offset.
 *
 * {{{
 * import rial.math.FuncKind._
 *
 * val layout = MathFuncConfig.standardFP32.romLayout
 *
 * layout.rows          // the rows of the merged ROM
 * layout.offset(Sqrt)  // the row the table of sqrt(x) starts at
 * layout.depth(Sqrt)   // the rows of that table
 * }}}
 *
 *  @param funcs  the functions whose tables the ROM holds
 *  @param coeffs the coefficients of the table of each of them
 */
class RomLayout(
  val funcs:  Seq[FuncKind.FuncKind],
  val coeffs: Map[FuncKind.FuncKind, PolynomialCoefficients],
) {
  import FuncKind.FuncKind

  /** The functions that own a table, in the order their tables are placed:
   *  the deepest table first, so that every table starts at a multiple of its
   *  own depth and an address is the offset with the address of the table
   *  itself below it. */
  val owners: Seq[FuncKind] = funcs.map(fn => TableAddrWidth.owner(fn))
    .distinct.sortBy(fn => (-coeffs(fn).coeffs.length, fn.id))

  private val offsets: Map[FuncKind, Int] = owners.zip(
    owners.scanLeft(0)((row, fn) => row + coeffs(fn).coeffs.length)).toMap

  /** the rows of the merged ROM. */
  val rows: Int = owners.map(fn => coeffs(fn).coeffs.length).sum

  /** the width of the address of the merged ROM, in bits. */
  val addrW: Int = log2Ceil(rows)

  /** the width of one row of the merged ROM, in bits. */
  val rowW: Int = coeffs(owners.head).bitW.coefBits.reduce(_+_)

  /** the rows of the merged ROM, in the order the tables are placed. */
  val packed: Seq[BigInt] = owners.flatMap(fn =>
    PolynomialGenerator.packTable(coeffs(fn)))

  /** Returns the row of the merged ROM the table `fn` reads starts at.
   *
   *  @param fn An enumerator of the function.
   */
  def offset(fn: FuncKind): Int = { offsets(TableAddrWidth.owner(fn)) }

  /** Returns the rows of the table `fn` reads.
   *
   *  @param fn An enumerator of the function.
   */
  def depth(fn: FuncKind): Int = {
    coeffs(TableAddrWidth.owner(fn)).coeffs.length
  }

  /** Returns the width of the address of the table `fn` reads, in bits: what
   *  the argument reduction of `fn` emits, below the offset of that table.
   *
   *  @param fn An enumerator of the function.
   */
  def adrW(fn: FuncKind): Int = { log2Ceil(depth(fn)) }

  owners.foreach(fn => {
    require(isPow2(depth(fn)),
      f"the table of ${FuncKind.getString(fn)} holds ${depth(fn)}%d rows, " +
      f"which is not a power of two")
    require(offset(fn) % depth(fn) == 0,
      f"the table of ${FuncKind.getString(fn)} starts at row " +
      f"${offset(fn)}%d, which its ${depth(fn)}%d rows do not divide")
    require(coeffs(fn).bitW.coefBits == coeffs(owners.head).bitW.coefBits,
      f"the table of ${FuncKind.getString(fn)} does not hold the rows of the " +
      f"same width as the others")
  })
}

// ============================================================================

/** A Config class for [[rial.math.MathFunctions]] Module.

 *  @constructor create a new MathFuncConfig.
 *  @param funcs the list of functions that should be supported.
 *  @param spec the floating point format the tables are generated for. A
 *              table that gives a region of its own to each binade of its
 *              argument places its grid on the mantissa of that format, so it
 *              has to know how many bits of it there are.
 *  @param addrW the width of the address of each table; a plain number gives
 *               every table the same one. A function that reads the table of
 *               another one reads it at the width of that function, so only
 *               the function that owns a table can be given a width; see
 *               [[TableAddrWidth]].
 *  @param inBitW the width of the position inside a sub-interval, shared by
 *                every table. A table whose address is narrower leaves more
 *                bits of the fraction below it than that width holds, and the
 *                ones it does not hold are dropped, so `inBitW` is what sets
 *                the precision of the position of every table alike.
 *  @param ScaleMixtureGaussianSigma sigmaA and sigmaB for ScaleMixtureGaussian.
 *  @param sinCos the range of the argument sin/cos are calculated for; `None`
 *                covers the whole format. Outside the range sin/cos return a
 *                meaningless, but finite, value in [-1, 1]. It restricts
 *                sinPi/cosPi to the arguments of the same sin/cos, i.e. to the
 *                x whose π x lies in that range.
 */
class MathFuncConfig(
  val funcs:   Seq[FuncKind.FuncKind],
  val spec:    RealSpec,
  val nOrder:  Int,
  val addrW:   TableAddrWidth,
  val inBitW:  Int,
  val outBitW: Int,
  val scaleMixtureGaussianSigma: Option[(Double, Double)] = None,
  val sinCos:  Option[SinCosDomain] = None,
) {

  /** mantissa width of the format the tables are generated for. */
  val manW = spec.manW

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

  // Two functions of one of these pairs read the same table, which the
  // configuration would then hold twice.
  Seq((Exp, Pow2), (Sin, SinPi), (Cos, CosPi)).foreach(pair => {
    val (lhs, rhs) = pair
    require(!(has(lhs) && has(rhs)),
      f"${FuncKind.getString(lhs)} and ${FuncKind.getString(rhs)} read the " +
      f"same table: configure one of them")
  })

  /** The width of UInt to represent function select signal.
   *
   *  signal contains: `None ++ [funcs..] ++ Invalid`.
   *  `io.in.bits.sel` of each functions corresponds to non-zero index in `funcs`.
   *  If `None` is given, it calculates nothing.
   *  If `Invalid` is given, it calculates something invalid.
   *
   *  `Invalid` is introduced to make implementation easy. There is no sense in
   *  sending `Invalid` as `sel` signal.
   */
  val signalW = log2Ceil(2 + funcs.length)

  /** Returns function select signal to be passed to [[rial.math.MathFunctions]].
   *
   * {{{
   * val mathfunc = Module(new MathFunction(...))
   *
   * when(doSqrt) {
   *   mathfunc.io.in.bits.sel := MathFuncConfig.signal(FuncKind.Sqrt)
   * }.otherwise {
   *   mathfunc.io.in.bits.sel := MathFuncConfig.signalNone
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
    case Sqrt        => (Sqrt        -> new SqrtTable       (nOrder, addrW.of(Sqrt),        inBitW, outBitW))
    case InvSqrt     => (InvSqrt     -> new InvSqrtTable    (nOrder, addrW.of(InvSqrt),     inBitW, outBitW))
    case Reciprocal  => (Reciprocal  -> new ReciprocalTable (nOrder, addrW.of(Reciprocal),  inBitW, outBitW))
    case Sin         => (Sin         -> new SinCosTable     (nOrder, addrW.of(Sin),         inBitW, outBitW))
    case Cos         => (Cos         -> new SinCosTable     (nOrder, addrW.of(Cos),         inBitW, outBitW))
    case Exp         => (Exp         -> new ExpTable        (nOrder, addrW.of(Exp),         inBitW, outBitW))
    case Log         => (Log         -> new LogTable        (nOrder, addrW.of(Log),         inBitW, outBitW))
    case ACosPhase1  => (ACosPhase1  -> new SqrtTable       (nOrder, addrW.of(ACosPhase1),  inBitW, outBitW))
    case ACosPhase2  => (ACosPhase2  -> new ACosPhase2Table (nOrder, addrW.of(ACosPhase2),  inBitW, outBitW))
    case ATan2Phase1 => (ATan2Phase1 -> new ReciprocalTable (nOrder, addrW.of(ATan2Phase1), inBitW, outBitW))
    case ATan2Phase2 => (ATan2Phase2 -> new ATan2Phase2Table(nOrder, addrW.of(ATan2Phase2), inBitW, outBitW))
    case Sigmoid     => (Sigmoid     -> new SigmoidTable    (nOrder, addrW.of(Sigmoid),     inBitW, outBitW))
    case SoftPlus    => (SoftPlus    -> new SoftPlusTable   (nOrder, addrW.of(SoftPlus),    inBitW, outBitW))
    case Pow2        => (Pow2        -> new ExpTable        (nOrder, addrW.of(Pow2),        inBitW, outBitW))
    case SinPi       => (SinPi       -> new SinCosTable     (nOrder, addrW.of(SinPi),       inBitW, outBitW))
    case CosPi       => (CosPi       -> new SinCosTable     (nOrder, addrW.of(CosPi),       inBitW, outBitW))
    case SinPiFixed  => (SinPiFixed  -> new SinPiFixedTable (nOrder, addrW.of(SinPiFixed),  inBitW, outBitW, manW))
    case CosPiFixed  => (CosPiFixed  -> new SinPiFixedTable (nOrder, addrW.of(CosPiFixed),  inBitW, outBitW, manW))
    case Log2Fixed   => (Log2Fixed   -> new Log2FixedTable  (nOrder, addrW.of(Log2Fixed),   inBitW, outBitW, manW))
    case ATanFixed   => (ATanFixed   -> new ATanFixedTable  (nOrder, addrW.of(ATanFixed),   inBitW, outBitW, manW))
    case ACosFixed   => (ACosFixed   -> new ACosFixedTable  (nOrder, addrW.of(ACosFixed),   inBitW, outBitW))
    case ScaleMixtureGaussian => {
      val (sA, sB) = scaleMixtureGaussianSigma.get
      (ScaleMixtureGaussian -> new ScaleMixtureGaussianTable(
        sA, sB, nOrder, addrW.of(ScaleMixtureGaussian), inBitW, outBitW))
    }
  }.toMap

  /** the widths the table of each function asks for, before they are merged. */
  private val tableBitWs = tables.map{ case (k,v) => (k -> v.getBitWidth()) }.toMap

  /** The widths the shared polynomial evaluator is built for: the widest
   *  coefficient and accumulator of any table of this configuration. Its
   *  `addrW` is the widest table address; the address of one table is
   *  [[bitWOf]]. */
  val bitW = tableBitWs.values.reduce((l, r) => { l.merge(r) } )

  /** The coefficients of the table of each function, at the widths the shared
   *  polynomial evaluator reads them at and at the address width of that
   *  table. */
  val coeffs = tables.map{ case (k,v) =>
    (k -> v.generateTable(bitW.copy(addrW = tableBitWs(k).addrW)))
  }.toMap

  // The exp table approximates 2^t-1 on [0, 1) with `padding` bits below the
  // mantissa LSB, and ExpPostProcess rounds that result to manW bits.  Should
  // the rounding carry into the hidden bit, the correction term in
  // ExpPostProcess would be built from the wrong mantissa, so require that the
  // table cannot reach the value at which the rounding carries.
  funcs.filter(fn => fn == Exp || fn == Pow2).foreach(fn => {
    val coeff   = coeffs(fn)
    val padding = coeff.bitW.outBitW - manW
    if (padding > 0) {
      val maxPoly   = ExpTable.maxPolynomial(coeff, spec)
      val threshold = (BigInt(1) << coeff.bitW.outBitW) -
                      (BigInt(1) << (padding - 1))

      require(maxPoly < threshold,
        f"ExpTable: the polynomial reaches $maxPoly, at or above the " +
        f"$threshold at which rounding to ${manW} bits carries into the " +
        f"hidden bit (nOrder=${coeff.bitW.nOrder}, " +
        f"addrW=${coeff.bitW.addrW}, inBitW=${coeff.bitW.inBitW}, " +
        f"outBitW=${coeff.bitW.outBitW}, manW=${manW}, padding=$padding). " +
        "Widen outBitW or refine the table so that the approximation of " +
        "2^t-1 stays below the threshold.")
    }
  })

  /** Returns the widths the table of `fn` is read at: those of [[bitW]], with
   *  the address width of that table.
   *
   *  @param fn An enumerator of the function.
   */
  def bitWOf(fn: FuncKind): PolynomialBitWidth = {
    coeffs(fn).bitW
  }

  /** Where the table of each function sits in the ROM they share; see
   *  [[RomLayout]]. */
  val romLayout = new RomLayout(funcs, coeffs)

  // --------------------------------------------------------------------------

  /** The same configuration with the argument range of sin/cos restricted.
   *
   *  @param domain the range to restrict the reduction to, `None` for the
   *                whole format
   */
  def restrictSinCos(domain: Option[SinCosDomain]): MathFuncConfig = {
    new MathFuncConfig(funcs, spec, nOrder, addrW, inBitW, outBitW,
      scaleMixtureGaussianSigma, domain)
  }

  def getString: String = {
    // Compact representation when funcs matches one of the predefined sets in
    // the companion object.  Keeps the test names (which include this string)
    // short enough to avoid hitting OS filename-length limits.
    val fs =
      if      (funcs == MathFuncConfig.all           ) { "all" }
      else if (funcs == MathFuncConfig.standard      ) { "standard" }
      else if (funcs == MathFuncConfig.standardScaled) { "standardScaled" }
      else if (funcs == MathFuncConfig.standardFixed ) { "standardFixed" }
      else if (funcs == MathFuncConfig.simple        ) { "simple" }
      else    { funcs.map(fn => FuncKind.getString(fn)).mkString("[", ", ", "]") }
    fs + sinCos.map(d => f"_sinCosEx${d.maxEx}").getOrElse("")
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
  val simple = Seq(Sqrt, InvSqrt, Reciprocal)

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
  val standard = Seq(Sqrt, InvSqrt, Reciprocal, Sin, Cos,
    ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2, Exp, Log)

  /** The [[standard]] set with the functions that need the pre-process
   * multiplier replaced by the ones that take a scaled argument.
   *
   * includes:
   * <ul>
   *   <li> sqrt       </li>
   *   <li> invsqrt    </li>
   *   <li> reciprocal </li>
   *   <li> acos       </li>
   *   <li> atan2      </li>
   *   <li> log        </li>
   *   <li> pow2       </li>
   *   <li> sinpi      </li>
   *   <li> cospi      </li>
   * </ul>
   */
  val standardScaled = Seq(Sqrt, InvSqrt, Reciprocal,
    ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2, Log, Pow2, SinPi, CosPi)

  /** The [[standard]] set with every function that needs a shared multiplier
   * replaced by the one that takes a scaled argument or reads the result
   * itself out of its table.
   */
  val standardFixed = normalize(Seq(Sqrt, InvSqrt, Reciprocal, Pow2,
    Log2Fixed, SinPiFixed, CosPiFixed, ACosFixed, ATanFixed))

  /** Defines all the supported math functions including too task-specific ones.
   *
   * includes:
   * <ul>
   *   <li> everything. </li>
   * </ul>
   */
  val all = Seq(Sqrt, InvSqrt, Reciprocal, Sin, Cos,
    ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2, Exp, Log,
    Sigmoid, SoftPlus, ScaleMixtureGaussian)

  /** Default (sigmaA, sigmaB) for [[FuncKind.ScaleMixtureGaussian]]. */
  val defaultSMGSigma: (Double, Double) = (exp(-1.0), exp(-6.0))

  // --------------------------------------------------------------------------
  // Default MathFuncConfig instances.  The HW table parameters (nOrder, addrW,
  // outBitW) match the typical per-function setups used elsewhere; users
  // wanting tighter precision should construct a custom MathFuncConfig.

  /** Default FP32 config supporting only the [[simple]] function set. */
  val simpleFP32 = new MathFuncConfig(
    funcs   = simple,
    spec    = RealSpec.Float32Spec,
    nOrder  = 2,
    addrW   = 8,
    inBitW  = RealSpec.Float32Spec.manW - 8,
    outBitW = RealSpec.Float32Spec.manW + 3,
  )
  /** Default BF16 config supporting only the [[simple]] function set. */
  val simpleBF16 = new MathFuncConfig(
    funcs   = simple,
    spec    = RealSpec.BFloat16Spec,
    nOrder  = 0,
    addrW   = RealSpec.BFloat16Spec.manW,
    inBitW  = 0,
    outBitW = RealSpec.BFloat16Spec.manW + 1,
  )

  /** Default FP32 config supporting the [[standard]] function set. */
  val standardFP32 = new MathFuncConfig(
    funcs   = standard,
    spec    = RealSpec.Float32Spec,
    nOrder  = 2,
    addrW   = 8,
    inBitW  = RealSpec.Float32Spec.manW - 8,
    outBitW = RealSpec.Float32Spec.manW + 3,
  )
  /** Default BF16 config supporting the [[standard]] function set. */
  val standardBF16 = new MathFuncConfig(
    funcs   = standard,
    spec    = RealSpec.BFloat16Spec,
    nOrder  = 0,
    addrW   = RealSpec.BFloat16Spec.manW,
    inBitW  = 0,
    outBitW = RealSpec.BFloat16Spec.manW + 1,
  )

  /** Default FP32 config supporting the [[standardScaled]] function set. */
  val standardScaledFP32 = new MathFuncConfig(
    funcs   = standardScaled,
    spec    = RealSpec.Float32Spec,
    nOrder  = 2,
    addrW   = 8,
    inBitW  = RealSpec.Float32Spec.manW - 8,
    outBitW = RealSpec.Float32Spec.manW + 3,
  )
  /** Default BF16 config supporting the [[standardScaled]] function set. */
  val standardScaledBF16 = new MathFuncConfig(
    funcs   = standardScaled,
    spec    = RealSpec.BFloat16Spec,
    nOrder  = 0,
    addrW   = RealSpec.BFloat16Spec.manW,
    inBitW  = 0,
    outBitW = RealSpec.BFloat16Spec.manW + 1,
  )

  /** Default FP32 config supporting the [[standardFixed]] function set.
   *
   *  The position inside a sub-interval is wider here than the mantissa
   *  leaves below the address: the tables of this set give a region of their
   *  own to each binade of the argument and divide that region into
   *  sub-intervals, so a position reaches further down the mantissa than a
   *  flat table's does. A fraction with fewer bits left is padded, so the
   *  functions whose table is flat read exactly what a narrower position
   *  would have given them.
   */
  val standardFixedFP32 = new MathFuncConfig(
    funcs   = standardFixed,
    spec    = RealSpec.Float32Spec,
    nOrder  = 2,
    addrW   = 8,
    inBitW  = 18,
    outBitW = RealSpec.Float32Spec.manW + 3,
  )
  /** Default BF16 config supporting the [[standardFixed]] function set. */
  val standardFixedBF16 = new MathFuncConfig(
    funcs   = standardFixed,
    spec    = RealSpec.BFloat16Spec,
    nOrder  = 0,
    addrW   = RealSpec.BFloat16Spec.manW,
    inBitW  = 0,
    outBitW = RealSpec.BFloat16Spec.manW + 1,
  )

  /** Default FP32 config supporting [[all]] math functions.  Uses
   *  `defaultSMGSigma` for ScaleMixtureGaussian.
   *
   *  `outBitW = manW + 8` is driven by [[FuncKind.ScaleMixtureGaussian]]:
   *  with the default sigmas SMG reconstructs the result by multiplying the
   *  polynomial by ~`2^15`, so fewer extra bits noticeably degrade the SMG
   *  output (at +3 the SMG result already carries ~8000 ULP of error).  The
   *  simpler functions (Sqrt / Sin / ...) work at this larger size but their
   *  ROMs become wider than strictly necessary -- use [[standardFP32]] (or a
   *  hand-rolled `MathFuncConfig`) when SMG is not needed.
   */
  val allFP32 = new MathFuncConfig(
    funcs   = all,
    spec    = RealSpec.Float32Spec,
    nOrder  = 2,
    addrW   = 8,
    inBitW  = RealSpec.Float32Spec.manW - 8,
    outBitW = RealSpec.Float32Spec.manW + 8,
    scaleMixtureGaussianSigma = Some(defaultSMGSigma),
  )
  /** Default BF16 config supporting [[all]] math functions.  Uses
   *  `defaultSMGSigma` for ScaleMixtureGaussian.
   *
   *  `outBitW = manW + 10` is driven by [[FuncKind.SoftPlus]] and
   *  [[FuncKind.ScaleMixtureGaussian]]: with the default sigmas SMG
   *  reconstructs the result by multiplying the polynomial by ~`2^15`, so
   *  fewer than ~10 extra bits collapses the SMG output entirely.  SoftPlus
   *  similarly needs ~1 BF16 ULP of headroom because its output range reaches
   *  the table's `rangeMax`.  The simpler functions (Sqrt / Sin / ...) still
   *  work at this larger size but their ROMs become wider than strictly
   *  necessary -- use [[standardBF16]] (or a hand-rolled `MathFuncConfig`)
   *  when SMG/SoftPlus are not needed and the extra precision is wasted area.
   */
  val allBF16 = new MathFuncConfig(
    funcs   = all,
    spec    = RealSpec.BFloat16Spec,
    nOrder  = 0,
    addrW   = RealSpec.BFloat16Spec.manW,
    inBitW  = 0,
    outBitW = RealSpec.BFloat16Spec.manW + 10,
    scaleMixtureGaussianSigma = Some(defaultSMGSigma),
  )
}
