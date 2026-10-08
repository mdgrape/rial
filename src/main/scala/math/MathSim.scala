//% @file MathSim.scala
//
// Stable software-simulation facade for the JNI binding (jni/, rialcpp).
// The C++ side looks up exactly this class and these method names.
//
// Unlike a plain `sqrt(x)` export, the approximation depends on the floating
// point format, on the table/polynomial parameters and on which functions the
// module supports -- the tables share one merged set of bit widths, so a
// sqrt-only module does not produce quite the same sqrt as an all-function one.
// All of that goes to the constructor: one MathSim instance == one
// MathFuncConfig. The tables are generated once here and reused by every call,
// and several instances with different configurations can coexist.
//
// The constructor takes only types plain JNI can build (ints, doubles, an
// int[] and a String[]), so the C++ side can instantiate it with `NewObject`
// without having to construct scala objects such as RealSpec / FuncKind /
// MathFuncConfig.
//
package rial.math

import rial.arith.{RealGeneric, RealSpec}
import rial.table.PolynomialCoefficients

import spire.math.SafeLong

/** Software simulation of the [[rial.math.MathFunctions]] module.
 *
 *  @constructor build the polynomial tables for the given configuration.
 *  @param exW       exponent width of the floating point format.
 *  @param exBias    exponent bias of the floating point format.
 *  @param manW      mantissa width of the floating point format.
 *  @param nOrder    order of the approximating polynomial.
 *  @param addrW     width of the table address (== log2 of the interval
 *                   count) of every table `addrWs` leaves out.
 *  @param addrWs    width of the table address of each function of
 *                   `funcNames`, in that order; empty reads every table at
 *                   `addrW`.  A name that reads the table of another function
 *                   (`cos` reads the one of `sin`) carries that table's width,
 *                   so the two have to agree.
 *  @param inBitW    width of the in-interval polynomial input.
 *  @param outBitW   width of the polynomial output ((1 << outBitW) means 1.0).
 *  @param funcNames the functions to support, named as in [[MathSim.funcNames]].
 *                   Dependencies are pulled in automatically (`acos` brings
 *                   `sqrt`, ...), exactly as [[FuncKind.normalize]] does it, so
 *                   the resulting table widths match a `MathFunctions` module
 *                   built from the same list.  [[supportedFuncs]] reports what
 *                   the instance ended up supporting.
 *  @param smgSigmaA sigmaA for ScaleMixtureGaussian; NaN by default.
 *  @param smgSigmaB sigmaB for ScaleMixtureGaussian; NaN by default.
 */
class MathSim(
  val exW:       Int,
  val exBias:    Int,
  val manW:      Int,
  val nOrder:    Int,
  val addrW:     Int,
  val addrWs:    Array[Int],
  val inBitW:    Int,
  val outBitW:   Int,
  val funcNames: Array[String],
  val smgSigmaA: Double,
  val smgSigmaB: Double,
) {
  import FuncKind._

  val spec = new RealSpec(exW, exBias, manW)

  private val kinds = MathSim.decodeFuncNames(funcNames.toIndexedSeq)

  private val funcs = FuncKind.normalize(kinds)

  /** The width each table is read at: the one `addrWs` gives for the function
   *  that owns it, and `addrW` for a table no name covers. */
  private val tableAddrW: TableAddrWidth = if (addrWs.isEmpty) {
    UniformAddrWidth(addrW)
  } else {
    require(addrWs.length == funcNames.length,
      s"MathSim: ${addrWs.length} table address widths for " +
      s"${funcNames.length} functions")

    val byOwner = kinds.zip(addrWs.toIndexedSeq).groupBy {
      case (fn, _) => TableAddrWidth.owner(fn)
    }
    byOwner.foreach { case (own, ws) =>
      val given = ws.map { case (fn, w) => s"${FuncKind.getString(fn)} -> $w" }
      require(ws.map(_._2).distinct.length == 1,
        s"MathSim: ${given.mkString(", ")} read the table of " +
        s"${FuncKind.getString(own)}, so they need one width")
    }
    PerFuncAddrWidth(addrW, byOwner.map { case (own, ws) => own -> ws.head._2 })
  }

  /** The underlying config; the same one a [[rial.math.MathFunctions]] module
   *  would be built from. */
  val config = new MathFuncConfig(
    funcs   = funcs,
    spec    = spec,
    nOrder  = nOrder,
    addrW   = tableAddrW,
    inBitW  = inBitW,
    outBitW = outBitW,
    scaleMixtureGaussianSigma = if (funcs.contains(ScaleMixtureGaussian)) {
      Some((
        if (smgSigmaA.isNaN) {MathFuncConfig.defaultSMGSigma._1} else {smgSigmaA},
        if (smgSigmaB.isNaN) {MathFuncConfig.defaultSMGSigma._2} else {smgSigmaB},
      ))
    } else { None },
  )

  /** The functions this instance actually supports, named as in
   *  [[MathSim.funcNames]].  A superset of what was requested: the dependencies
   *  a function pulls in become callable too, just as they are selectable in
   *  the corresponding hardware module. */
  def supportedFuncs: Array[String] = {
    funcs.flatMap(fn => MathSim.kindNames.get(fn)).distinct.toArray
  }

  // --------------------------------------------------------------------------
  // Coefficient tables.  `lazy` so that a function left out of `funcs` only
  // complains when it is actually called.

  private def coeffs(fn: FuncKind, name: String): PolynomialCoefficients =
    config.coeffs.getOrElse(fn, throw new IllegalArgumentException(
      s"MathSim: ${name} is not supported by this configuration " +
      s"(supported: ${supportedFuncs.mkString(", ")})"))

  private def table[T](fn: FuncKind, name: String): T =
    config.tables.getOrElse(fn, throw new IllegalArgumentException(
      s"MathSim: ${name} is not supported by this configuration " +
      s"(supported: ${supportedFuncs.mkString(", ")})")).asInstanceOf[T]

  private lazy val cSqrt       = coeffs(Sqrt,        "sqrt")
  private lazy val cInvSqrt    = coeffs(InvSqrt,     "invsqrt")
  private lazy val cReciprocal = coeffs(Reciprocal,  "reciprocal")
  private lazy val cSin        = coeffs(Sin,         "sin")
  private lazy val cCos        = coeffs(Cos,         "cos")
  private lazy val cACos       = coeffs(ACosPhase2,  "acos")
  private lazy val cATan2      = coeffs(ATan2Phase2, "atan2")
  private lazy val cExp        = coeffs(Exp,         "exp")
  private lazy val cSigmoid    = coeffs(Sigmoid,     "sigmoid")
  private lazy val cSMG        = coeffs(ScaleMixtureGaussian, "scaleMixtureGaussian")

  // acos/atan2 are two-phase; phase 1 reuses the sqrt/reciprocal tables.
  private lazy val cACosPhase1  = coeffs(Sqrt,       "acos")
  private lazy val cATan2Phase1 = coeffs(Reciprocal, "atan2")

  // log and softplus split their domain over several sub-tables, which the
  // shared `coeffs` map only carries in concatenated form; go via the table.
  private lazy val logTable = table[LogTable](Log, "log")
  private lazy val cLogNormal   = logTable.genNormalCoeffs  (config.bitWOf(Log))
  private lazy val cLogSmallPos = logTable.genSmallPosCoeffs(config.bitWOf(Log))
  private lazy val cLogSmallNeg = logTable.genSmallNegCoeffs(config.bitWOf(Log))

  private lazy val softPlusTable = table[SoftPlusTable](SoftPlus, "softplus")
  private lazy val cSoftPlusPos  = softPlusTable.genPosCoeffs(config.bitWOf(SoftPlus))
  private lazy val cSoftPlusNeg  = softPlusTable.genNegCoeffs(config.bitWOf(SoftPlus))

  // --------------------------------------------------------------------------

  private def real(x: Float): RealGeneric = { RealGeneric.fromFloat(spec, x) }

  def sqrt      (x: Float): Float = { SqrtSimGeneric      (cSqrt,       real(x)).toFloat }
  def invsqrt   (x: Float): Float = { InvSqrtSimGeneric   (cInvSqrt,    real(x)).toFloat }
  def reciprocal(x: Float): Float = { ReciprocalSimGeneric(cReciprocal, real(x)).toFloat }
  def sin       (x: Float): Float = { SinCosSimGeneric    (true,  cSin, real(x)).toFloat }
  def cos       (x: Float): Float = { SinCosSimGeneric    (false, cCos, real(x)).toFloat }
  def exp       (x: Float): Float = { ExpSimGeneric       (cExp,        real(x)).toFloat }
  def sigmoid   (x: Float): Float = { SigmoidSimGeneric   (cSigmoid,    real(x)).toFloat }

  def acos(x: Float): Float = {
    ACosPhase2SimGeneric(cACos, ACosPhase1SimGeneric(cACosPhase1, real(x))).toFloat
  }

  def atan2(y: Float, x: Float): Float = {
    ATan2Phase2SimGeneric(cATan2,
      ATan2Phase1SimGeneric(cATan2Phase1, real(y), real(x))).toFloat
  }

  def log(x: Float): Float = {
    LogSimGeneric(cLogNormal, cLogSmallPos, cLogSmallNeg, real(x)).toFloat
  }

  def softplus(x: Float): Float = {
    SoftPlusSimGeneric(cSoftPlusPos, cSoftPlusNeg, real(x)).toFloat
  }

  def scaleMixtureGaussian(x: Float): Float = {
    val (sgmA, sgmB) = config.scaleMixtureGaussianSigma.get
    ScaleMixtureGaussianSimGeneric(cSMG, real(x), sgmA, sgmB).toFloat
  }

  // --------------------------------------------------------------------------
  // The same functions on raw bit patterns of `spec`, for formats a Float
  // cannot carry (manW > 23, exW > 8, ...). Bits above `spec.W` are ignored.

  private def realOf(x: Long): RealGeneric = {
    new RealGeneric(spec, SafeLong(BigInt(x) & MathSim.maskOf(spec.W)))
  }
  private def bitsOf(r: RealGeneric): Long = {
    (r.value.toBigInt & MathSim.maskOf(spec.W)).toLong
  }

  def sqrtBits      (x: Long): Long = { bitsOf(SqrtSimGeneric      (cSqrt,       realOf(x))) }
  def invsqrtBits   (x: Long): Long = { bitsOf(InvSqrtSimGeneric   (cInvSqrt,    realOf(x))) }
  def reciprocalBits(x: Long): Long = { bitsOf(ReciprocalSimGeneric(cReciprocal, realOf(x))) }
  def sinBits       (x: Long): Long = { bitsOf(SinCosSimGeneric    (true,  cSin, realOf(x))) }
  def cosBits       (x: Long): Long = { bitsOf(SinCosSimGeneric    (false, cCos, realOf(x))) }
  def expBits       (x: Long): Long = { bitsOf(ExpSimGeneric       (cExp,        realOf(x))) }
  def sigmoidBits   (x: Long): Long = { bitsOf(SigmoidSimGeneric   (cSigmoid,    realOf(x))) }

  def acosBits(x: Long): Long = {
    bitsOf(ACosPhase2SimGeneric(cACos, ACosPhase1SimGeneric(cACosPhase1, realOf(x))))
  }

  def atan2Bits(y: Long, x: Long): Long = {
    bitsOf(ATan2Phase2SimGeneric(cATan2,
      ATan2Phase1SimGeneric(cATan2Phase1, realOf(y), realOf(x))))
  }

  def logBits(x: Long): Long = {
    bitsOf(LogSimGeneric(cLogNormal, cLogSmallPos, cLogSmallNeg, realOf(x)))
  }

  def softplusBits(x: Long): Long = {
    bitsOf(SoftPlusSimGeneric(cSoftPlusPos, cSoftPlusNeg, realOf(x)))
  }

  def scaleMixtureGaussianBits(x: Long): Long = {
    val (sgmA, sgmB) = config.scaleMixtureGaussianSigma.get
    bitsOf(ScaleMixtureGaussianSimGeneric(cSMG, realOf(x), sgmA, sgmB))
  }
}

/** Function names and default configurations for [[MathSim]]. */
object MathSim {
  import FuncKind._

  /** The function names the binding selects functions by.
   *
   *  `acos` / `atan2` name their second phase: the first phase reuses the
   *  sqrt / reciprocal table and [[FuncKind.normalize]] pulls it in.
   */
  val funcNames: Map[String, FuncKind] = Map(
    "sqrt"                 -> Sqrt,
    "invsqrt"              -> InvSqrt,
    "reciprocal"           -> Reciprocal,
    "sin"                  -> Sin,
    "cos"                  -> Cos,
    "acos"                 -> ACosPhase2,
    "atan2"                -> ATan2Phase2,
    "exp"                  -> Exp,
    "log"                  -> Log,
    "sigmoid"              -> Sigmoid,
    "softplus"             -> SoftPlus,
    "scaleMixtureGaussian" -> ScaleMixtureGaussian,
  )

  private val kindNames: Map[FuncKind, String] = funcNames.map(_.swap)

  private def maskOf(w: Int): BigInt = { (BigInt(1) << w) - 1 }

  /** Translates the binding's function names into [[FuncKind]]s.  The result is
   *  not normalized; [[MathSim]] normalizes it. */
  def decodeFuncNames(names: Seq[String]): Seq[FuncKind] = {
    require(names.nonEmpty, "MathSim: at least one function must be selected")
    names.map { n =>
      funcNames.getOrElse(n, throw new IllegalArgumentException(
        s"MathSim: unknown function name '${n}' " +
        s"(known: ${funcNames.keys.toSeq.sorted.mkString(", ")})"))
    }
  }

  /** The configuration of [[MathFuncConfig.standardFP32]]: a second order
   *  polynomial over the [[MathFuncConfig.standard]] function set, each table
   *  read at the narrowest address that keeps its maximum error -- sqrt 5
   *  bits, 1/x and acos 7, the rest 6. */
  def standardFP32: MathSim = { of(MathFuncConfig.standardFP32) }

  /** The parameters of `cfg`, as [[MathSim]] takes them. */
  private def of(cfg: MathFuncConfig): MathSim = {
    val names = cfg.funcs.flatMap(kindNames.get).distinct.toArray
    new MathSim(
      exW       = cfg.spec.exW,
      exBias    = cfg.spec.exBias,
      manW      = cfg.spec.manW,
      nOrder    = cfg.nOrder,
      addrW     = defaultAddrW(cfg.addrW),
      addrWs    = names.map(n => cfg.addrW.of(funcNames(n))),
      inBitW    = cfg.inBitW,
      outBitW   = cfg.outBitW,
      funcNames = names,
      smgSigmaA = cfg.scaleMixtureGaussianSigma.map(_._1).getOrElse(Double.NaN),
      smgSigmaB = cfg.scaleMixtureGaussianSigma.map(_._2).getOrElse(Double.NaN),
    )
  }

  /** The width of a table address the per-function widths leave out. */
  private def defaultAddrW(addrW: TableAddrWidth): Int = {
    addrW match {
      case UniformAddrWidth(w)    => w
      case PerFuncAddrWidth(w, _) => w
    }
  }
}
