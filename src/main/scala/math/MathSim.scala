//% @file MathSim.scala
//
// Stable software-simulation facade for the JNI binding (jni/, rialcpp).
//
// The C++ side looks up exactly this class and these method names, so the
// binding does not depend on the internal per-function Sim/Table object names
// (which change as math/ is refactored).
//
// Unlike a plain `sqrt(x)` export, the approximation depends on the floating
// point format and on the table/polynomial parameters, so those are given to
// the constructor: one MathSim instance == one configuration.  The polynomial
// tables are generated once here and reused by every call, and several
// instances with different configurations can coexist.
//
// The constructor takes only `Int`s so that the C++ side can build an instance
// through plain JNI (`NewObject` with a "(IIIIIII)V" signature) without having
// to construct Scala objects such as RealSpec / MathFuncConfig.
//
package rial.math

import rial.arith.{RealGeneric, RealSpec}
import rial.table.PolynomialCoefficients

/** Software simulation of the [[rial.math.MathFunctions]] module.
 *
 *  Covers [[MathFuncConfig.standard]], i.e. every function the C++ binding
 *  exposes.  The tables are shared exactly as the hardware shares them (one
 *  merged set of bit widths over all the functions), so a MathSim built from
 *  the same numbers as a `MathFunctions` module reproduces it bit-exactly.
 *
 *  @constructor build the polynomial tables for the given configuration.
 *  @param exW     exponent width of the floating point format.
 *  @param exBias  exponent bias of the floating point format.
 *  @param manW    mantissa width of the floating point format.
 *  @param nOrder  order of the approximating polynomial.
 *  @param addrW   width of the table address (== log2 of the interval count).
 *  @param inBitW  width of the in-interval polynomial input.
 *  @param outBitW width of the polynomial output ((1 << outBitW) means 1.0).
 */
class MathSim(
  val exW:     Int,
  val exBias:  Int,
  val manW:    Int,
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) {
  import FuncKind._

  val spec = new RealSpec(exW, exBias, manW)

  /** The underlying config; the same one a [[rial.math.MathFunctions]] module
   *  would be built from. */
  val config = new MathFuncConfig(
    funcs   = MathFuncConfig.standard,
    nOrder  = nOrder,
    addrW   = addrW,
    inBitW  = inBitW,
    outBitW = outBitW,
  )

  // --------------------------------------------------------------------------
  // Coefficient tables, resolved once.

  private val cSqrt       = config.coeffs(Sqrt)
  private val cInvSqrt    = config.coeffs(InvSqrt)
  private val cReciprocal = config.coeffs(Reciprocal)
  private val cSin        = config.coeffs(Sin)
  private val cCos        = config.coeffs(Cos)
  private val cACos       = config.coeffs(ACosPhase2)
  private val cATan2      = config.coeffs(ATan2Phase2)
  private val cExp        = config.coeffs(Exp)

  // acos/atan2 are two-phase; phase 1 reuses the sqrt/reciprocal tables.
  private val cACosPhase1  = config.coeffs(Sqrt)
  private val cATan2Phase1 = config.coeffs(Reciprocal)

  // log uses three tables (normal / small-positive / small-negative) that the
  // shared `coeffs` map only carries in concatenated form, so go via the table.
  private val logTable: LogTable = config.tables(Log).asInstanceOf[LogTable]
  private val cLogNormal:   PolynomialCoefficients = logTable.genNormalCoeffs  (config.bitW)
  private val cLogSmallPos: PolynomialCoefficients = logTable.genSmallPosCoeffs(config.bitW)
  private val cLogSmallNeg: PolynomialCoefficients = logTable.genSmallNegCoeffs(config.bitW)

  // --------------------------------------------------------------------------

  private def real(x: Float): RealGeneric = RealGeneric.fromFloat(spec, x)

  def sqrt      (x: Float): Float = SqrtSimGeneric      (cSqrt,       real(x)).toFloat
  def invsqrt   (x: Float): Float = InvSqrtSimGeneric   (cInvSqrt,    real(x)).toFloat
  def reciprocal(x: Float): Float = ReciprocalSimGeneric(cReciprocal, real(x)).toFloat
  def sin       (x: Float): Float = SinCosSimGeneric(true,  cSin, real(x)).toFloat
  def cos       (x: Float): Float = SinCosSimGeneric(false, cCos, real(x)).toFloat
  def exp       (x: Float): Float = ExpSimGeneric       (cExp,        real(x)).toFloat

  def acos(x: Float): Float =
    ACosPhase2SimGeneric(cACos, ACosPhase1SimGeneric(cACosPhase1, real(x))).toFloat

  def atan2(y: Float, x: Float): Float =
    ATan2Phase2SimGeneric(cATan2,
      ATan2Phase1SimGeneric(cATan2Phase1, real(y), real(x))).toFloat

  def log(x: Float): Float =
    LogSimGeneric(cLogNormal, cLogSmallPos, cLogSmallNeg, real(x)).toFloat
}

/** Default configurations for [[MathSim]]. */
object MathSim {
  /** The FP32 configuration the binding used before it became configurable:
   *  the parameters of [[MathFuncConfig.standardFP32]]. */
  def standardFP32: MathSim = new MathSim(
    exW     = RealSpec.Float32Spec.exW,
    exBias  = RealSpec.Float32Spec.exBias,
    manW    = RealSpec.Float32Spec.manW,
    nOrder  = 2,
    addrW   = 8,
    inBitW  = RealSpec.Float32Spec.manW - 8,
    outBitW = RealSpec.Float32Spec.manW + 3,
  )
}
