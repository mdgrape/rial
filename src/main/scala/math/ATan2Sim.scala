//% @file ATan2Sim.scala
//
// Simulators for atan2(y, x). Two-phase computation:
//
//   phase 1:  (y, x) -> t = min(|x|,|y|) / max(|x|,|y|)    (uses Reciprocal table)
//   phase 2:   t     -> z = atan2(y, x)                    (own table; PostMul)
//
// Phase-1 packs four 2-bit status bits + the original ysgn into the result
// using NaN-boxing (high mantissa bit for special values, low mantissa bit
// for `xsgn ^ yIsLarger`, MSB of exponent for `yIsLarger`). Phase-2 decodes
// these and emits one of the six pi-related constants when needed, or
// otherwise computes one of the four (zsgn, atan(t), pi, pi/2) combinations.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._

import spire.math.SafeLong
import spire.math.Real
import spire.implicits._

import rial.table._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding._

// Phase-2 table: approximates atan(t)/t for t in [0, 1).
//   t -> 0+ : atan(t)/t -> 1, clamped to 1 - eps so the value fits in [0, 1).
//   t -> 1  : atan(1)/1 = pi/4 ~ 0.785.
class ATan2Phase2Table(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  private val func = (t: Double) => {
    val eps = pow(2.0, -outBitW)
    if (t == 0.0) { 1.0 - eps }
    else          { scala.math.min(atan(t) / t, 1.0 - eps) }
  }

  private val nIntervals = 1 << addrW
  private val ft  = FuncTable(func, nOrder, nIntervals, 0.0, 1.0)
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW)

  def getBitWidth(): PolynomialBitWidth = { gen.getBitWidth() }
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    gen.generateTable(bitW)
  }
}

// =========================================================================
// Phase 1: (y, x) -> min/max ratio, with state bits and special-value NaN box
// =========================================================================

object ATan2Phase1SimGeneric {

  /** Phase 1 of atan2.
   *
   *  @param tReciprocal  Reciprocal table coefficients (shared with Reciprocal)
   *  @param y, x         the two atan2 inputs
   *  @return  RealGeneric encoding (zsgn=ysgn, zex with status MSB, zman with
   *           ratio + status LSB), or NaN-boxed special-value flag.
   */
  def apply(tReciprocal: PolynomialCoefficients, y: RealGeneric, x: RealGeneric): RealGeneric = {
    assert(x.spec == y.spec)

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val outPrec   = tReciprocal.bitW.outPrec
    val fracW     = outPrec
    val extraBits = outPrec - manW

    val xsgn = x.sgn
    val ysgn = y.sgn

    val xnan  = x.isNaN
    val xinf  = x.isInfinite
    val xzero = x.isZero
    val ynan  = y.isNaN
    val yinf  = y.isInfinite
    val yzero = y.isZero

    val xySameMan = x.man == y.man
    val yIsLarger = slice(0, spec.W - 1, x.value) < slice(0, spec.W - 1, y.value)

    val minxy = if (yIsLarger) { x } else { y }
    val maxxy = if (yIsLarger) { y } else { x }

    val minex = minxy.ex
    val maxex = maxxy.ex

    val diffexDec = if (yIsLarger) {
      if (y.man > x.man) 1 else 0
    } else {
      if (x.man > y.man) 1 else 0
    }

    val zex0 = if (maxex == maskI(exW) && minex == maskI(exW)) {
      exBias
    } else if ((maxex == maskI(exW) && minex != maskI(exW)) || minex == 0) {
      0
    } else {
      minex - maxex - diffexDec + exBias
    }

    val zEx    = if (zex0 < 0) { 0 } else { zex0 }
    val zeroed = (zEx == 0)

    val tooLargeX = (zeroed && !yIsLarger) || ( xinf && !yinf) || (!xzero &&  yzero)
    val tooLargeY = (zeroed &&  yIsLarger) || (!xinf &&  yinf) || ( xzero && !yzero)

    val atan2Status  = if (yIsLarger) { 2 + x.sgn } else { x.sgn }
    val atan2Status0 = (atan2Status >> 0) & 1
    val atan2Status1 = (atan2Status >> 1) & 1

    // ---- special values (NaN-boxed in mantissa) ----
    val zerozero = xzero && yzero
    val xpos = x.sgn == 0
    val xneg = x.sgn == 1
    val znan       = xnan || ynan
    val zzero      = (tooLargeX && xpos) || (zerozero && xpos)
    val zpi        = (tooLargeX && xneg) || (zerozero && xneg)
    val zhalfpi    = tooLargeY
    val z1piover4  = (xinf && yinf && xpos) || (x.ex == y.ex && x.man == y.man && xpos)
    val z3piover4  = (xinf && yinf && xneg) || (x.ex == y.ex && x.man == y.man && xneg)

    if (znan)      { return new RealGeneric(spec, ysgn, maskI(spec.exW), SafeLong(1) << (spec.manW - 1)) }
    if (zzero)     { return new RealGeneric(spec, ysgn, maskI(spec.exW), SafeLong(1)) }
    if (zpi)       { return new RealGeneric(spec, ysgn, maskI(spec.exW), SafeLong(2)) }
    if (zhalfpi)   { return new RealGeneric(spec, ysgn, maskI(spec.exW), SafeLong(3)) }
    if (z1piover4) { return new RealGeneric(spec, ysgn, maskI(spec.exW), SafeLong(4)) }
    if (z3piover4) { return new RealGeneric(spec, ysgn, maskI(spec.exW), SafeLong(5)) }

    // ---- Reciprocal table lookup at max-mantissa address ----
    // New ReciprocalTable approximates 2/(1+t)-1 for t in [0, 1); adding the
    // hidden 1 gives `denomW1 = (2/(1+max_frac)) * 2^fracW`, the inverse of
    // the max value's mantissa (with the factor of 2 absorbed in the
    // multiplication below).
    val addrW = tReciprocal.bitW.addrW
    val xW    = manW - addrW
    val addr  = slice(xW, addrW, maxxy.man).toInt
    val xiv   = slice(0,  xW,    maxxy.man).toBigInt

    val recMan = if (maxxy.man == SafeLong(0)) {
      SafeLong(0)
    } else {
      SafeLong(PolynomialEval.eval(addr, xiv, tReciprocal))
    }

    // ---- z = min/max ratio ----
    val denomW1 = (SafeLong(1) << fracW) + recMan
    val numerW1 = (SafeLong(1) << manW)  + minxy.man

    val zProd = denomW1 * numerW1
    val bp    = fracW + manW

    val zProdMoreThan2           = bit(fracW + 1 + manW + 1 - 1, zProd)
    val roundBits                = bp - manW + zProdMoreThan2
    val zProdRounded             = roundBySpec(RoundSpec.roundToEven, roundBits, zProd)
    val zProdMoreThan2AfterRound = bit(manW + 1, zProdRounded)
    val zMan0                    = zProdRounded >> zProdMoreThan2AfterRound

    val zMan = if (zEx == 0 || xySameMan) {
      SafeLong(0)
    } else if (maxxy.man == SafeLong(0)) {
      minxy.man
    } else {
      zMan0
    }

    // ---- encode status bits + ysgn ----
    assert(zEx <= spec.exBias)

    val zSgnEnc = ysgn
    val zExEnc  = zEx | (atan2Status1 << (spec.exW - 1))
    val zManEnc = (zMan & ((SafeLong(1) << spec.manW) - 2)) + atan2Status0

    new RealGeneric(spec, zSgnEnc, zExEnc, zManEnc)
  }
}

// =========================================================================
// Phase 2: t -> atan2(y, x)
// =========================================================================

object ATan2Phase2SimGeneric {

  def calcLinearThreshold(manW: Int): Int = { -math.ceil(manW / 2.0 + 1.0).toInt }

  def apply(tATan2: PolynomialCoefficients, x: RealGeneric): RealGeneric = {
    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val addrW     = tATan2.bitW.addrW
    val outPrec   = tATan2.bitW.outPrec
    val fracW     = outPrec
    val extraBits = outPrec - manW

    // ---- special values (NaN-boxed from phase 1) ----
    if (x.ex == maskI(exW)) {
      val sgnCoef = if (x.sgn == 0) { 1 } else { -1 }
      if      (x.man == SafeLong(1)) { return new RealGeneric(spec, x.sgn, 0, 0) }
      else if (x.man == SafeLong(2)) { return new RealGeneric(spec, sgnCoef * Pi) }
      else if (x.man == SafeLong(3)) { return new RealGeneric(spec, sgnCoef * Pi * 0.5) }
      else if (x.man == SafeLong(4)) { return new RealGeneric(spec, sgnCoef * Pi * 0.25) }
      else if (x.man == SafeLong(5)) { return new RealGeneric(spec, sgnCoef * Pi * 0.75) }
      else                           { return RealGeneric.nan(spec) }
    }

    // ---- decode status ----
    val xsgn = 0
    val xex  =  x.ex  & ((1 << (exW - 1)) - 1) // strip MSB (status1)
    val xman = (x.man & ((SafeLong(1) << manW) - 2)).toLong // strip LSB (status0)
    val xmanW1 = SafeLong(xman) + (SafeLong(1) << manW)

    val status = (((x.ex >> (exW - 1)) & 1) << 1) | (x.man & 1).toInt
    val ysgn   = x.sgn

    // ---- atan(t)/t table ----
    val xShift0 = exBias - xex
    val xShift  = Seq(xShift0, manW + 1).min
    val xFixed  = (xmanW1 | SafeLong(1)) >> xShift
    assert(0 < xShift0, f"xShift0 = ${xShift0}")

    val (resManW1, resEx) = if (xex <= exBias + calcLinearThreshold(manW)) {
      // Very small t: linear region; atan(t)/t ~ 1.
      (SafeLong(1) << fracW, 0)
    } else {
      val xW   = manW - addrW
      val addr = slice(xW, addrW, xFixed).toInt
      val xiv  = slice(0,  xW,    xFixed).toBigInt
      val res0 = PolynomialEval.eval(addr, xiv, tATan2)
      assert(res0 <= (SafeLong(1) << fracW),
        f"xFixed = ${xFixed}, res0 = ${res0} > 2^${fracW}")
      (SafeLong(res0) << 1, -1)
    }

    // ---- compute t * atan(t)/t = atan(t) ----
    val atanProd           = resManW1 * xmanW1
    val atanProdMoreThan2  = bit((fracW + 1) + (manW + 1) - 1, atanProd).toInt
    val atanRoundBit       = fracW + atanProdMoreThan2
    val atanRound          = roundBySpec(RoundSpec.roundToEven, atanRoundBit, atanProd)
    val atanRoundMoreThan2 = bit(manW + 1, atanRound).toInt
    val atanMan            = slice(0, manW, atanRound >> atanRoundMoreThan2)
    val atanEx             = xex + resEx + atanProdMoreThan2 + atanRoundMoreThan2

    // ---- combine with pi/pi-2 depending on status ----
    val zsgn        = ysgn
    val piManW1     = Real.pi(manW + 2)
    val halfpiManW1 = (Real.pi / Real.two)(manW + 2)
    val atanManW1   = (atanMan + (SafeLong(1) << manW)) << 2

    assert(0 <= atanEx && atanEx < exBias)
    val atanShift   = exBias - atanEx
    val atanAligned = atanManW1 >> atanShift

    if (status == 0) {
      //   ysgn *  atan(|y|/|x|)        (x>0, |x|>|y|)
      new RealGeneric(spec, zsgn, atanEx, atanMan)
    } else if (status == 1) {
      //   ysgn * (pi - atan(|y|/|x|))  (x<0, |x|>|y|)
      val zman0 = piManW1 - atanAligned
      val zman  = roundBySpec(RoundSpec.round, 3, zman0)
      new RealGeneric(spec, zsgn, 1 + exBias, zman - (SafeLong(1) << manW))
    } else if (status == 2) {
      //   ysgn * (pi/2 - atan(|x|/|y|)) (x>0, |x|<|y|)
      val zman0          = halfpiManW1 - atanAligned
      val zman0LessThan1 = if (bit(manW + 2, zman0) == 0) 1 else 0
      val zmanRound      = roundBySpec(RoundSpec.round, 2 - zman0LessThan1, zman0)
      val zman           = zmanRound - (SafeLong(1) << manW)
      val zex            = exBias - zman0LessThan1
      new RealGeneric(spec, zsgn, zex, zman)
    } else {
      //   ysgn * (pi/2 + atan(|x|/|y|)) (x<0, |x|<|y|)
      assert(status == 3)
      val zman0         = halfpiManW1 + atanAligned
      val zmanRound     = roundBySpec(RoundSpec.round, 2, zman0)
      val zmanMoreThan2 = bit(manW + 1, zmanRound)
      val zman          = slice(zmanMoreThan2, manW, zmanRound)
      val zex           = zmanMoreThan2 + exBias
      new RealGeneric(spec, zsgn, zex, zman)
    }
  }
}

/** Software entry point: atan2(y, x) = phase2(phase1(y, x)) with the `standard`
 *  FP32 configuration. For other configurations, see [[rial.math.MathSim]]. */
object ATan2Phase2Sim {
  import FuncKind.{Reciprocal, ATan2Phase2}

  private lazy val cfg   = MathFuncConfig.standardFP32
  private lazy val tRec  = cfg.coeffs(Reciprocal)    // phase 1 reuses the reciprocal table
  private lazy val tATan = cfg.coeffs(ATan2Phase2)

  def atan2DefaultFP32(y: Float, x: Float): Float = {
    val yg = RealGeneric.fromFloat(RealSpec.Float32Spec, y)
    val xg = RealGeneric.fromFloat(RealSpec.Float32Spec, x)
    ATan2Phase2SimGeneric(tATan, ATan2Phase1SimGeneric(tRec, yg, xg)).toFloat
  }
}

// =========================================================================
// atan(t) for |t| <= 1, without either shared multiplier
// =========================================================================

// The table approximates atan(t) itself for t in [0, 1).
//   atan(t) is in [0, pi/4), so the value is stored as it is.
class ATanFixedTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  private val func = (t: Double) => { atan(t) }

  private val nIntervals = 1 << addrW
  private val ft  = FuncTable(func, nOrder, nIntervals, 0.0, 1.0)
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW)

  def getBitWidth(): PolynomialBitWidth = { gen.getBitWidth() }
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    gen.generateTable(bitW)
  }
}

/** Software simulator for the hardware that computes atan(t) for |t| <= 1
 *  without either shared multiplier.
 *
 *  Its table holds atan(t) itself, so the polynomial output is a fixed point
 *  value and the result assembly normalises it instead of multiplying it by
 *  the argument.
 */
object ATanFixedSimGeneric {

  /** Compute atan(t).
   *
   *  @param t polynomial coefficients (from ATanFixedTable)
   *  @param x the argument, whose absolute value must be at most one
   */
  def apply(t: PolynomialCoefficients, x: RealGeneric): RealGeneric = {

    val spec   = x.spec
    val manW   = spec.manW
    val exBias = spec.exBias

    if (x.isNaN) { return RealGeneric.nan(spec) }

    // |t| >= 1 gives atan(1) = pi/4, the value the accurate atan2 gives for
    // |y| == |x|; above 1 it is finite but meaningless.
    if (x.ex >= exBias) {
      val quarterPi = new RealGeneric(spec, Pi * 0.25)
      return new RealGeneric(spec, x.sgn, quarterPi.ex, quarterPi.man)
    }
    if (x.isZero) { return new RealGeneric(spec, x.sgn, 0, SafeLong(0)) }

    // atan(t) is t to the working precision below this exponent; the accurate
    // second phase of atan2 shortens the same range the same way.
    if (x.ex <= exBias + ATan2Phase2SimGeneric.calcLinearThreshold(manW)) {
      return new RealGeneric(spec, x.sgn, x.ex, x.man)
    }

    // -----------------------------------------------------------------------
    // table lookup at |t| aligned at the binary point

    val addrW   = t.bitW.addrW
    val outPrec = t.bitW.outPrec

    val xShift = exBias - x.ex
    val xFixed = x.manW1 >> xShift

    val xW = manW - addrW
    val (addr, xiv) = if (t.bitW.nOrder == 0) {
      (xFixed.toInt, BigInt(0))
    } else {
      (slice(xW, addrW, xFixed).toInt, slice(0, xW, xFixed).toBigInt)
    }

    val zres = SafeLong(PolynomialEval.eval(addr, xiv, t))
    assert(SafeLong(0) < zres && zres < (SafeLong(1) << outPrec),
      f"xFixed = ${xFixed}, zres = ${zres} is not in (0, 2^${outPrec})")

    // -----------------------------------------------------------------------
    // the table gave atan(t) in [0, 1) as a fixed point value: its leading
    // zeros give the exponent, and what is left of it the mantissa.

    val count   = outPrec - binaryWidthSL(zres)
    val shifted = zres << count

    val extraBits            = outPrec - 1 - manW
    val zRounded             = roundBySpec(RoundSpec.roundToEven, extraBits, shifted)
    val zMoreThan2AfterRound = bit(manW + 1, zRounded)

    val zMan = slice(0, manW, zRounded)
    val zEx  = exBias - 1 - count + zMoreThan2AfterRound

    new RealGeneric(spec, x.sgn, zEx, zMan)
  }
}
