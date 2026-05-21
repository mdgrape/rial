//% @file ExpSim.scala
//
// Simulator for exp(x).
//
// Algorithm:
//   1. multiply x by log2(e)              (emulated PreProcMultiplier)
//   2. split the product into integer + fractional parts; the integer part
//      becomes the result exponent, the fractional part feeds the table
//   3. table approximates 2^t - 1 for t in [0, 1)
//   4. a small correction term (ln(2) * shifted-out bits) is added back to
//      compensate for the rounding done at step 2.
//
// The log2(e) precision is exposed as `log2eExtraW` (default = 14, matches the
// empirical setting from the original math/ port) so the eventual hardware
// can use the same value and stay bit-identical.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._

import spire.math.SafeLong
import spire.math.Real
import spire.implicits._

import rial.approx._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding
import rial.arith.Rounding._

// The table approximates  f(t) = 2^t - 1  for t in [0, 1).
class ExpTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  private val func = (t: Double) => pow(2.0, t) - 1.0

  private val nIntervals = 1 << addrW
  private val ft  = FuncTable(func, nOrder, nIntervals, 0.0, 1.0)
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW)

  def getBitWidth(): PolynomialBitWidth = {
    gen.getBitWidth()
  }
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    gen.generateTable(bitW)
  }
}

object ExpTable {
  def fp32(): ExpTable = {
    val addrW     = 8
    val extraBits = 3
    new ExpTable(
      nOrder  = 2,
      addrW   = addrW,
      inBitW  = (RealSpec.Float32Spec.manW - addrW),
      outBitW = (RealSpec.Float32Spec.manW + extraBits),
    )
  }
  def bf16(): ExpTable = {
    val addrW     = RealSpec.BFloat16Spec.manW
    val extraBits = 1
    new ExpTable(
      nOrder  = 0,
      addrW   = addrW,
      inBitW  = (RealSpec.BFloat16Spec.manW - addrW),
      outBitW = (RealSpec.BFloat16Spec.manW + extraBits),
    )
  }

  private lazy val fp32Table = fp32().generateTable(fp32().getBitWidth())
  def expDefaultFP32(x: Float): Float = {
    ExpSimGeneric(fp32Table,
      RealGeneric.fromFloat(RealSpec.Float32Spec, x)).toFloat
  }
}

object ExpSimGeneric {

  /** Default precision-pad bits for the log2(e) constant. The HW must use
   *  the same value (matches the original math/ port's empirical choice).
   */
  def defaultLog2eExtraW(spec: RealSpec): Int = { 14 }

  /** Compute exp(x).
   *
   *  @param t            polynomial coefficients (from ExpTable)
   *  @param x            input
   *  @param log2eExtraW  precision pad for the log2(e) constant. If `None`,
   *                      uses `defaultLog2eExtraW(x.spec)`.
   */
  def apply(
    t:            PolynomialCoefficients,
    x:            RealGeneric,
    log2eExtraW:  Option[Int] = None,
  ): RealGeneric = {

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val xsgn = x.sgn

    if (x.isNaN)                   { return RealGeneric.nan (spec)    }
    if (x.isInfinite && xsgn == 0) { return RealGeneric.inf (spec, 0) }
    if (x.isInfinite && xsgn == 1) { return new RealGeneric(spec, 1.0) }

    val log2eW = log2eExtraW.getOrElse(defaultLog2eExtraW(spec))

    val addrW   = t.bitW.addrW
    val outPrec = t.bitW.outPrec
    val inBitW  = t.bitW.inBitW
    val extraBits = outPrec - manW

    // -----------------------------------------------------------------------
    // determine widths

    val log2 = (a: Double) => { log(a) / log(2.0) }
    val xExOvfLimit = scala.math.ceil(log2(maskL(exW) - exBias)).toInt // e.g. 7 for FP32
    val xExUdfLimit = scala.math.ceil(log2(abs(0 - exBias))).toInt     // e.g. 7 for FP32

    val xIntW    = max(xExOvfLimit, xExUdfLimit).toInt
    val padding  = extraBits
    val extraMan = padding + xIntW

    // -----------------------------------------------------------------------
    // emulated pre-mult: x * log2(e)

    val log2eManW  = manW + log2eW
    val log2eFixed = (Real.one / Real.log2)(log2eManW)

    val xprod          = x.manW1 * log2eFixed
    val xprodW         = (1 + manW) + (1 + log2eManW)
    val xprodMoreThan2 = bit(xprodW - 1, xprod)

    assert(log2eManW + xprodMoreThan2 - extraMan > 0)
    val xprodRoundedW1 = Rounding.roundToEven(
      log2eManW + xprodMoreThan2 - extraMan, xprod)

    assert(bit(manW + 1 + extraMan, xprodRoundedW1) == 1 ||
           bit(manW +     extraMan, xprodRoundedW1) == 1)
    val xprodMoreThan2AfterRounded = bit(manW + 1 + extraMan, xprodRoundedW1)

    val xman = slice(0, manW + extraMan, xprodRoundedW1)
    val xex  = x.ex + xprodMoreThan2AfterRounded + xprodMoreThan2

    // -----------------------------------------------------------------------
    // pow2

    val xexNobias = xex - exBias

    if (xsgn == 0 && xExOvfLimit <= xexNobias) { return RealGeneric.inf (spec, 0) }
    if (xsgn == 1 && xExUdfLimit <= xexNobias) { return RealGeneric.zero(spec)    }
    // (the result may still over/underflow after the table+correction, but
    //  the early exit above catches the bulk of out-of-range inputs.)

    val xFracW = manW + padding
    val xValW  = xIntW + xFracW + 1
    val xVal   = (SafeLong(1) << (manW + extraMan)) + xman

    assert(xVal < maskSL(xValW))
    assert(xexNobias < xIntW)

    val xshift      = xIntW - xexNobias
    val xValShifted = (xVal >> xshift) + bit(xshift - 1, xVal) // round-to-nearest
    val xint0       = slice(xFracW, xIntW,  xValShifted)
    val xfrac0      = slice(0,      xFracW, xValShifted)

    // for negative x, 2^(-x) = 2^(-xInt - xFrac) = 2^(-xInt - 1 + (1 - xFrac))
    val (xint, xfrac) = if (xsgn == 0) {
      (xint0, xfrac0)
    } else if (xfrac0 != 0) {
      (-xint0 - 1L, (SafeLong(1) << xFracW) - xfrac0)
    } else {
      (-xint0, SafeLong(0))
    }

    // -----------------------------------------------------------------------
    // table lookup of 2^t - 1

    val dxbp = inBitW - 1
    val addr = slice(padding + dxbp + 1, addrW,    xfrac).toInt
    val xiv  = slice(padding,            inBitW,   xfrac).toBigInt

    val zmanRaw = if (spec == RealSpec.Float64Spec && xfrac == SafeLong(0)) {
      // FP64 + x << 1 path: polynomial result can dip slightly negative; pin to 0.
      SafeLong(0)
    } else {
      PolynomialEval.eval(addr, xiv, t).toSafeLong
    }
    assert(0 <= zmanRaw && zmanRaw < (SafeLong(1) << outPrec),
      f"x = ${x.toDouble}, zmanRaw = ${zmanRaw.toLong.toBinaryString} not in [0, 2^${outPrec})")

    val zmanRound = if (extraBits > 0) {
      (zmanRaw >> extraBits) + bit(extraBits - 1, zmanRaw)
    } else {
      zmanRaw
    }

    // -----------------------------------------------------------------------
    // small correction: add ln(2) * (bits of xfrac shifted out during rounding)

    val corrTermW = min(6, 1 + manW + padding)
    val coeffW    = min(8, manW + 1)

    val coefficient        = (Real.log(Real.two))(coeffW)
    val xShiftedOut        = slice(0, padding, xfrac)
    val zCorrectionCoef0   = coefficient * xShiftedOut

    val corrShift          = coeffW + padding - corrTermW
    assert(corrShift >= 0)
    val zCorrectionCoef    = zCorrectionCoef0 >> corrShift

    val zCorrectionTerm0   = ((SafeLong(1) << manW) + zmanRound) * zCorrectionCoef
    val zCorrectionTerm    = zCorrectionTerm0 >> (manW + 1)
    val zCorrectionShifted = bit(corrTermW - 1, zCorrectionTerm) +
                             bit(corrTermW - 2, zCorrectionTerm)
    assert(0 <= zCorrectionShifted && zCorrectionShifted <= 1)

    val zmanCorrected                = zmanRound + zCorrectionShifted
    val zmanMoreThan2AfterCorrection = bit(manW, zmanCorrected)
    val z                            = slice(0, manW, zmanCorrected)

    val zex = xint + exBias + zmanMoreThan2AfterCorrection

    if (zex >= maskI(exW)) { return RealGeneric.inf (spec, 0) }
    if (zex <= 0)          { return RealGeneric.zero(spec)    }

    new RealGeneric(spec, 0, zex.toInt, z)
  }
}
