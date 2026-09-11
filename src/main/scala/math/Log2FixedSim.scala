//% @file Log2FixedSim.scala
//
// Simulator for log2(x) without either shared multiplier.
//
// One table, log2(1+t) for t in [0, 1): log2(x) = e + log2(1.man), so the
// polynomial output is the fraction of the result and the exponent of the
// argument its integer part.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._

import spire.math.SafeLong
import spire.implicits._

import rial.table._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding._

// The table approximates  z(t) = log2(1 + t)  for t in [0, 1)
// z is in [0, 1), so the value is stored as it is.
class Log2FixedTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  private val func = (t: Double) => {
    log(1.0 + t) / log(2.0)
  }

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

object Log2FixedTable {
  def fp32(): Log2FixedTable = {
    val addrW     = 8
    val extraBits = 3
    new Log2FixedTable(
      nOrder  = 2,
      addrW   = addrW,
      inBitW  = (RealSpec.Float32Spec.manW - addrW),
      outBitW = (RealSpec.Float32Spec.manW + extraBits),
    )
  }
  def bf16(): Log2FixedTable = {
    val addrW     = RealSpec.BFloat16Spec.manW
    val extraBits = 1
    new Log2FixedTable(
      nOrder  = 0,
      addrW   = addrW,
      inBitW  = (RealSpec.BFloat16Spec.manW - addrW),
      outBitW = (RealSpec.BFloat16Spec.manW + extraBits),
    )
  }
}

/** Software simulator for the hardware that computes log2(x) without either
 *  shared multiplier.
 *
 *  Its table holds log2(1 + t) itself, so the polynomial output is the
 *  fraction of the result and the result assembly normalises it together with
 *  the exponent of the argument instead of multiplying it by a constant.
 */
object Log2FixedSimGeneric {

  /** Compute log2(x).
   *
   *  @param t polynomial coefficients (from Log2FixedTable)
   *  @param x input
   */
  def apply(t: PolynomialCoefficients, x: RealGeneric): RealGeneric = {

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    //   log2(nan) -> nan      log2(inf) -> inf
    //   log2(0)   -> -inf     log2(neg) -> nan      log2(1) -> +0
    if (x.isNaN)      { return RealGeneric.nan(spec)           }
    if (x.isZero)     { return RealGeneric.inf(spec, sgn = 1)  }
    if (x.sgn == 1)   { return RealGeneric.nan(spec)           }
    if (x.isInfinite) { return RealGeneric.inf(spec, sgn = 0)  }
    if (x.ex == exBias && x.man == SafeLong(0)) {
      return RealGeneric.zero(spec)
    }

    // -----------------------------------------------------------------------
    // table lookup at the mantissa

    val addrW = t.bitW.addrW
    val fracW = t.bitW.outPrec

    val xW   = manW - addrW
    val addr = slice(xW, addrW, x.man).toInt
    val xiv  = slice(0,  xW,    x.man).toBigInt

    val res0 = PolynomialEval.eval(addr, xiv, t)
    val res = if (res0 < 0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
      SafeLong(0)
    } else if (res0 >= (SafeLong(1) << fracW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
      maskSL(fracW)
    } else {
      SafeLong(res0)
    }

    // -----------------------------------------------------------------------
    // the magnitude of the result as a fixed point value
    //
    //  x >= 1 (e >= 0):  log2(x) =    e     + frac
    //  x <  1 (e <  0):  log2(x) = -((|e|-1) + (1 - frac))
    //
    // A frac of zero leaves 1 - frac at one, an increment of the integer part.

    val exNobias = x.ex - exBias
    val zsgn     = if (exNobias < 0) { 1 } else { 0 }
    val zint0    = if (exNobias >= 0) {
      SafeLong(exNobias)
    } else {
      SafeLong(-exNobias - 1)
    }

    val (zint, zfrac) = if (zsgn == 0) {
      (zint0, res)
    } else if (res == SafeLong(0)) {
      (zint0 + 1, SafeLong(0))
    } else {
      (zint0, (SafeLong(1) << fracW) - res)
    }

    val mag = (zint << fracW) + zfrac
    assert(mag > 0, f"log2 magnitude is zero at x = ${x.toDouble}")

    // -----------------------------------------------------------------------
    // the leading zeros of the magnitude give the exponent, and what is left
    // of it the mantissa

    val magW    = exW + fracW
    val count   = magW - binaryWidthSL(mag)
    val shifted = mag << count

    val extraBits            = magW - 1 - manW
    val zRounded             = roundBySpec(RoundSpec.roundToEven, extraBits, shifted)
    val zMoreThan2AfterRound = bit(manW + 1, zRounded)

    val zMan = slice(0, manW, zRounded)
    val zEx  = exBias + exW - 1 - count + zMoreThan2AfterRound

    new RealGeneric(spec, zsgn, zEx, zMan)
  }
}

/** Software entry point: log2(x) with the default FP32 table.
 *  For other configurations, see [[rial.math.MathSim]]. */
object Log2FixedSim {

  private lazy val fp32Table =
    Log2FixedTable.fp32().generateTable(Log2FixedTable.fp32().getBitWidth())

  /** log2 simulation with the default fp32 configuration */
  def log2DefaultFP32(x: Float): Float = {
    Log2FixedSimGeneric(fp32Table,
      RealGeneric.fromFloat(RealSpec.Float32Spec, x)).toFloat
  }
}
