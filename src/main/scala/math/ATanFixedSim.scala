//% @file ATanFixedSim.scala
//
// Simulator for atan(t), |t| <= 1, without either shared multiplier.
//
// atan(t) vanishes with t: a table of atan(t) indexed by t itself would
// answer a small t with the few bits its sub-interval happens to hold there.
// The table below is indexed by the binade of |t| instead, one region per
// binade, and each region stores atan(t) scaled by the power of two that
// brings it into the same range as every other region, so the result keeps
// its precision however small t is. Once t is small enough for atan(t) to be
// t to the working precision, the module hands the argument back and no
// region is needed for it.
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math._
import java.lang.Math.scalb

import spire.math.SafeLong
import spire.implicits._

import rial.table._

import rial.util.ScalaUtil._

import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding._

/** The table of atan(t) for |t| in (0, 1), one region per binade of |t|.
 *
 *  The binade e, i.e. |t| = (1 + m) * 2^-e with m in [0, 1), stores
 *  atan(t) * 2^(e-1), which is in [1/4, 1) for every e. The widest binades
 *  are halved by the top bit of m: they are the ones atan bends over, and a
 *  halved region leaves room in the address for every bit of the mantissa.
 *  The result assembly undoes the scale in the exponent.
 *
 *  @param manW mantissa width of the argument
 */
class ATanFixedTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
  val manW:    Int,
) extends TableBase {

  /** The regions of the table and the address that selects one. */
  val layout = ATanFixedTable.layout(nOrder, addrW, manW)

  private val gen = new PolynomialGenerator(
    layout.funcTable(nOrder, manW, inBitW)(ATanFixedTable.value),
    inBitW, outBitW)

  /** The width of the sub-interval index is the one of a region; the address
   *  the hardware forms carries the region above it. */
  def getBitWidth(): PolynomialBitWidth = {
    gen.getBitWidth().copy(addrW = addrW)
  }

  /** Concatenated form: the regions in address order, zero-filled to the
   *  power of two the address is. Consumed by the HW (ATanFixedTableCoeff)
   *  and by [[rial.math.ATanFixedSimGeneric]]. */
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    layout.generateTable(gen, bitW)
  }
}

/** The layout of [[rial.math.ATanFixedTable]] and the value it stores. */
object ATanFixedTable {

  /** The widest binade of the argument: |t| is below one. */
  val firstBinade: Int = 1

  /** The narrowest binade the table holds; below it atan(t) is t itself,
   *  which the module hands back unchanged.
   *
   *  @param manW mantissa width of the argument
   */
  def lastBinade(manW: Int): Int = {
    -ATan2Phase2SimGeneric.calcLinearThreshold(manW) - 1
  }

  /** How many of the widest binades are halved by the top bit of the
   *  mantissa; every other binade the table holds is one region.
   *
   *  They are the binades atan bends over, and a halved region leaves room in
   *  the address for every bit of the mantissa.
   */
  val nHalvedBinade: Int = 4

  /** Sub-intervals a region is divided into, as a power of two.
   *
   *  A pure lookup indexes with every mantissa bit it has; a polynomial needs
   *  only a few sub-intervals per region, and spends the rest of the address
   *  on the region.
   */
  def subW(nOrder: Int, addrW: Int): Int = {
    if (nOrder == 0) {
      addrW
    } else {
      addrW - 4
    }
  }

  /** The regions of the table and the address that selects one.
   *
   *  A pure lookup reads the mantissa itself, so it needs no region narrower
   *  than a binade; a polynomial halves [[nHalvedBinade]] of them, which
   *  fills the address with as many regions as it holds.
   *
   *  @param manW mantissa width of the argument
   */
  def layout(nOrder: Int, addrW: Int, manW: Int): GradedLayout = {
    val splitW = (e: Int) => {
      if (nOrder != 0 && e < firstBinade + nHalvedBinade) {
        1
      } else {
        0
      }
    }
    new GradedLayout(
      (firstBinade to lastBinade(manW)).map(e => GradedBinade(e, splitW(e))),
      subW(nOrder, addrW))
  }

  /** The value a region stores for the mantissa `m` of its binade:
   *  atan(t) of t = (1 + m) * 2^-e, scaled into [1/4, 1). */
  def value(r: GradedRegion, m: Double): Double = {
    atan(scalb(1.0 + m, -r.binade)) * scalb(1.0, r.binade - 1)
  }
}

/** Software simulator for the hardware that computes atan(t) for |t| <= 1
 *  without either shared multiplier.
 *
 *  Its table holds atan(t) itself, one region per binade of |t|, so the
 *  polynomial output is the result as a fixed point value and the result
 *  assembly normalises it instead of multiplying it by the argument.
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
    // the region the binade of |t| falls in, and the fraction that addresses
    // the sub-interval inside it

    val nOrder = t.bitW.nOrder
    val inBitW = t.bitW.inBitW
    val fracW  = t.bitW.outPrec
    val layout = ATanFixedTable.layout(nOrder, t.bitW.addrW, manW)

    val binade = exBias - x.ex
    val region = layout.regionIndex(binade, x.man, manW)
    val key    = layout.key(binade, x.man, manW)

    val addr = (region << layout.subW) +
               TableAddress.adr(key, manW, layout.subW)
    val dx   = TableAddress.dx(key, manW, layout.subW, inBitW)

    val res0 = PolynomialEval.eval(addr, dx, t)

    // the hardware clamps whatever leaves the polynomial to [0, 1)
    val zres = if (res0 < 0 || res0 >= (BigInt(1) << fracW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial out of range " +
              f"at x = ${x.toDouble}, region = ${region}, res = ${res0}")
      maskSL(fracW)
    } else {
      SafeLong(res0)
    }

    // -----------------------------------------------------------------------
    // every region holds at least a quarter, so the value has at most two
    // leading zeros, and atan(t) is it shifted by 2^(1 - binade)

    val lead = if (bit(fracW - 1, zres) == 1) {
      0
    } else if (bit(fracW - 2, zres) == 1) {
      1
    } else {
      2
    }

    val shifted   = zres << lead
    val extraBits = fracW - 1 - manW

    val zRounded             = roundBySpec(RoundSpec.roundToEven, extraBits, shifted)
    val zMoreThan2AfterRound = bit(manW + 1, zRounded)

    val zMan = slice(0, manW, zRounded)
    val zEx  = x.ex - lead + zMoreThan2AfterRound

    new RealGeneric(spec, x.sgn, zEx, zMan)
  }
}
