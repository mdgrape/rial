//% @file SinPiFixedSim.scala
//
// Simulator for sin(pi x) and cos(pi x) without either shared multiplier.
//
// The reduction leaves an argument y in [0, 1/2), and sin(pi y) vanishes with
// y: a table of sin(pi y) indexed by y itself would answer a small y with the
// few bits its sub-interval happens to hold there. The table below is indexed
// by the binade of y instead, one region per binade, and each region stores
// sin(pi y) scaled by the power of two that brings it into the same range as
// every other region, so the result keeps its precision however small y is.
// Once y is small enough for sin(pi y) to be pi y to the working precision,
// every binade stores the same straight line, and one region answers them all.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._
import java.lang.Math.scalb

import spire.math.SafeLong
import spire.implicits._

import rial.table._

import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding._

/** The table of sin(pi y) for the reduced argument y in (0, 1/2), one region
 *  per binade of y.
 *
 *  The binade e, i.e. y = (1 + m) * 2^-e with m in [0, 1), stores
 *  sin(pi y) * 2^(e-3), which is in [0.35, 0.79) for every e. The two widest
 *  binades bend too much for one region and are halved by the top bit of m;
 *  once sin(pi y) has settled into the straight line pi y, every narrower
 *  binade holds the same values, so one region answers them all. The result
 *  assembly undoes the scale in the exponent.
 *
 *  @param manW mantissa width of the argument
 */
class SinPiFixedTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
  val manW:    Int,
) extends TableBase {

  /** The regions of the table and the address that selects one. */
  val layout = SinPiFixedTable.layout(nOrder, addrW)

  // the reduction leaves a mantissa of `manW` bits, which is what a region
  // is centred on; an argument the reduction cancelled down to fewer of them
  // reads the same region
  private val gen = new PolynomialGenerator(
    layout.funcTable(nOrder, _ => manW, inBitW)(SinPiFixedTable.value),
    inBitW, outBitW)

  /** The width of the sub-interval index is the one of a region; the address
   *  the hardware forms carries the region above it. */
  def getBitWidth(): PolynomialBitWidth = {
    gen.getBitWidth().copy(addrW = addrW)
  }

  /** Concatenated form: the regions in address order, zero-filled to the
   *  power of two the address is. Consumed by the HW (SinPiFixedTableCoeff)
   *  and by [[rial.math.SinPiFixedSimGeneric]]. */
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    layout.generateTable(gen, bitW)
  }
}

/** The layout of [[rial.math.SinPiFixedTable]] and the value it stores. */
object SinPiFixedTable {

  /** The widest binade of the reduced argument: y is below a half. */
  val firstBinade: Int = 2

  /** The narrowest binade the table holds; every narrower one reads its
   *  region.
   *
   *  By this binade sin(pi y) has settled into the straight line pi y to well
   *  within the resolution of the table, and a narrower binade stores that
   *  same line scaled by the same power of two, so one region answers them
   *  all.
   */
  val lastBinade: Int = 15

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
   *  than a binade; a polynomial halves the two widest binades, over which
   *  one quadratic per sub-interval would not hold sin(pi y).
   */
  def layout(nOrder: Int, addrW: Int): GradedLayout = {
    val splitW = (e: Int) => {
      if (nOrder != 0 && (e == firstBinade || e == firstBinade + 1)) {
        1
      } else {
        0
      }
    }
    new GradedLayout(
      (firstBinade to lastBinade).map(e => GradedBinade(e, splitW(e))),
      subW(nOrder, addrW))
  }

  /** The value a region stores for the mantissa `m` of its binade:
   *  sin(pi y) of y = (1 + m) * 2^-e, scaled into [1/4, 1). */
  def value(r: GradedRegion, m: Double): Double = {
    sin(Pi * scalb(1.0 + m, -r.binade)) * scalb(1.0, r.binade - 3)
  }
}

/** Software simulator for the hardware that computes sin(pi x) and cos(pi x)
 *  without either shared multiplier.
 *
 *  Its table holds sin(pi y) itself, one region per binade of the reduced
 *  argument, so the polynomial output is the result as a fixed point value
 *  and the result assembly normalises it instead of multiplying it by that
 *  argument.
 */
object SinPiFixedSimGeneric {

  /** Compute sin(pi x) or cos(pi x).
   *
   *  @param isSin true for sin(pi x), false for cos(pi x)
   *  @param t     polynomial coefficients (from SinPiFixedTable)
   *  @param x     input
   */
  def apply(isSin: Boolean, t: PolynomialCoefficients, x: RealGeneric):
      RealGeneric = {

    val spec   = x.spec
    val manW   = spec.manW
    val exBias = spec.exBias

    if (x.isNaN)      { return RealGeneric.nan(spec) }
    if (x.isInfinite) { return RealGeneric.nan(spec) }

    val y = SinCosReduction.reduce(isSin, x, None, scaled = true)

    // y == 1/2 exactly  ->  result is +-1
    if (y.yex == exBias - 1) {
      return new RealGeneric(spec, y.zsgn, exBias, SafeLong(0))
    }
    // y == 0  ->  result is +-0
    if (y.yex == 0) {
      return new RealGeneric(spec, y.zsgn, 0, SafeLong(0))
    }

    // -----------------------------------------------------------------------
    // the region the binade of y falls in, and the fraction that addresses
    // the sub-interval inside it

    val nOrder = t.bitW.nOrder
    val inBitW = t.bitW.inBitW
    val fracW  = t.bitW.outPrec
    val layout = SinPiFixedTable.layout(nOrder, t.bitW.addrW)

    val binade = exBias - y.yex
    val region = layout.regionIndex(binade, y.yman, manW)
    val key    = layout.key(binade, y.yman, manW)

    val addr = (region << layout.subW) +
               TableAddress.adr(key, manW, layout.subW)
    val dx   = if (nOrder == 0) {
      BigInt(0)
    } else {
      TableAddress.dx(key, manW, layout.subW, inBitW)
    }

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
    // leading zeros, and sin(pi y) is it shifted by 2^(3 - binade)

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
    val zEx  = y.yex + 2 - lead + zMoreThan2AfterRound

    new RealGeneric(spec, y.zsgn, zEx, zMan)
  }
}
