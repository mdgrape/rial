//% @file ACosFixedSim.scala
//
// Simulator for acos(x) without either shared multiplier.
//
// acos has a square-root singularity at x = +-1, so no polynomial in x alone
// holds it near the ends. What holds it is the binade of u = 1 - |x|:
// acos(1 - u) is sqrt(2u) times a smooth function of u, so over a binade a
// few sub-intervals of a quadratic in the mantissa of u are enough, and the
// same is true of acos(x) in x itself while |x| stays below a half. The table
// is therefore one region per binade plus one for the argument itself, and
// the value a binade stores is scaled by a power of two so that every region
// fills the same range and the polynomial output keeps its precision as the
// angle approaches zero.
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math._
import java.lang.Math.scalb

import spire.math.Real
import spire.math.SafeLong
import spire.implicits._

import chisel3.util.log2Ceil

import rial.table._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding._

/** The table of acos(x), one region after another: the argument itself over
 *  a in [0, 1/2), then the binades of u = 1 - |x|, each over the mantissa m
 *  of u in [0, 1).
 *
 *  Region 0 holds acos(a) halved, in (0.52, 0.79). The binade e holds
 *  acos(1 - u) of u = (1 + m) * 2^-e, scaled by 2^floor(e/2) / 4, which puts
 *  every binade in [0.25, 0.55). The result assembly undoes the scale in the
 *  exponent.
 */
class ACosFixedTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  /** sub-intervals a region is divided into, as a power of two. */
  val subW = ACosFixedTable.subW(nOrder, addrW)

  /** regions: the argument itself, and one per binade of the difference. */
  val nRegion = ACosFixedTable.nRegion(outBitW)

  private val nSub = 1 << subW

  private val funcArgument = (a: Double) => { acos(a) * 0.5 }

  // acos(1 - u) is 2 asin(sqrt(u/2)), which no cancellation spoils as u -> 0
  private def funcBinade(e: Int) = (m: Double) => {
    val u = scalb(1.0 + m, -e)
    2.0 * asin(sqrt(u * 0.5)) * scalb(0.25, e / 2)
  }

  private val tables =
    Seq(FuncTable(funcArgument, nOrder, nSub, 0.0, 0.5)) ++
    (1 until nRegion).map(e => FuncTable(funcBinade(e), nOrder, nSub, 0.0, 1.0))

  private val gens = tables.map(ft => new PolynomialGenerator(ft, inBitW, outBitW))

  /** The width of the sub-interval index is the one of a region; the address
   *  the hardware forms carries the region above it. */
  def getBitWidth(): PolynomialBitWidth = {
    gens.map(_.getBitWidth()).reduce((l, r) => l.merge(r)).copy(addrW = addrW)
  }

  /** Concatenated form: the regions in address order, zero-filled to the
   *  power of two the address is. Consumed by the coefficient ROM of the
   *  HW and by [[rial.math.ACosFixedSimGeneric]]. */
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    val cs    = gens.map(_.generateTable(bitW).coeffs).reduce(_ ++ _)
    val zeros = Seq.fill(nOrder + 1)(0L)
    val nRow  = 1 << (ACosFixedTable.regionW(outBitW) + subW)
    require(cs.length <= nRow,
      f"${cs.length}%d table entries do not fit the ${nRow}%d the address has")
    PolynomialCoefficients(bitW, cs ++ Seq.fill(nRow - cs.length)(zeros))
  }
}

/** The layout of [[rial.math.ACosFixedTable]] and the constants the result
 *  assembly needs at its scale. */
object ACosFixedTable {

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
      addrW - 3
    }
  }

  /** Regions of the table: the argument itself, and one per binade of
   *  1 - |x|.
   *
   *  The difference has at most manW + 1 fraction bits, so its binades reach
   *  2^-(manW + 1); `outBitW`, which holds a hidden one above a mantissa, is
   *  at least that many and covers them all.
   */
  def nRegion(outBitW: Int): Int = 1 + outBitW

  /** Width of the region index the table address carries above the
   *  sub-interval. */
  def regionW(outBitW: Int): Int = log2Ceil(nRegion(outBitW))

  /** Guard bits below the polynomial output at which the angle of a negative
   *  argument is subtracted from pi. */
  val guardW: Int = 2

  /** pi in the units the result assembly counts the angle in, 2^-fracW. */
  def pi(fracW: Int): SafeLong = {
    Real.pi(fracW)
  }
}

/** Software simulator for the hardware that computes acos(x) without either
 *  shared multiplier.
 *
 *  Its table holds the angle itself, one region per binade of 1 - |x|, so the
 *  polynomial output is the angle as a fixed point value and the result
 *  assembly normalises it instead of multiplying it by the argument.
 */
object ACosFixedSimGeneric {

  /** Compute acos(x).
   *
   *  @param t polynomial coefficients (from ACosFixedTable)
   *  @param x the argument, whose absolute value must be at most one
   */
  def apply(t: PolynomialCoefficients, x: RealGeneric): RealGeneric = {

    val spec   = x.spec
    val manW   = spec.manW
    val exBias = spec.exBias

    // acos is defined on [-1, 1] alone; a finite argument beyond an end gets
    // the value of that end, an infinity or a NaN gives NaN
    if (x.isNaN || x.isInfinite) {
      return RealGeneric.nan(spec)
    }
    if (x.ex >= exBias) {
      if (x.sgn == 0) {
        return RealGeneric.zero(spec)
      } else {
        return new RealGeneric(spec, Pi)
      }
    }
    if (x.isZero) { return new RealGeneric(spec, Pi * 0.5) }

    // -----------------------------------------------------------------------
    // the region the argument falls in, and the fraction that addresses it

    val nOrder = t.bitW.nOrder
    val inBitW = t.bitW.inBitW
    val fracW  = t.bitW.outPrec
    val subW   = ACosFixedTable.subW(nOrder, t.bitW.addrW)

    val (region, key) = if (x.ex < exBias - 1) {
      // |x| < 1/2, where acos is smooth in the argument itself
      val shift = min(exBias - 1 - x.ex, manW + 1)
      (0, slice(0, manW, x.manW1 >> shift))
    } else {
      // |x| in [1/2, 1), where u = 1 - |x| is exact and its binade answers
      val u    = (SafeLong(1) << (manW + 1)) - x.manW1
      val lead = manW + 1 - binaryWidthSL(u)
      (lead + 1, slice(0, manW, u << lead))
    }

    val addr = (region << subW) + TableAddress.adr(key, manW, subW)
    val dx   = if (nOrder == 0) {
      BigInt(0)
    } else {
      TableAddress.dx(key, manW, subW, inBitW)
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
    // the angle: the value of the region undone by the power of two it was
    // stored with, and pi minus it when x is negative

    val guardW = ACosFixedTable.guardW
    val angleW = 2 + fracW + guardW

    val lead = if (bit(fracW - 1, zres) == 1) {
      0
    } else if (bit(fracW - 2, zres) == 1) {
      1
    } else {
      2
    }

    // the region stores its angle scaled up by this power of two
    val thetaShift = if (region == 0) { 1 } else { region / 2 }

    val (shifted, zex0) = if (x.sgn == 0) {
      (slice(0, angleW, zres << (lead + 2 + guardW)),
       exBias + 2 - thetaShift - 1 - lead)
    } else {
      val theta    = (zres << (2 + guardW)) >> thetaShift
      val angle    = ACosFixedTable.pi(fracW + guardW) - theta
      val atLeast2 = bit(angleW - 1, angle)
      (slice(0, angleW, angle << (1 - atLeast2)), exBias + atLeast2)
    }

    // -----------------------------------------------------------------------
    // the leading one of the angle sits at the top of the word: what is left
    // of it is the mantissa

    val extraBits            = angleW - 1 - manW
    val zRounded             = roundBySpec(RoundSpec.roundToEven, extraBits, shifted)
    val zMoreThan2AfterRound = bit(manW + 1, zRounded)

    val zMan = slice(0, manW, zRounded)
    val zEx  = zex0 + zMoreThan2AfterRound

    new RealGeneric(spec, 0, zEx, zMan)
  }
}
