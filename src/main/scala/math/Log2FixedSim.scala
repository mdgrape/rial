//% @file Log2FixedSim.scala
//
// Simulator for log2(x) without either shared multiplier.
//
// log2 vanishes at x = 1, and a table indexed by the mantissa of x alone
// loses the arguments that approach it: what the result keeps there is the
// few leading bits of a fraction that carries a fixed number of them. The
// table below is indexed by the binade of the distance of x from one
// instead -- x - 1 on one side of it, 1 - x on the other, both exact in
// fixed point -- and each binade stores |log2(x)| scaled by the power of two
// that brings it into the same range as every other binade, so the result
// keeps its precision however close x comes to one.
//
// Away from one the result is the unbiased exponent E of x plus log2(1 + t)
// of its mantissa t, and |E + log2(1 + t)| is at least one. The regions above
// x = 1 hold log2(1 + d) of the mantissa d, which is that same t, so the same
// table answers there too: its value shifted back by the scale of the binade
// it was read at is the fraction of the result, and the exponent carries the
// integer part.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._
import java.lang.Math.scalb

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

/** The table of |log2(x)| near x = 1, one region per binade of the distance
 *  from it, on each of the two sides.
 *
 *  Above one the distance is d = x - 1, which is the mantissa of x itself;
 *  below one it is d = 1 - x, which is exact in fixed point. Either way the
 *  binade e, i.e. d = (1 + m) * 2^-e with m in [0, 1), stores
 *  |log2(x)| * 2^(e-2), which is in [0.29, 1) for every e on either side.
 *  The binades over which log2 bends too much for one region are split by the
 *  top bits of m; the rest hold the same straight line scaled by the same
 *  power of two. The result assembly undoes the scale in the exponent.
 *
 *  The address carries the side of x = 1 above the region index, so the two
 *  sides are one table, each of them zero-filled to the rows that index
 *  reaches.
 *
 *  @param manW mantissa width of the argument
 */
class Log2FixedTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
  val manW:    Int,
) extends TableBase {

  /** The regions above x = 1, indexed by the binade of x - 1. */
  val above = Log2FixedTable.above(nOrder, addrW, manW)

  /** The regions below x = 1, indexed by the binade of 1 - x. */
  val below = Log2FixedTable.below(nOrder, addrW, manW)

  require(above.regionW == below.regionW,
    f"the two sides of x = 1 share the region index of the address, so they " +
    f"must be the same width: ${above.regionW}%d above and " +
    f"${below.regionW}%d below")

  private val gens = Seq(
    new PolynomialGenerator(
      above.funcTable(nOrder, manW, inBitW)(Log2FixedTable.aboveValue),
      inBitW, outBitW),
    new PolynomialGenerator(
      below.funcTable(nOrder, manW, inBitW)(Log2FixedTable.belowValue),
      inBitW, outBitW))

  /** The width of the sub-interval index is the one of a region; the address
   *  the hardware forms carries the region and the side above it. */
  def getBitWidth(): PolynomialBitWidth = {
    gens.map(_.getBitWidth()).reduce((l, r) => { l.merge(r) }).copy(addrW = addrW)
  }

  /** Concatenated form: the regions above x = 1 in address order, zero-filled
   *  to the rows the region index reaches, then the ones below it. Consumed
   *  by the HW (Log2FixedTableCoeff) and by
   *  [[rial.math.Log2FixedSimGeneric]]. */
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    val cs = Seq(above, below).zip(gens).map({ case (layout, gen) => {
      layout.generateTable(gen, bitW).coeffs
    }}).reduce(_ ++ _)
    PolynomialCoefficients(bitW, cs)
  }
}

/** The layout of [[rial.math.Log2FixedTable]] and the value it stores. */
object Log2FixedTable {

  /** Sub-intervals a region is divided into, as a power of two.
   *
   *  A pure lookup indexes with every mantissa bit it has; a polynomial needs
   *  only a few sub-intervals per region, and spends the rest of the address
   *  on the region and the side.
   */
  def subW(nOrder: Int, addrW: Int): Int = {
    if (nOrder == 0) {
      addrW
    } else {
      addrW - 5
    }
  }

  /** Width of the binade a region is looked up on: the distance from one
   *  reaches 2^-(manW + 1) below it. */
  def binadeW(manW: Int): Int = {
    log2Ceil(manW + 2)
  }

  /** The regions above x = 1, one per binade of d = x - 1.
   *
   *  d is the mantissa of x, so its binades run from 1 to `manW`. A pure
   *  lookup reads the mantissa of d itself and needs no region narrower than
   *  a binade; a polynomial splits the widest binades, over which one
   *  quadratic per sub-interval would not hold log2(1 + d).
   */
  def above(nOrder: Int, addrW: Int, manW: Int): GradedLayout = {
    val splitW = (e: Int) => {
      if (nOrder == 0) {
        0
      } else if (e <= 2) {
        2
      } else if (e <= 4) {
        1
      } else {
        0
      }
    }
    new GradedLayout(
      (1 to manW).map(e => GradedBinade(e, splitW(e))),
      subW(nOrder, addrW))
  }

  /** The regions below x = 1, one per binade of d = 1 - x.
   *
   *  x in (1/2, 1) carries a mantissa of `manW` bits one exponent below one,
   *  so the difference has `manW + 1` fraction bits and its binades run from
   *  2 to `manW + 1`. -log2(1 - d) reaches one as d reaches a half, so the
   *  widest binade here is split further than the one above x = 1.
   */
  def below(nOrder: Int, addrW: Int, manW: Int): GradedLayout = {
    val splitW = (e: Int) => {
      if (nOrder == 0) {
        0
      } else if (e <= 2) {
        3
      } else if (e <= 4) {
        1
      } else {
        0
      }
    }
    new GradedLayout(
      (2 to manW + 1).map(e => GradedBinade(e, splitW(e))),
      subW(nOrder, addrW))
  }

  /** The value a region above x = 1 stores for the mantissa `m` of its
   *  binade: log2(1 + d) of d = (1 + m) * 2^-e, scaled into [1/4, 1). */
  def aboveValue(r: GradedRegion, m: Double): Double = {
    val d = scalb(1.0 + m, -r.binade)
    log(1.0 + d) / log(2.0) * scalb(1.0, r.binade - 2)
  }

  /** The value a region below x = 1 stores for the mantissa `m` of its
   *  binade: -log2(1 - d) of d = (1 + m) * 2^-e, scaled into [1/4, 1). */
  def belowValue(r: GradedRegion, m: Double): Double = {
    val d = scalb(1.0 + m, -r.binade)
    -log(1.0 - d) / log(2.0) * scalb(1.0, r.binade - 2)
  }

  /** The default FP32 table. */
  def fp32(): Log2FixedTable = {
    new Log2FixedTable(
      nOrder  = 2,
      addrW   = 8,
      inBitW  = 18,
      outBitW = RealSpec.Float32Spec.manW + 3,
      manW    = RealSpec.Float32Spec.manW,
    )
  }

  /** The default BFloat16 table. */
  def bf16(): Log2FixedTable = {
    new Log2FixedTable(
      nOrder  = 0,
      addrW   = RealSpec.BFloat16Spec.manW,
      inBitW  = 0,
      outBitW = RealSpec.BFloat16Spec.manW + 1,
      manW    = RealSpec.BFloat16Spec.manW,
    )
  }
}

/** Software simulator for the hardware that computes log2(x) without either
 *  shared multiplier.
 *
 *  Its table holds |log2(x)| itself, one region per binade of the distance of
 *  x from one, so near one the polynomial output is the result as a fixed
 *  point value and the result assembly normalises it. Away from one the same
 *  regions give log2(1 + t) of the mantissa t, and the exponent of the
 *  argument carries the integer part of the result.
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

    val nOrder = t.bitW.nOrder
    val inBitW = t.bitW.inBitW
    val fracW  = t.bitW.outPrec
    val above  = Log2FixedTable.above(nOrder, t.bitW.addrW, manW)
    val below  = Log2FixedTable.below(nOrder, t.bitW.addrW, manW)
    val subW   = above.subW

    // -----------------------------------------------------------------------
    // the distance of x from one, and the region of the binade it lies in
    //
    // x in (1/2, 2) lies as close to one as the table holds |log2 x| of.
    // Below one the distance 1 - x is exact in fixed point (Sterbenz), and is
    // the complement of the mantissa there; away from one the mantissa itself
    // is the t of log2(1 + t).

    val nearOne  = (x.ex == exBias || x.ex == exBias - 1) && x.man != SafeLong(0)
    val belowOne = nearOne && x.ex == exBias - 1
    val whole    = x.man == SafeLong(0)

    val diff = if (belowOne) {
      (SafeLong(1) << manW) - x.man
    } else {
      x.man
    }

    // a distance in [2^-e, 2^-(e-1)) has e - 1 leading zeros in a fraction of
    // manW bits, so the binade it lies in is that count plus one, and one
    // more below x = 1, where the distance is that fraction halved
    val lead = if (diff == SafeLong(0)) {
      manW - 1
    } else {
      manW - binaryWidthSL(diff)
    }
    val binade = lead + 1 + (if (belowOne) { 1 } else { 0 })
    val man    = slice(0, manW, diff << (lead + 1))

    val layout = if (belowOne) { below } else { above }
    val side   = if (belowOne) { 1 } else { 0 }

    val region = layout.regionIndex(binade, man, manW)
    val key    = layout.key(binade, man, manW)

    val addr = (side << (layout.regionW + subW)) + (region << subW) +
               TableAddress.adr(key, manW, subW)
    val dx   = TableAddress.dx(key, manW, subW, inBitW)

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
    // the result as a normalized magnitude of exW + fracW bits, and the
    // exponent its leading one weighs

    val magW      = exW + fracW
    val extraBits = magW - 1 - manW

    val (shifted, zex0, zsgn) = if (nearOne) {
      // the region of the binade e holds |log2 x| scaled up by 2^(e-2), and
      // its value with `zlead` leading zeros is a fraction in [1/2, 1) times
      // 2^-(zlead + 1), so the result is that fraction times
      // 2^(2 - e - zlead)
      val zlead = if (bit(fracW - 1, zres) == 1) {
        0
      } else if (bit(fracW - 2, zres) == 1) {
        1
      } else {
        2
      }
      (slice(0, fracW, zres << zlead) << exW,
       exBias + 1 - binade - zlead,
       if (belowOne) { 1 } else { 0 })
    } else {
      // log2 x = E + log2(1 + t), and the region read at the binade of t
      // holds log2(1 + t) scaled up by 2^(e-2). A power of two leaves the
      // table nothing to add.
      val frac = if (whole) {
        SafeLong(0)
      } else {
        slice(0, fracW, (zres << 2) >> binade)
      }

      //  x >= 1 (E >= 0):  log2(x) =    E      + frac
      //  x <  1 (E <  0):  log2(x) = -((|E|-1) + (1 - frac))
      //
      // A frac of zero leaves 1 - frac at one, an increment of the integer
      // part.
      val exNobias = x.ex - exBias
      val zsgn     = if (exNobias < 0) { 1 } else { 0 }
      val zint0    = if (exNobias >= 0) {
        SafeLong(exNobias)
      } else {
        SafeLong(-exNobias - 1)
      }

      val (zint, zfrac) = if (zsgn == 0) {
        (zint0, frac)
      } else if (frac == SafeLong(0)) {
        (zint0 + 1, SafeLong(0))
      } else {
        (zint0, (SafeLong(1) << fracW) - frac)
      }

      val mag = (zint << fracW) + zfrac
      assert(mag > 0, f"log2 magnitude is zero at x = ${x.toDouble}")

      val count = magW - binaryWidthSL(mag)
      (mag << count, exBias + exW - 1 - count, zsgn)
    }

    // -----------------------------------------------------------------------

    val zRounded             = roundBySpec(RoundSpec.roundToEven, extraBits, shifted)
    val zMoreThan2AfterRound = bit(manW + 1, zRounded)

    val zMan = slice(0, manW, zRounded)
    val zEx  = zex0 + zMoreThan2AfterRound

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
