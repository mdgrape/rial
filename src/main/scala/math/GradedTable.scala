//% @file GradedTable.scala
//
// A table that gives a region of its own to each binade of its argument.
//
// A table indexed by a fixed point argument loses the small arguments: the
// sub-intervals nearest to zero are the ones that cover the widest range of
// the argument's own magnitude, which is where a function that vanishes there
// needs them narrowest. What the layout below does instead is to take the
// argument as a float and give each of its binades a region of the table,
// divided into the same number of sub-intervals as every other region, so
// that the sub-interval is always a fixed fraction of the argument itself.
// A binade over which the function bends too much for one region is split
// further, by the top bits of the mantissa; the binades past the last region
// share it, which is exact when the function has settled into a straight line
// by then.
//
// Each region stores its values scaled by a power of two, chosen so that
// every region fills the same range: the result assembly then needs only a
// short leading-zero count and undoes the scale in the exponent.
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math._
import java.lang.Math.scalb

import spire.math.SafeLong
import spire.implicits._

import chisel3.util.log2Ceil

import rial.table._
import rial.util.ScalaUtil._

/** One binade of the argument of a graded table.
 *
 *  @param binade the binade `e`, i.e. the arguments in `[2^-e, 2^-(e-1))`
 *  @param splitW the binade is divided into `2^splitW` regions, chosen by
 *                that many top bits of the mantissa
 */
case class GradedBinade(binade: Int, splitW: Int)

/** One region of a graded table: the `top`-th `2^splitW`-th of the binade. */
case class GradedRegion(binade: Int, splitW: Int, top: Int) {

  /** The mantissas the region covers, as a range of `[0, 1)`. */
  def manRange: (Double, Double) = {
    val width = scalb(1.0, -splitW)
    (top * width, (top + 1) * width)
  }
}

/** The regions of a graded table, in the order the address puts them in.
 *
 *  The address is the region index above a sub-interval index of `subW` bits.
 *  The region index is a lookup on the binade of the argument plus, for a
 *  binade that is split, the top bits of the mantissa; the sub-interval index
 *  and the position inside it are the bits below those. A binade outside the
 *  list is answered by the nearest region the list has.
 *
 *  @param binades the binades the table holds, in increasing order
 *  @param subW    sub-intervals of a region, as a power of two
 */
class GradedLayout(val binades: Seq[GradedBinade], val subW: Int) {

  require(binades.nonEmpty, "a graded table holds at least one binade")
  require(binades.map(_.binade) == binades.map(_.binade).sorted.distinct,
    "the binades of a graded table must be listed in increasing order, " +
    "each of them once")

  /** The regions, in address order. */
  val regions: Seq[GradedRegion] = binades.flatMap(b => {
    (0 until (1 << b.splitW)).map(t => GradedRegion(b.binade, b.splitW, t))
  })

  val nRegion = regions.length

  /** Width of the region index the address carries above the sub-interval. */
  val regionW = log2Ceil(nRegion)

  /** Width of the whole address. */
  val adrW = regionW + subW

  /** The smallest binade the table holds; a smaller one reads its region. */
  val minBinade = binades.head.binade

  /** The largest binade the table holds; a larger one reads its region. */
  val maxBinade = binades.last.binade

  /** The most bits of the mantissa any binade is split by. */
  val splitMaxW = binades.map(_.splitW).max

  /** The binade whose region answers `e`. */
  def binadeOf(e: Int): Int = {
    min(max(e, minBinade), maxBinade)
  }

  /** The first region of the binade that answers `e`, indexed by `e` itself.
   *  The hardware reads it as a lookup on the binade. */
  val baseOf: Seq[Int] = (0 to maxBinade).map(e => {
    regions.indexWhere(_.binade == binadeOf(e))
  })

  /** How many top bits of the mantissa choose the region within the binade
   *  that answers `e`, indexed by `e` itself. */
  val splitWOf: Seq[Int] = (0 to maxBinade).map(e => {
    binades.find(_.binade == binadeOf(e)).get.splitW
  })

  /** The region that answers an argument of binade `e` and mantissa `man`.
   *
   *  @param e    binade of the argument
   *  @param man  mantissa of the argument
   *  @param manW width of the mantissa
   */
  def regionIndex(e: Int, man: SafeLong, manW: Int): Int = {
    val splitW = splitWOf(binadeOf(e))
    val top    = if (splitW == 0) {
      0
    } else {
      slice(manW - splitW, splitW, man).toInt
    }
    baseOf(binadeOf(e)) + top
  }

  /** The fraction the sub-interval index and the position inside it are taken
   *  from: the mantissa with the bits that chose the region shifted out.
   *
   *  @param e    binade of the argument
   *  @param man  mantissa of the argument
   *  @param manW width of the mantissa
   */
  def key(e: Int, man: SafeLong, manW: Int): SafeLong = {
    slice(0, manW, man << splitWOf(binadeOf(e)))
  }

  /** Half of the span of the mantissas the address and the position read as
   *  one, and zero when they read every mantissa of `manW` bits as itself.
   *
   *  They keep the top bits of the mantissa and drop the rest, so the
   *  mantissas they read as one run from the grid point up to one unit of
   *  their own last place short of the next. A polynomial generated at the
   *  argument plus this much is centred on them, so that dropping those bits
   *  costs as much either way instead of always the same way; a mantissa the
   *  grid keeps whole is read as itself, and asks for no offset at all.
   *
   *  @param r      the region the mantissa falls in
   *  @param manW   bits the mantissa of that region's argument carries
   *  @param inBitW width of the position inside a sub-interval
   */
  def gridOffset(r: GradedRegion, manW: Int, inBitW: Int): Double = {
    val keptW = r.splitW + subW + inBitW
    if (manW > keptW) {
      (scalb(1.0, -keptW) - scalb(1.0, -manW)) * 0.5
    } else {
      0.0
    }
  }

  /** The sub-intervals of every region, in address order, of the value `f`
   *  gives for a region and a mantissa of its binade.
   *
   *  They are one table, so that the widths its coefficients are stored at
   *  are the ones the whole table needs.
   *
   *  @param nOrder order of the polynomial of a sub-interval
   *  @param manWOf bits the mantissa of a region's argument carries, which a
   *                caller that shifts the leading one of its argument out
   *                gives one region at a time
   *  @param inBitW width of the position inside a sub-interval
   *  @param f      the value a region stores for a mantissa in `[0, 1)`
   */
  def funcTable(nOrder: Int, manWOf: GradedRegion => Int, inBitW: Int)(
    f: (GradedRegion, Double) => Double): FuncTable = {
    val table = new FuncTable(nOrder)
    regions.foreach(r => {
      val (manMin, manMax) = r.manRange
      val offset = gridOffset(r, manWOf(r), inBitW)
      table.addEqualInterval((m: Double) => { f(r, m + offset) }, 1 << subW,
        manMin, manMax)
    })
    table
  }

  /** The regions in address order, zero-filled to the power of two the
   *  address is. */
  def generateTable(gen: PolynomialGenerator,
    bitW: PolynomialBitWidth): PolynomialCoefficients = {
    val cs    = gen.generateTable(bitW).coeffs
    val zeros = Seq.fill(bitW.nOrder + 1)(0L)
    val nRow  = 1 << adrW
    require(cs.length <= nRow,
      f"${cs.length}%d table entries do not fit the ${nRow}%d the address has")
    PolynomialCoefficients(bitW, cs ++ Seq.fill(nRow - cs.length)(zeros))
  }
}
