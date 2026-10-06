package rial.tests

import org.scalatest.funspec.AnyFunSpec

import org.apache.commons.math3.special.Erf

import rial.arith.RealSpec
import rial.rng._
import rial.table.PolynomialEval

import scala.math._

// ULP precision of the probit table of GenRandomGaussian, measured on the
// table's own rows rather than on random words: which row a word reads is a
// design freedom of the decode, and that mapping is tested on its own in
// GenRandomGaussianSimMatchTest.
//
// Every row of the table (one sub-interval of one interval) is evaluated at
// the positions a Chebyshev fit peaks at -- around the extrema of
// T_(nOrder+1) -- plus a uniform grid, and the result is carried through the
// arithmetic the generator applies after the table. The reference is the
// function the row itself approximates, at the x the position stands for in
// the row's own convention: position 0 is the xMin end, t = -1.
class GenGaussianPrecisionTest extends AnyFunSpec {

  val spec   = RealSpec.Float32Spec
  val manW   = spec.manW                  // 23
  val rndW   = 64
  val nOrder = 3
  val lump   = 2
  val one     = BigInt(1) << manW
  val manMask = one - 1
  val sqrt2   = sqrt(2.0)

  // the intervals the module asks the table for at rndW bits: the central page,
  // one per octave of the tail, and the lumped deepest one
  val mainOctaves = rndW - 2 - manW            // 39
  val nOctaves    = mainOctaves + manW - lump  // 60
  val nIntervals  = nOctaves + 2               // 62

  // |z| of the central page, erfInv polished with two Newton steps on erf
  private def probitCentral(x: Double): Double = {
    var z = sqrt2 * Erf.erfInv(x)
    for (_ <- 0 until 2) {
      val g  = Erf.erf(z / sqrt2) - x
      val gp = sqrt(2.0 / Pi) * exp(-0.5 * z * z)
      z = z - g / gp
    }
    z
  }

  private def bitsOf(z: Double): Int = java.lang.Float.floatToIntBits(z.toFloat)

  /** The generator's arithmetic after the table, copied from GenGaussian.scala
   *  (`polyMan` onwards): the guard bits of the polynomial are rounded away
   *  half up, with the carry into 1.0 clamped; 1.y is multiplied by the FP32
   *  `mul`, the octave's scale or, on the central page, x itself; the product
   *  is truncated to manW bits and normalised by a selection.
   */
  private def afterTable(poly: BigInt, mul: Int, guardBits: Int): Int = {
    val resMax = (BigInt(1) << (manW + guardBits)) - 1
    val res    = poly.max(BigInt(0)).min(resMax)
    val polyMan = if (guardBits == 0) { res } else {
      val rounded = (res >> guardBits) + ((res >> (guardBits - 1)) & 1)
      if (rounded > manMask) { manMask } else { rounded }
    }
    val prod = (one | polyMan) * (one | (BigInt(mul) & manMask))
    val more = ((prod >> (2 * manW + 1)) & 1).toInt
    val zman = (prod >> (manW + more)) & manMask
    ((((BigInt(mul) >> manW) + more) << manW) | zman).toInt
  }

  /** the positions of a row the error is sampled at: eight around each extremum
   *  of T_(nOrder+1) and a grid of sixteen, about fifty in all
   */
  private def positions(inBitW: Int): Seq[Int] = {
    val len = 1 << inBitW
    if (len <= 64) { 0 until len } else {
      val clusters = (0 to nOrder + 1).flatMap(k => {
        val o  = round((cos(k * Pi / (nOrder + 1)) + 1.0) * 0.5 * len).toInt
        val lo = max(0, min(len - 8, o - 4))
        lo until lo + 8
      })
      val grid = (0 until 16).map(j => {
        min(len - 1, round((j + 0.5) * len / 16.0).toInt)
      })
      (clusters ++ grid).distinct.sorted
    }
  }

  private val regionName = Seq("central", "octaves", "lump")

  /** the largest |ULP| error of each region of the table at that configuration */
  private def measure(splitW: Int, guardBits: Int): Seq[Long] = {
    val nSub    = 1 << splitW
    val inBitW  = manW - splitW
    val outBitW = manW + guardBits
    val table   = new ProbitTable(spec, nOrder, nIntervals, inBitW, outBitW,
      lump, splitW)
    val coeffs = table.generateTable(table.getBitWidth())
    val dxs    = positions(inBitW)
    val worst  = Array.fill(3)(0L)

    for (idx <- table.ft.intervals.indices) {
      val interval = table.ft.intervals(idx)
      val addr     = idx / nSub
      val region   = if (addr == 0) { 0 } else if (addr <= nOctaves) { 1 } else { 2 }
      val scale    = table.scaleTable(addr).toInt
      for (dx <- dxs) {
        // the row's own convention, the one FuncTable.addInterval fitted with
        val t = 2.0 * dx / (1 << inBitW) - 1.0
        val x = (t + 1.0) * 0.5 * (interval.xMax - interval.xMin) + interval.xMin
        assert(abs(interval.mapToDomain(x) - t) < 1e-9,
          f"row $idx: mapToDomain disagrees at dx=$dx")

        // x = 0 is the one position with no |z| to compare against
        if (x > 0.0) {
          // the central page multiplies by x, which carries at most manW
          // significant bits and is therefore an exact FP32
          val mul = if (addr == 0) {
            require(x.toFloat.toDouble == x, f"row $idx: x=$x is not an FP32")
            bitsOf(x)
          } else { scale }

          val z = afterTable(PolynomialEval.eval(idx, BigInt(dx), coeffs), mul,
            guardBits)
          val zExact = if (addr == 0) {
            probitCentral(x)
          } else if (addr <= nOctaves) {
            table.probitHalf(x)
          } else {
            min(table.probitHalf(x), table.fMax)
          }

          val err = z.toLong - bitsOf(zExact).toLong
          if (abs(err) > abs(worst(region))) { worst(region) = err }
        }
      }
    }

    worst.toSeq.zip(regionName).map({ case (e, name) => {
      info(f"splitW=$splitW guardBits=$guardBits $name%-8s max ${e}%+8d ULP")
      e
    }})
  }

  it("splitW = 3 with 3 guard bits is within 3 ULP of the exact probit") {
    val errs = measure(3, 3)
    assert(errs.map(abs(_)).max <= 3, "max |ULP| per region: " + errs.mkString(", "))
  }

  it("one polynomial per interval and no guard bits is thousands of ULP off") {
    val errs = measure(0, 0)
    // the assertion documents the error the split and the guard bits remove
    assert(errs.map(abs(_)).max >= 1000, "max |ULP| per region: " + errs.mkString(", "))
  }
}
