package rial.tests

import org.scalatest.funspec.AnyFunSpec

import org.apache.commons.math3.special.Erf

import rial.arith.RealSpec
import rial.rng._

import scala.math._

// ULP precision of GenRandomGaussianSim (the bit-exact model of
// GenRandomGaussian) against the exact probit, over the windows that carry the
// largest error: the central page, the shallowest tail octave, a deep-ish
// individual octave and the first deep-tail octave.
//
// The reference is the inverse CDF of the x the input bits stand for, with the
// table's own convention inside an interval (position 0 is the xMin, shallow,
// end and x runs linearly in the position). The sampled positions are those
// where a Chebyshev fit peaks -- around the extrema of T_(nOrder+1) of every
// sub-interval -- plus a uniform grid.
class GenGaussianPrecisionTest extends AnyFunSpec {

  val spec   = RealSpec.Float32Spec
  val manW   = spec.manW                  // 23
  val rndW   = 64
  val nOrder = 3
  val lump   = 2
  val mainOctaves = rndW - 2 - manW       // 39
  val one   = 1 << manW                   // the position of x = xMax
  val page1 = BigInt(1) << (rndW - 2)
  val sqrt2 = sqrt(2.0)

  // |z| of the tail page, solved exactly as the table itself solves it. The
  // solver keeps state, so the reference gets a (tiny) table of its own.
  private val ref = new ProbitTable(spec, 1, 4, 8, 8, lump, 0)
  private def probitHalf(x: Double): Double = ref.probitHalf(x)

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

  /** one window of the input space: the inputs that read one interval */
  private case class Window(name: String, gran: Int,
    rndOf: Int => BigInt, xOf: Int => Double, refOf: Double => Double)

  /** the central page at one octave of x: the position drops the lowest
   *  xex+1 bits of the mantissa, so the mantissa itself is the sample
   */
  private def central(xex: Int) = Window(f"central xex=$xex", 1 << (xex + 1),
    (dx: Int) => (BigInt(1) << (manW + xex)) | BigInt(dx),
    (dx: Int) => (one | dx).toDouble * pow(2.0, -(manW + 2 + xex)),
    probitCentral)

  /** the main tail octave k: x in [2^-(k+1), 2^-(k+2)], the mantissa is the
   *  position itself
   */
  private def main(k: Int) = Window(f"main xex=$k", 1,
    (dx: Int) => page1 | (BigInt(1) << (manW + k)) | BigInt(dx),
    (dx: Int) => {
      val xMin = pow(2.0, -(k + 1))
      xMin + (pow(2.0, -(k + 2)) - xMin) * dx / one
    },
    probitHalf)

  /** the deep tail octave c (xex0 == 0): the leading one of the mantissa sits
   *  at manW-1-c and the position is what is left below it
   */
  private def deep(c: Int) = Window(f"deep clz=$c", 1 << (c + 1),
    (dx: Int) => page1 | (BigInt(1) << (manW - 1 - c)) | BigInt(dx >> (c + 1)),
    (dx: Int) => {
      val addr = mainOctaves + 1 + c
      val xMin = pow(2.0, -addr)
      xMin + (pow(2.0, -(addr + 1)) - xMin) * dx / one
    },
    probitHalf)

  /** the positions a Chebyshev fit of `splitW` sub-intervals peaks at, rounded
   *  down to the grid the window can be read at
   */
  private def positions(splitW: Int, gran: Int): Seq[Int] = {
    val nSub = 1 << splitW
    val len  = 1 << (manW - splitW)
    (0 until nSub).flatMap(s => {
      val clusters = (0 to nOrder + 1).flatMap(k => {
        val o  = round((cos(k * Pi / (nOrder + 1)) + 1.0) * 0.5 * len).toInt
        val lo = max(0, min(len - 8, o - 4))
        lo until lo + 8
      })
      val grid = (0 until 16).map(j => min(len - 1, round((j + 0.5) * len / 16.0).toInt))
      (clusters ++ grid).map(o => ((s * len + o) / gran) * gran)
    }).distinct.filter(_ < one)
  }

  /** the largest |ULP| error of the window at that configuration */
  private def maxUlp(w: Window, splitW: Int, guardBits: Int): Long = {
    var worst = 0L
    for (dx <- positions(splitW, w.gran)) {
      val z   = GenRandomGaussianSim.genGaussian(w.rndOf(dx), rndW, spec, nOrder,
        lump, splitW, guardBits).value.toBigInt.toInt
      val err = (z & 0x7fffffff).toLong - bitsOf(w.refOf(w.xOf(dx))).toLong
      if (abs(err) > abs(worst)) { worst = err }
    }
    worst
  }

  private val windows = Seq(central(0), main(0), main(19), deep(0))

  private def report(splitW: Int, guardBits: Int): Seq[Long] = {
    windows.map(w => {
      val e = maxUlp(w, splitW, guardBits)
      info(f"splitW=$splitW guardBits=$guardBits ${w.name}%-14s max ${e}%+8d ULP")
      e
    })
  }

  it("splitW = 3 with 3 guard bits is within 3 ULP of the exact probit") {
    val errs = report(3, 3)
    assert(errs.map(abs(_)).max <= 3, "max |ULP| per window: " + errs.mkString(", "))
  }

  it("one polynomial per interval and no guard bits is thousands of ULP off") {
    val errs = report(0, 0)
    // the shallowest tail octave and the central page are the worst two; the
    // assertion documents the error the split and the guard bits remove
    assert(errs.map(abs(_)).max >= 1000, "max |ULP| per window: " + errs.mkString(", "))
  }
}
