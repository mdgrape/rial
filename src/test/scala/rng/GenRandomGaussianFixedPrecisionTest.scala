package rial.tests

import org.scalatest.funspec.AnyFunSpec

import org.apache.commons.math3.special.Erf

import rial.fixed.FixedSpec
import rial.rng._
import rial.table.PolynomialEval

import scala.math._

// Precision of the table of GenRandomGaussianFixed in units of the output LSB
// (2^-fracW), measured on the table's own rows rather than on random words:
// which row a word reads is tested on its own in
// GenRandomGaussianFixedSimMatchTest.
//
// Every sub-interval of every row is evaluated at the positions a Chebyshev
// fit peaks at -- around the extrema of T_(nOrder+1) -- plus a uniform grid,
// and rounded as the generator rounds it. The reference is the capped |z| the
// row approximates, at the x the position stands for in the row's own
// convention: position 0 is the xMin end of the interval, t = -1.
object GenRandomGaussianFixedPrecision {

  val nOrder = 3
  val sqrt2  = sqrt(2.0)

  // |z| of the central page, erfInv polished with two Newton steps on erf
  def probitCentral(x: Double): Double = {
    var z = sqrt2 * Erf.erfInv(x)
    for (_ <- 0 until 2) {
      val g  = Erf.erf(z / sqrt2) - x
      val gp = sqrt(2.0 / Pi) * exp(-0.5 * z * z)
      z = z - g / gp
    }
    z
  }

  /** the positions of a sub-interval the error is sampled at: eight around
   *  each extremum of T_(nOrder+1) and a grid of sixteen, and both ends */
  def positions(inBitW: Int): Seq[Int] = {
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

  /** worst signed error, mean error, and the largest magnitude output, of one
   *  region of the table */
  case class Stat(worst: Double, mean: Double, maxOut: Double)

  val regionName = Seq("central", "octaves", "last")

  /** the error of each region of the table, in units of 2^-fracW */
  def measure(rndW: Int, zSpec: FixedSpec, splitW: Int, guardBits: Int,
    posW: Option[Int] = None): Seq[Stat] = {
    val positionW = GenRandomGaussianFixed.positionWidth(zSpec, posW)
    val nOctaves  = GenRandomGaussianFixed.octaves(rndW, zSpec, posW)
    val table     = GenRandomGaussianFixed.table(zSpec, nOrder, nOctaves,
      positionW, splitW, guardBits)
    val inBitW    = table.inBitW
    val nSub      = 1 << splitW
    val eps       = pow(2.0, -zSpec.fracW)
    val dxs       = positions(inBitW)

    val worst  = Array.fill(3)(0.0)
    val sum    = Array.fill(3)(0.0)
    val count  = Array.fill(3)(0L)
    val maxOut = Array.fill(3)(0.0)

    for (idx <- table.ft.intervals.indices) {
      val interval = table.ft.intervals(idx)
      val row      = idx / nSub
      val region   = if (row == 0) { 0 } else if (row <= nOctaves) { 1 } else { 2 }
      for (dx <- dxs) {
        val t = 2.0 * dx / (1 << inBitW) - 1.0
        val x = (t + 1.0) * 0.5 * (interval.xMax - interval.xMin) + interval.xMin

        val polyOut = PolynomialEval.clamp(
          PolynomialEval.eval(idx, BigInt(dx), table.coeffs),
          table.coeffs.bitW, clampNegToZero = true)
        val mag = GenRandomGaussianFixed.roundGuard(polyOut, zSpec, guardBits)

        val zExact = if (row == 0) {
          min(probitCentral(x), table.zCap)
        } else {
          table.exact(row, x)
        }
        val err = mag.toDouble - zExact / eps
        if (abs(err) > abs(worst(region))) { worst(region) = err }
        sum(region)    += err
        count(region)  += 1
        maxOut(region)  = max(maxOut(region), mag.toDouble * eps)
      }
    }
    (0 until 3).map(r => Stat(worst(r), sum(r) / count(r), maxOut(r)))
  }
}

class GenRandomGaussianFixedPrecisionTest extends AnyFunSpec {
  import GenRandomGaussianFixedPrecision._

  private def check(rndW: Int, zSpec: FixedSpec, bound: Double): Unit = {
    it(f"Q${zSpec.intW}.${zSpec.fracW} at rndW=$rndW is within $bound%.1f LSB of the exact probit") {
      val stats = measure(rndW, zSpec, 3, 3)
      stats.zip(regionName).foreach({ case (s, name) => {
        info(f"$name%-8s max ${s.worst}%+7.3f mean ${s.mean}%+7.4f LSB")
      }})
      assert(stats.map(s => abs(s.worst)).max <= bound,
        "max error per region: " + stats.map(_.worst).mkString(", "))
    }
  }

  check(64, FixedSpec(16, 12), 4.0)
  check(64, FixedSpec(16, 13), 4.0)
  check(32, FixedSpec(16, 13), 4.0)
  check(64, FixedSpec(32, 24), 4.0)

  it("no sampled position goes beyond the reach, and the reach is the cap of the last row") {
    for ((rndW, zSpec) <- Seq((64, FixedSpec(16, 12)), (32, FixedSpec(16, 13)),
                              (32, FixedSpec(16, 12)), (64, FixedSpec(32, 24)))) {
      val reach = GenRandomGaussianFixed.reach(rndW, zSpec)
      val out   = measure(rndW, zSpec, 3, 3).map(_.maxOut).max
      val nOct  = GenRandomGaussianFixed.octaves(rndW, zSpec)
      val cap   = GenRandomGaussianFixed.table(zSpec, nOrder, nOct,
        GenRandomGaussianFixed.positionWidth(zSpec), 3, 3).lastCap
      val eps   = pow(2.0, -zSpec.fracW)
      info(f"Q${zSpec.intW}.${zSpec.fracW} rndW=$rndW: largest sampled $out%.7f, " +
        f"reach $reach%.7f, cap of the last row $cap%.7f")
      assert(out <= reach)
      assert(abs(reach - cap) <= 4 * eps)
    }
  }
}
