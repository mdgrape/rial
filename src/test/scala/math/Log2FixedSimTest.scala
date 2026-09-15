package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._
import rial.arith._

import spire.math.SafeLong
import spire.implicits._

import scala.math._
import scala.util.Random

//
// log2Fixed keeps its relative precision on both sides of x = 1: its table
// holds |log2 x| itself over the binades of the distance of x from one, so
// the result carries the same number of significant bits however close x
// comes to it. The tests below measure that over the bands of |x - 1| on
// either side, and over the binades of x away from one, where the exponent
// of the argument carries the integer part of the result.
//
class Log2FixedSimTest extends AnyFunSpec {

  val spec = RealSpec.Float32Spec
  val r    = new Random(20260911)

  /** how many arguments each band of |x - 1| is measured over. */
  val n = 20000

  /** how many arguments each binade of x is measured over. */
  val nBinade = 2000

  /** the tolerance the result is held to, in ulp. */
  val tolerance = 2

  /** the binades of x the bands are printed for; every one is measured. */
  val printedBinade = 20

  val cfg = MathFuncConfig.standardFixedFP32

  private def fl32(x: Double): RealGeneric = {
    new RealGeneric(spec, x)
  }

  private def log2(x: Double): Double = {
    log(x) / log(2.0)
  }

  /** the distance between two results, in ulp. */
  private def ulp(z: RealGeneric, ref: RealGeneric): Double = {
    (z.value - ref.value).abs.toDouble
  }

  private def log2Fixed(x: RealGeneric): RealGeneric = {
    Log2FixedSimGeneric(cfg.coeffs(Log2Fixed), x)
  }

  /** the band of |x - 1| the argument falls in: |x - 1| in [2^-k, 2^-k+1)
   *  has the band k. */
  private def bandOf(x: RealGeneric): Int = {
    val d = abs(1.0 - x.toDouble)
    if (d == 0.0) {
      0
    } else {
      ceil(-log2(d)).toInt
    }
  }

  /** The errors log2Fixed makes over one band of arguments. */
  private class Band {

    var count  = 0
    var maxUlp = 0.0
    var sumUlp = 0.0
    var maxRel = 0.0

    def add(x: RealGeneric): Unit = {
      val exact = log2(x.toDouble)
      val z     = log2Fixed(x)

      assert(!z.isNaN && !z.isInfinite,
        f"x = ${x.toDouble}, log2Fixed = ${z.toDouble} is not finite")
      assert(z.toDouble <= 0.0 || x.toDouble > 1.0,
        f"x = ${x.toDouble} <= 1, but log2Fixed = ${z.toDouble} is positive")
      assert(z.toDouble >= 0.0 || x.toDouble < 1.0,
        f"x = ${x.toDouble} >= 1, but log2Fixed = ${z.toDouble} is negative")

      val err = ulp(z, fl32(exact))
      count  += 1
      maxUlp  = max(maxUlp, err)
      sumUlp += err
      if (exact != 0.0) {
        maxRel = max(maxRel, abs(z.toDouble - exact) / abs(exact))
      }
    }

    def print(label: String): Unit = {
      println(f"  ${label}%-16s ${count}%8d ${maxUlp}%8.0f " +
              f"${sumUlp / count}%10.2f ${maxRel}%13.3e")
    }
  }

  private val header =
    "  band                    n  max ulp   mean ulp  max rel err"

  // --------------------------------------------------------------------------

  it("the error of log2Fixed(x) over the bands of x - 1 above one") {
    val bands = Array.fill(spec.manW + 2)(new Band)
    for (k <- 1 to spec.manW) {
      for (i <- 1 to n) {
        val x = fl32(1.0 + pow(2.0, -k) * (1.0 + r.nextDouble()))
        val b = bandOf(x)
        if (b >= 1 && b < bands.length) {
          bands(b).add(x)
        }
      }
    }

    println(f"log2Fixed(x) against log(x)/log(2) for x - 1 in each band")
    println(header)
    for (k <- 1 until bands.length) {
      if (bands(k).count > 0) {
        bands(k).print(f"[2^${-k}%3d, 2^${1 - k}%3d)")
      }
    }
    for (k <- 1 until bands.length) {
      assert(bands(k).maxUlp <= tolerance,
        f"x - 1 in [2^${-k}%d, 2^${1 - k}%d): " +
        f"${bands(k).maxUlp}%.0f ulp is more than the ${tolerance}%d allowed")
    }
  }

  it("the error of log2Fixed(x) over the bands of 1 - x below one") {
    val bands = Array.fill(spec.manW + 3)(new Band)
    for (k <- 2 to spec.manW + 1) {
      for (i <- 1 to n) {
        val x = fl32(1.0 - pow(2.0, -k) * (1.0 + r.nextDouble()))
        val b = bandOf(x)
        if (b >= 2 && b < bands.length) {
          bands(b).add(x)
        }
      }
    }

    println(f"log2Fixed(x) against log(x)/log(2) for 1 - x in each band")
    println(header)
    for (k <- 2 until bands.length) {
      if (bands(k).count > 0) {
        bands(k).print(f"[2^${-k}%3d, 2^${1 - k}%3d)")
      }
    }
    for (k <- 2 until bands.length) {
      assert(bands(k).maxUlp <= tolerance,
        f"1 - x in [2^${-k}%d, 2^${1 - k}%d): " +
        f"${bands(k).maxUlp}%.0f ulp is more than the ${tolerance}%d allowed")
    }
  }

  it("the error of log2Fixed(x) over the binades of x outside [1/2, 2)") {
    println(f"log2Fixed(x) against log(x)/log(2) in each binade of x")
    println(header)

    var maxUlpAll = 0.0
    for (e <- spec.exBias to (1 - spec.exBias) by -1) {
      // |log2 x| >= 1 everywhere but in [1/2, 2), which the bands of the
      // distance from one cover instead
      if (e != 0 && e != -1) {
        val band = new Band
        for (i <- 1 to nBinade) {
          band.add(new RealGeneric(spec, 0, e + spec.exBias,
            SafeLong(BigInt(spec.manW, r))))
        }
        if (abs(e) <= printedBinade) {
          band.print(f"[2^$e%4d, 2^${e + 1}%4d)")
        }
        assert(band.maxUlp <= tolerance,
          f"x in [2^$e%d, 2^${e + 1}%d): ${band.maxUlp}%.0f ulp is more " +
          f"than the ${tolerance}%d allowed")
        maxUlpAll = max(maxUlpAll, band.maxUlp)
      }
    }
    println(f"  the largest error over every binade: ${maxUlpAll}%.0f ulp")
  }

  it("the error of log2Fixed(x) at the exact powers of two") {
    for (e <- spec.exBias to (1 - spec.exBias) by -1) {
      val x   = fl32(pow(2.0, e))
      val z   = log2Fixed(x)
      val ref = fl32(e.toDouble)
      assert(ulp(z, ref) == 0.0,
        f"x = 2^${e}%d, log2Fixed = ${z.toDouble}, ref = ${ref.toDouble}, " +
        f"error = ${ulp(z, ref)}%.0f ulp")
    }
    assert(log2Fixed(fl32(1.0)).value == fl32(0.0).value,
      "log2Fixed(1) is not +0")
  }
}
