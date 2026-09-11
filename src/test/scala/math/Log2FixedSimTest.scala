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
// The precision of log2Fixed is absolute, not relative: its table holds
// log2(1 + t) itself, so the fraction of the result carries a fixed number of
// bits whatever the result is. Where |log2 x| >= 1 the integer part carries
// the precision and the result is good to the ulp; as x approaches one the
// result approaches zero while its absolute error stays ~2^-outBitW, so its
// relative error grows. The first test below measures that growth over the
// bands of |x - 1| and prints it; it asserts nothing there beyond the result
// being finite and of the right sign.
//
class Log2FixedSimTest extends AnyFunSpec {

  val spec = RealSpec.Float32Spec
  val r    = new Random(20260911)

  /** how many arguments each band of |x - 1| is measured over. */
  val n = 20000

  /** how many arguments each binade of x is measured over. */
  val nBinade = 2000

  /** the tolerance LogSimTest holds log to, in ulp. */
  val tolerance = 3

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

  /** a random argument with |x - 1| in [2^-k, 2^-k+1), on either side of 1. */
  private def aroundOne(k: Int): RealGeneric = {
    val d = pow(2.0, -k) * (1.0 + r.nextDouble())
    fl32(if (r.nextBoolean()) { 1.0 + d } else { 1.0 - d })
  }

  /** a random argument with x in [2^e, 2^(e+1)). */
  private def inBinade(e: Int): RealGeneric = {
    new RealGeneric(spec, 0, e + spec.exBias, SafeLong(BigInt(spec.manW, r)))
  }

  // --------------------------------------------------------------------------

  it("the error of log2Fixed(x) over the bands of |x - 1|") {
    println(f"log2Fixed(x) against log(x)/log(2), ${n}%d arguments per band")
    println("  |x - 1|           max ulp   mean ulp   max relative error")
    for (k <- 24 to 1 by -1) {
      var maxUlp = 0.0
      var sumUlp = 0.0
      var maxRel = 0.0
      for (i <- 1 to n) {
        val x     = aroundOne(k)
        val exact = log2(x.toDouble)
        val z     = log2Fixed(x)

        assert(!z.isNaN && !z.isInfinite,
          f"x = ${x.toDouble}, log2Fixed = ${z.toDouble} is not finite")
        assert(z.toDouble <= 0.0 || x.toDouble > 1.0,
          f"x = ${x.toDouble} <= 1, but log2Fixed = ${z.toDouble} is positive")
        assert(z.toDouble >= 0.0 || x.toDouble < 1.0,
          f"x = ${x.toDouble} >= 1, but log2Fixed = ${z.toDouble} is negative")

        maxUlp = max(maxUlp, ulp(z, fl32(exact)))
        sumUlp += ulp(z, fl32(exact))
        if (exact != 0.0) {
          maxRel = max(maxRel, abs(z.toDouble - exact) / abs(exact))
        }
      }
      println(f"  [2^${-k}%3d, 2^${1 - k}%3d)  ${maxUlp}%10.0f ${sumUlp / n}%10.1f " +
              f"${maxRel}%20.3e")
    }
  }

  it("the error of log2Fixed(x) over the binades of x where |log2 x| >= 1") {
    println(f"log2Fixed(x) against log(x)/log(2), ${nBinade}%d arguments per binade")
    println("  x                 max ulp   mean ulp")
    var maxUlpAll = 0.0
    for (e <- 127 to -126 by -1) {
      // |log2 x| >= 1 everywhere in the binade except [0.5, 2), which the
      // bands of |x - 1| cover instead.
      if (e != 0 && e != -1) {
        var maxUlp = 0.0
        var sumUlp = 0.0
        for (i <- 1 to nBinade) {
          val x   = inBinade(e)
          val ref = fl32(log2(x.toDouble))
          val z   = log2Fixed(x)

          assert(ulp(z, ref) <= tolerance,
            f"x = ${x.toDouble}, log2Fixed = ${z.toDouble}, " +
            f"ref = ${ref.toDouble}, error = ${ulp(z, ref)}%.0f ulp")

          maxUlp = max(maxUlp, ulp(z, ref))
          sumUlp += ulp(z, ref)
        }
        maxUlpAll = max(maxUlpAll, maxUlp)
        println(f"  [2^$e%4d, 2^${e + 1}%4d)  ${maxUlp}%8.0f ${sumUlp / nBinade}%10.2f")
      }
    }
    println(f"  the largest error over every binade: ${maxUlpAll}%.0f ulp")
  }

  it("the error of log2Fixed(x) at the exact powers of two") {
    for (e <- 127 to -126 by -1) {
      val x = fl32(pow(2.0, e))
      val z = log2Fixed(x)
      val ref = fl32(e.toDouble)
      assert(ulp(z, ref) <= tolerance,
        f"x = 2^${e}%d, log2Fixed = ${z.toDouble}, ref = ${ref.toDouble}, " +
        f"error = ${ulp(z, ref)}%.0f ulp")
    }
    assert(log2Fixed(fl32(1.0)).value == fl32(0.0).value,
      "log2Fixed(1) is not +0")
  }
}
