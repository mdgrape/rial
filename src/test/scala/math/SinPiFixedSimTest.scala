package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._
import rial.arith._

import scala.math._
import scala.util.Random

//
// The precision of sinPiFixed is absolute, not relative: its table holds
// sin(pi y) itself, so the ulp error of the result grows as the result
// approaches zero. The first test below measures that growth over the bands
// of |x| and prints it; it asserts nothing about the small bands beyond the
// result staying a finite value of a sine. The tolerance of SinPiSimTest is
// asserted only where the result is of the order of one.
//
class SinPiFixedSimTest extends AnyFunSpec {

  val spec = RealSpec.Float32Spec
  val r    = new Random(20260911)

  /** how many arguments each band is measured over. */
  val n = 20000

  /** the tolerance SinPiSimTest holds sinPi to, in ulp. */
  val tolerance = 3

  val cfg = MathFuncConfig.standardFixedFP32

  private def fl32(x: Double): RealGeneric = {
    new RealGeneric(spec, x)
  }

  /** sin(pi x), with the argument reduced so that an integral x gives zero. */
  private def sinPi(x: Double): Double = {
    val k = rint(x)
    val s = if (rint(k / 2.0) * 2.0 == k) { 1.0 } else { -1.0 }
    s * sin(Pi * (x - k))
  }

  /** the distance between two results, in ulp. */
  private def ulp(z: RealGeneric, ref: RealGeneric): Double = {
    (z.value - ref.value).abs.toDouble
  }

  /** a random argument with |x| in [2^k, 2^(k+1)), of either sign. */
  private def inBand(k: Int): RealGeneric = {
    val mag = pow(2.0, k) * (1.0 + r.nextDouble())
    fl32(if (r.nextBoolean()) { mag } else { -mag })
  }

  private def sinPiFixed(x: RealGeneric): RealGeneric = {
    SinPiFixedSimGeneric(true, cfg.coeffs(SinPiFixed), x)
  }

  // --------------------------------------------------------------------------

  it("the error of sinPiFixed(x) over the bands of |x| below 2^-2") {
    println(f"sinPiFixed(x) against sin(pi x), ${n}%d arguments per band")
    println("  band            max ulp   mean ulp   max relative error")
    for (k <- -24 to -2) {
      var maxUlp = 0.0
      var sumUlp = 0.0
      var maxRel = 0.0
      for (i <- 1 to n) {
        val x     = inBand(k)
        val exact = sinPi(x.toDouble)
        val z     = sinPiFixed(x)

        assert(!z.isNaN && !z.isInfinite && abs(z.toDouble) <= 1.0,
          f"x = ${x.toDouble}, sinPiFixed = ${z.toDouble} is not a sine")

        maxUlp = max(maxUlp, ulp(z, fl32(exact)))
        sumUlp += ulp(z, fl32(exact))
        if (exact != 0.0) {
          maxRel = max(maxRel, abs(z.toDouble - exact) / abs(exact))
        }
      }
      println(f"  [2^$k%3d, 2^${k + 1}%3d)  ${maxUlp}%10.0f ${sumUlp / n}%10.1f " +
              f"${maxRel}%20.3e")
    }
  }

  it("the error of sinPiFixed(x) where the result is of the order of one") {
    // The band of |x| does not bound the result: sin(pi x) approaches zero at
    // every integral x as well. The tolerance is asserted where the result
    // itself is of the order of one, which is what an absolute precision of
    // the table resolves to a few ulp of it.
    var used = 0
    for (k <- -2 to 0) {
      for (i <- 1 to n) {
        val x   = inBand(k)
        val ref = fl32(sinPi(x.toDouble))
        if (abs(ref.toDouble) >= 0.25) {
          val z = sinPiFixed(x)
          assert(ulp(z, ref) <= tolerance,
            f"x = ${x.toDouble}, sinPiFixed = ${z.toDouble}, " +
            f"ref = ${ref.toDouble}, error = ${ulp(z, ref)}%.0f ulp")
          used += 1
        }
      }
    }
    assert(used > n, f"only ${used}%d of the arguments were of the order of one")
  }
}
