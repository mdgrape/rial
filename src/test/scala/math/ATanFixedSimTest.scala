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
// The precision of atanFixed is absolute, not relative: its table holds
// atan(t) itself, so the ulp error of the result grows as the result
// approaches zero. Below the linear threshold atan(t) is t, which the module
// returns unchanged, so the error falls back to zero there; the first test
// below measures the growth over the bands of |t| and prints it. The
// tolerance of ATan2Phase2SimTest is asserted where the result is of the
// order of one and in the linear range, twice that where the result is a
// quarter of one, and the bands in between are only held to a finite result
// of the right sign.
//
class ATanFixedSimTest extends AnyFunSpec {

  val spec = RealSpec.Float32Spec
  val r    = new Random(20260911)

  /** how many arguments each band is measured over. */
  val n = 20000

  /** the tolerance ATan2Phase2SimTest holds atan2 to, in ulp. */
  val tolerance = 4

  val cfg = MathFuncConfig.standardFixedFP32

  /** the first band of |t| in which the module returns its argument unchanged. */
  val firstLinearBand = -ATan2Phase2SimGeneric.calcLinearThreshold(spec.manW)

  private def fl32(x: Double): RealGeneric = {
    new RealGeneric(spec, x)
  }

  /** the distance between two results, in ulp. */
  private def ulp(z: RealGeneric, ref: RealGeneric): Double = {
    (z.value - ref.value).abs.toDouble
  }

  /** a random argument with |t| in [2^-k, 2^-k+1), of either sign. */
  private def inBand(k: Int): RealGeneric = {
    val mag = pow(2.0, -k) * (1.0 + r.nextDouble())
    fl32(if (r.nextBoolean()) { mag } else { -mag })
  }

  private def atanFixed(x: RealGeneric): RealGeneric = {
    ATanFixedSimGeneric(cfg.coeffs(ATanFixed), x)
  }

  // --------------------------------------------------------------------------

  it("the error of atanFixed(t) over the bands of |t|") {
    println(f"atanFixed(t) against atan(t), ${n}%d arguments per band")
    println("  band            max ulp   mean ulp   max relative error")
    println("  (the band of 2^0 is |t| >= 1, outside the domain, where the " +
            "result is atan(1))")
    for (k <- 24 to 0 by -1) {
      var maxUlp = 0.0
      var sumUlp = 0.0
      var maxRel = 0.0
      for (i <- 1 to n) {
        val t     = inBand(k)
        val exact = atan(t.toDouble)
        val z     = atanFixed(t)

        assert(!z.isNaN && !z.isInfinite && abs(z.toDouble) <= Pi / 4.0 * (1.0 + 1e-6),
          f"t = ${t.toDouble}, atanFixed = ${z.toDouble} is not an arc tangent " +
          f"of an argument in [-1, 1]")
        assert(z.sgn == t.sgn,
          f"t = ${t.toDouble}, atanFixed = ${z.toDouble} has the wrong sign")

        maxUlp = max(maxUlp, ulp(z, fl32(exact)))
        sumUlp += ulp(z, fl32(exact))
        maxRel = max(maxRel, abs(z.toDouble - exact) / abs(exact))
      }
      println(f"  [2^${-k}%3d, 2^${1 - k}%3d)  ${maxUlp}%10.0f ${sumUlp / n}%10.1f " +
              f"${maxRel}%20.3e")
    }
  }

  it("the error of atanFixed(t) where the result is of the order of one") {
    // |t| in [1/2, 1], where atan(t) is in [pi/8, pi/4) and the quantum of the
    // polynomial output, 2^-outBitW, is half an ulp of the result.
    for (i <- 1 to n) {
      val t   = inBand(1)
      val ref = fl32(atan(t.toDouble))
      val z   = atanFixed(t)
      assert(ulp(z, ref) <= tolerance,
        f"t = ${t.toDouble}, atanFixed = ${z.toDouble}, " +
        f"ref = ${ref.toDouble}, error = ${ulp(z, ref)}%.0f ulp")
    }
    for (i <- 1 to n) {
      val t   = fl32(if (r.nextBoolean()) { 1.0 } else { -1.0 })
      val ref = fl32(atan(t.toDouble))
      val z   = atanFixed(t)
      assert(ulp(z, ref) <= tolerance,
        f"t = ${t.toDouble}, atanFixed = ${z.toDouble}, ref = ${ref.toDouble}")
    }
  }

  it("the error of atanFixed(t) where the result is a quarter of one") {
    // |t| in [1/4, 1/2), where atan(t) falls below 1/4 and the same quantum is
    // a whole ulp of the result, so the same absolute error counts twice.
    for (i <- 1 to n) {
      val t   = inBand(2)
      val ref = fl32(atan(t.toDouble))
      val z   = atanFixed(t)
      assert(ulp(z, ref) <= 2 * tolerance,
        f"t = ${t.toDouble}, atanFixed = ${z.toDouble}, " +
        f"ref = ${ref.toDouble}, error = ${ulp(z, ref)}%.0f ulp")
    }
  }

  it("atanFixed(t) is its argument below the linear threshold") {
    for (k <- 126 to firstLinearBand by -1) {
      for (i <- 1 to 100) {
        val t   = inBand(k)
        val ref = fl32(atan(t.toDouble))
        val z   = atanFixed(t)
        assert(z.value == t.value,
          f"t = ${t.toDouble}, atanFixed = ${z.toDouble} is not the argument")
        assert(ulp(z, ref) <= tolerance,
          f"t = ${t.toDouble}, atanFixed = ${z.toDouble}, " +
          f"ref = ${ref.toDouble}, error = ${ulp(z, ref)}%.0f ulp")
      }
    }
  }

  it("the special values of atanFixed(t)") {
    assert(atanFixed(fl32(+0.0)).value == fl32(+0.0).value)
    assert(atanFixed(fl32(-0.0)).value == fl32(-0.0).value)
    assert(atanFixed(fl32(Double.NaN)).isNaN)

    // |t| >= 1 gives atan(1) = pi/4, which is the accurate value at 1 and a
    // finite, meaningless one above it.
    val quarterPi = fl32(Pi * 0.25)
    for (t <- Seq(1.0, 2.0, 1e30, Double.PositiveInfinity)) {
      assert(atanFixed(fl32( t)).value ==  quarterPi.value,
        f"atanFixed(${t}) = ${atanFixed(fl32(t)).toDouble} != ${quarterPi.toDouble}")
      assert(atanFixed(fl32(-t)).value == (quarterPi.value + (SafeLong(1) << (spec.W - 1))),
        f"atanFixed(${-t}) = ${atanFixed(fl32(-t)).toDouble} != ${-quarterPi.toDouble}")
    }
  }
}
