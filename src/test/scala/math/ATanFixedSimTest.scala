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
// atanFixed keeps its precision relative, not absolute: the table holds
// atan(t) itself, but one region per binade of |t|, each scaled by a power of
// two, so the result stays good to a few ulp however close to zero it is.
// What limits it is the address: a region is divided into 16 sub-intervals
// and the position inside one carries inBitW bits, which leaves a bit of the
// mantissa unread over the binades a whole region covers. Below the binade
// the table ends at, atan(t) is t itself and the module hands the argument
// back. The tests below measure the error over the binades of |t| and print
// it; what they assert of every binade is the tolerance below.
//
class ATanFixedSimTest extends AnyFunSpec {

  val spec = RealSpec.Float32Spec
  val r    = new Random(20260911)

  /** how many arguments each binade is measured over. */
  val n = 20000

  /** the error the dropped position bit and the table leave, in ulp. */
  val tolerance = 4

  val cfg = MathFuncConfig.standardFixedFP32

  /** the first binade of |t| in which the module returns its argument
   *  unchanged. */
  val firstLinearBinade = -ATan2Phase2SimGeneric.calcLinearThreshold(spec.manW)

  private def fl32(x: Double): RealGeneric = {
    new RealGeneric(spec, x)
  }

  /** the distance between two results, in ulp. */
  private def ulp(z: RealGeneric, ref: RealGeneric): Double = {
    (z.value - ref.value).abs.toDouble
  }

  /** a random argument with |t| in [2^-e, 2^(1-e)), of either sign. */
  private def inBinade(e: Int): RealGeneric = {
    val mag = pow(2.0, -e) * (1.0 + r.nextDouble())
    fl32(if (r.nextBoolean()) { mag } else { -mag })
  }

  private def atanFixed(x: RealGeneric): RealGeneric = {
    ATanFixedSimGeneric(cfg.coeffs(ATanFixed), x)
  }

  // --------------------------------------------------------------------------

  it("the error of atanFixed(t) over the binades of |t|") {
    println(f"atanFixed(t) against atan(t), ${n}%d arguments per binade of t")
    println(f"  (the argument is handed back unchanged from the binade " +
            f"${firstLinearBinade}%d on)")
    println(f"  binade   band of |t|            max ulp   mean ulp   " +
            f"max relative error")
    for (e <- firstLinearBinade + 1 to 1 by -1) {
      var maxUlp = 0.0
      var sumUlp = 0.0
      var maxRel = 0.0
      for (i <- 1 to n) {
        val t   = inBinade(e)
        val ref = atan(t.toDouble)
        val z   = atanFixed(t)

        assert(!z.isNaN && !z.isInfinite &&
               abs(z.toDouble) <= Pi / 4.0 * (1.0 + 1e-6),
          f"t = ${t.toDouble}, z = ${z.toDouble} is not an arc tangent of an " +
          f"argument in [-1, 1]")
        assert(z.sgn == t.sgn,
          f"t = ${t.toDouble}, z = ${z.toDouble} has the wrong sign")
        assert(ulp(z, fl32(ref)) <= tolerance,
          f"t = ${t.toDouble}, z = ${z.toDouble}, ref = ${ref}, " +
          f"error = ${ulp(z, fl32(ref))}%.0f ulp")

        maxUlp = max(maxUlp, ulp(z, fl32(ref)))
        sumUlp += ulp(z, fl32(ref))
        maxRel = max(maxRel, abs(z.toDouble - ref) / abs(ref))
      }
      println(f"  e = ${e}%2d   [2^${-e}%3d, 2^${1 - e}%3d)   ${maxUlp}%10.0f " +
              f"${sumUlp / n}%10.2f ${maxRel}%20.3e")
    }
  }

  it("the error of atanFixed(t) at the ends of every binade") {
    // the largest argument of a binade is where the value its region stores
    // comes closest to one, which the polynomial output has to stay below
    val manMax = (SafeLong(1) << spec.manW) - 1
    for (e <- 1 to firstLinearBinade) {
      for (k <- 0 to 255) {
        for (man <- Seq(SafeLong(k), manMax - k)) {
          val t   = new RealGeneric(spec, 0, spec.exBias - e, man)
          val ref = fl32(atan(t.toDouble))
          val z   = atanFixed(t)
          assert(ulp(z, ref) <= tolerance,
            f"t = ${t.toDouble}, atanFixed = ${z.toDouble}, " +
            f"ref = ${ref.toDouble}, error = ${ulp(z, ref)}%.0f ulp")
        }
      }
    }
  }

  it("atanFixed(t) is its argument below the linear threshold") {
    for (e <- 126 to firstLinearBinade by -1) {
      for (i <- 1 to 100) {
        val t   = inBinade(e)
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
