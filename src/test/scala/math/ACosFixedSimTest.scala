package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._
import rial.arith._

import scala.math._
import scala.util.Random

//
// acosFixed keeps its precision relative, not absolute: its table holds the
// angle itself, but one region per binade of 1 - |x|, each scaled by a power
// of two, so the result stays good to a few ulp as it approaches zero next to
// x = 1. What limits it is the address: a region is split into 32
// sub-intervals and the position inside one carries inBitW bits, three fewer
// than the mantissa of the argument leaves, which costs about 2^-21 of the
// result. The tests below measure the error over the bands of 1 - x and of
// 1 + x and print it; what they assert of every band is the tolerance below.
//
class ACosFixedSimTest extends AnyFunSpec {

  val spec = RealSpec.Float32Spec
  val r    = new Random(20260911)

  /** how many arguments each band is measured over. */
  val n = 20000

  /** the error the dropped position bits and the table leave, in ulp. */
  val tolerance = 8

  val cfg = MathFuncConfig.standardFixedFP32

  private def fl32(x: Double): RealGeneric = {
    new RealGeneric(spec, x)
  }

  /** the distance between two results, in ulp. */
  private def ulp(z: RealGeneric, ref: RealGeneric): Double = {
    (z.value - ref.value).abs.toDouble
  }

  private def acosFixed(x: RealGeneric): RealGeneric = {
    ACosFixedSimGeneric(cfg.coeffs(ACosFixed), x)
  }

  /** a random argument whose distance from `edge` lies in [2^k, 2^(k+1)). */
  private def inBand(edge: Double, k: Int): RealGeneric = {
    val d = pow(2.0, k) * (1.0 + r.nextDouble())
    fl32(edge - signum(edge) * d)
  }

  /** measure and print the error over the bands of the distance from `edge`. */
  private def bands(edge: Double): Unit = {
    println(f"  band              max ulp   mean ulp   max relative error")
    for (k <- -24 to -1) {
      var maxUlp = 0.0
      var sumUlp = 0.0
      var maxRel = 0.0
      for (i <- 1 to n) {
        val x     = inBand(edge, k)
        val exact = acos(x.toDouble)
        val z     = acosFixed(x)

        assert(!z.isNaN && !z.isInfinite && 0.0 <= z.toDouble && z.toDouble <= Pi,
          f"x = ${x.toDouble}, acos = ${z.toDouble} is not an angle of acos")
        assert(ulp(z, fl32(exact)) <= tolerance,
          f"x = ${x.toDouble}, acos = ${z.toDouble}, ref = ${exact}, " +
          f"error = ${ulp(z, fl32(exact))}%.0f ulp")

        maxUlp = max(maxUlp, ulp(z, fl32(exact)))
        sumUlp += ulp(z, fl32(exact))
        if (exact != 0.0) {
          maxRel = max(maxRel, abs(z.toDouble - exact) / abs(exact))
        }
      }
      println(f"  [2^$k%3d, 2^${k + 1}%3d)   ${maxUlp}%10.0f ${sumUlp / n}%10.1f " +
              f"${maxRel}%20.3e")
    }
  }

  // --------------------------------------------------------------------------

  it("the error of acos(x) over the bands of 1 - x, where the result vanishes") {
    println(f"acos(x) against math.acos, ${n}%d arguments per band of 1 - x")
    bands(1.0)
  }

  it("the error of acos(x) over the bands of 1 + x, where the result is pi") {
    println(f"acos(x) against math.acos, ${n}%d arguments per band of 1 + x")
    bands(-1.0)
  }

  it("the error of acos(x) where the argument is below a half") {
    println(f"acos(x) against math.acos, ${n}%d arguments of |x| < 1/2")
    println(f"  band              max ulp   mean ulp   max relative error")
    var maxUlp = 0.0
    var sumUlp = 0.0
    var maxRel = 0.0
    for (i <- 1 to n) {
      val x     = fl32(r.nextDouble() - 0.5)
      val exact = acos(x.toDouble)
      val z     = acosFixed(x)

      assert(ulp(z, fl32(exact)) <= tolerance,
        f"x = ${x.toDouble}, acos = ${z.toDouble}, ref = ${exact}, " +
        f"error = ${ulp(z, fl32(exact))}%.0f ulp")

      maxUlp = max(maxUlp, ulp(z, fl32(exact)))
      sumUlp += ulp(z, fl32(exact))
      maxRel = max(maxRel, abs(z.toDouble - exact) / abs(exact))
    }
    println(f"  |x| < 1/2        ${maxUlp}%10.0f ${sumUlp / n}%10.1f " +
            f"${maxRel}%20.3e")
  }

  it("the special values of acos(x)") {
    assert(acosFixed(fl32(+1.0)).value == fl32(0.0).value)
    assert(acosFixed(fl32(-1.0)).value == fl32(Pi).value)
    assert(acosFixed(fl32(+0.0)).value == fl32(Pi * 0.5).value)
    assert(acosFixed(fl32(-0.0)).value == fl32(Pi * 0.5).value)
    assert(acosFixed(fl32(Double.NaN)).isNaN)

    // a finite argument beyond an end gets the value of that end, so an
    // argument that overshoots it by rounding needs no clamp
    for (x <- Seq(1.5, 2.0, 1e30)) {
      assert(acosFixed(fl32( x)).value == fl32(0.0).value, f"acos(${x}) is not 0")
      assert(acosFixed(fl32(-x)).value == fl32(Pi).value,  f"acos(${-x}) is not pi")
    }
    assert(acosFixed(fl32(Double.PositiveInfinity)).isNaN)
    assert(acosFixed(fl32(Double.NegativeInfinity)).isNaN)
  }
}
