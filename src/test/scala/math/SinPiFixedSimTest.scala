package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._
import rial.arith._

import scala.math._
import scala.util.Random

//
// sinPiFixed and cosPiFixed keep their precision relative, not absolute: the
// table holds sin(pi y) of the reduced argument y, but one region per binade
// of y, each scaled by a power of two, so the result stays good to a few ulp
// however close it is to zero -- at a tiny argument, and next to every x at
// which the sine or the cosine vanishes. What limits it is the address: a
// region is split into 16 sub-intervals and the position inside one carries
// inBitW bits, one fewer than the mantissa of the reduced argument leaves.
// The tests below measure the error over the binades of y and print it; what
// they assert of every binade is the tolerance below.
//
class SinPiFixedSimTest extends AnyFunSpec {

  val spec = RealSpec.Float32Spec
  val r    = new Random(20260911)

  /** how many arguments each binade is measured over. */
  val n = 20000

  /** the error the dropped position bit and the table leave, in ulp. */
  val tolerance = 4

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

  /** cos(pi x), reduced so that a half-odd-integral x gives zero. */
  private def cosPi(x: Double): Double = {
    sinPi(x + 0.5)
  }

  /** the distance between two results, in ulp. */
  private def ulp(z: RealGeneric, ref: RealGeneric): Double = {
    (z.value - ref.value).abs.toDouble
  }

  private def sinPiFixed(x: RealGeneric): RealGeneric = {
    SinPiFixedSimGeneric(true, cfg.coeffs(SinPiFixed), x)
  }

  private def cosPiFixed(x: RealGeneric): RealGeneric = {
    SinPiFixedSimGeneric(false, cfg.coeffs(CosPiFixed), x)
  }

  /** a random argument whose distance from `edge` lies in [2^k, 2^(k+1)). */
  private def inBand(edge: Double, k: Int): RealGeneric = {
    val d = pow(2.0, k) * (1.0 + r.nextDouble())
    fl32(if (r.nextBoolean()) { edge + d } else { edge - d })
  }

  /** Measure and print the error over the binades of the reduced argument,
   *  which is the distance of the argument from `edge`.
   *
   *  @param edge  the argument at which the function vanishes
   *  @param f     the function under test
   *  @param exact its reference
   */
  private def bands(edge: Double, f: RealGeneric => RealGeneric,
    exact: Double => Double): Unit = {
    println(f"  binade   band of y              max ulp   mean ulp   " +
            f"max relative error")
    for (k <- -24 to -2) {
      var maxUlp = 0.0
      var sumUlp = 0.0
      var maxRel = 0.0
      for (i <- 1 to n) {
        val x   = inBand(edge, k)
        val ref = exact(x.toDouble)
        val z   = f(x)

        assert(!z.isNaN && !z.isInfinite && abs(z.toDouble) <= 1.0,
          f"x = ${x.toDouble}, z = ${z.toDouble} is not a sine")
        assert(ulp(z, fl32(ref)) <= tolerance,
          f"x = ${x.toDouble}, z = ${z.toDouble}, ref = ${ref}, " +
          f"error = ${ulp(z, fl32(ref))}%.0f ulp")

        maxUlp = max(maxUlp, ulp(z, fl32(ref)))
        sumUlp += ulp(z, fl32(ref))
        if (ref != 0.0) {
          maxRel = max(maxRel, abs(z.toDouble - ref) / abs(ref))
        }
      }
      println(f"  e = ${-k}%2d   [2^$k%3d, 2^${k + 1}%3d)   ${maxUlp}%10.0f " +
              f"${sumUlp / n}%10.2f ${maxRel}%20.3e")
    }
  }

  // --------------------------------------------------------------------------

  it("the error of sinPiFixed(x) over the binades of a small |x|") {
    println(f"sinPiFixed(x) against sin(pi x), ${n}%d arguments per binade of x")
    bands(0.0, sinPiFixed(_), sinPi(_))
  }

  it("the error of sinPiFixed(x) over the binades of the distance from 1") {
    println(f"sinPiFixed(x) against sin(pi x), ${n}%d arguments per binade " +
            f"of 1 - x, where the result vanishes")
    bands(1.0, sinPiFixed(_), sinPi(_))
  }

  it("the error of cosPiFixed(x) over the binades of the distance from 1/2") {
    println(f"cosPiFixed(x) against cos(pi x), ${n}%d arguments per binade " +
            f"of x - 1/2, where the result vanishes")
    bands(0.5, cosPiFixed(_), cosPi(_))
  }

  it("the error of sinPiFixed(x) and cosPiFixed(x) over the whole period") {
    println(f"sinPiFixed(x) and cosPiFixed(x) against sin/cos(pi x), " +
            f"${n}%d arguments in [-2, 2)")
    var maxUlpSin = 0.0
    var maxUlpCos = 0.0
    for (i <- 1 to n) {
      val x = fl32(4.0 * (r.nextDouble() - 0.5))

      val refSin = fl32(sinPi(x.toDouble))
      val zSin   = sinPiFixed(x)
      assert(ulp(zSin, refSin) <= tolerance,
        f"x = ${x.toDouble}, sinPiFixed = ${zSin.toDouble}, " +
        f"ref = ${refSin.toDouble}, error = ${ulp(zSin, refSin)}%.0f ulp")

      val refCos = fl32(cosPi(x.toDouble))
      val zCos   = cosPiFixed(x)
      assert(ulp(zCos, refCos) <= tolerance,
        f"x = ${x.toDouble}, cosPiFixed = ${zCos.toDouble}, " +
        f"ref = ${refCos.toDouble}, error = ${ulp(zCos, refCos)}%.0f ulp")

      maxUlpSin = max(maxUlpSin, ulp(zSin, refSin))
      maxUlpCos = max(maxUlpCos, ulp(zCos, refCos))
    }
    println(f"  sinPiFixed       ${maxUlpSin}%10.0f")
    println(f"  cosPiFixed       ${maxUlpCos}%10.0f")
  }

  it("the special values of sinPiFixed(x) and cosPiFixed(x)") {
    assert(sinPiFixed(fl32(+0.0)).value == fl32(+0.0).value)
    assert(sinPiFixed(fl32(-0.0)).value == fl32(-0.0).value)
    assert(sinPiFixed(fl32(+0.5)).value == fl32(+1.0).value)
    assert(sinPiFixed(fl32(-0.5)).value == fl32(-1.0).value)
    assert(sinPiFixed(fl32(+1.5)).value == fl32(-1.0).value)
    assert(sinPiFixed(fl32(Double.NaN)).isNaN)
    assert(sinPiFixed(fl32(Double.PositiveInfinity)).isNaN)
    assert(sinPiFixed(fl32(Double.NegativeInfinity)).isNaN)

    assert(cosPiFixed(fl32(+0.0)).value == fl32(+1.0).value)
    assert(cosPiFixed(fl32(-0.0)).value == fl32(+1.0).value)
    assert(cosPiFixed(fl32(+1.0)).value == fl32(-1.0).value)
    assert(cosPiFixed(fl32(Double.NaN)).isNaN)
  }
}
