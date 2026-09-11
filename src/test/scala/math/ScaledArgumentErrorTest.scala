package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._
import rial.arith._

import scala.math._
import scala.util.Random

//
// Not an assertion test: it measures what scaling the argument costs the
// caller.  Over bands of |x| it compares
//
//   (a) exp(x) and sin(x) as the functions that scale the argument themselves
//   (b) pow2(x * log2(e)) and sinPi(x / pi), the product rounded to FP32 as a
//       caller would have it
//
// against exp(x) / sin(x) in double, in ulp of the result, and prints the
// relative error of both around the multiples of pi, where the rounding of
// x / pi costs the most.
//
// Run with:  sbt "testOnly rial.tests.ScaledArgumentErrorTest"
//
class ScaledArgumentErrorTest extends AnyFunSpec {

  val spec = RealSpec.Float32Spec
  val r    = new Random(20260911)

  /** how many arguments each band is measured over. */
  val n = 20000

  val cfg       = MathFuncConfig.standardFP32
  val scaledCfg = MathFuncConfig.standardScaledFP32

  val log2e = 1.0 / log(2.0)

  private def fl32(x: Double): RealGeneric = {
    new RealGeneric(spec, x)
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

  /** whether a result is one the distance in ulp says something about. */
  private def measurable(z: RealGeneric): Boolean = {
    !z.isZero && !z.isInfinite && !z.isNaN && !z.isSubnormal
  }

  private def band(k: Int, error: RealGeneric => Option[(Double, Double)]) = {
    var maxA = 0.0
    var maxB = 0.0
    var sumA = 0.0
    var sumB = 0.0
    var used = 0
    for (i <- 1 to n) {
      error(inBand(k)).foreach { case (a, b) => {
        maxA = max(maxA, a)
        maxB = max(maxB, b)
        sumA += a
        sumB += b
        used += 1
      }}
    }
    val div = max(1, used)
    println(f"  [2^$k%3d, 2^${k + 1}%3d)  ${maxA}%10.0f ${sumA / div}%10.3f " +
            f"${maxB}%10.0f ${sumB / div}%10.3f  ${used}%7d")
  }

  // --------------------------------------------------------------------------

  it("the error of exp(x) and of pow2(x * log2(e)) over the bands of |x|") {
    println(f"exp(x) against pow2(fl32(x * log2(e))), ${n}%d arguments per band")
    println("  band            max ulp   mean ulp    max ulp   mean ulp  " +
            "samples")
    println("                    exp        exp        pow2       pow2")
    for (k <- -4 to 6) {
      band(k, x => {
        val ref = fl32(exp(x.toDouble))
        if (measurable(ref)) {
          val a = ulp(ExpSimGeneric(cfg.coeffs(Exp), x), ref)
          val b = ulp(Pow2SimGeneric(scaledCfg.coeffs(Pow2),
                                     fl32(x.toDouble * log2e)), ref)
          Some((a, b))
        } else {
          None
        }
      })
    }
  }

  it("the error of sin(x) and of sinPi(x / pi) over the bands of |x|") {
    println(f"sin(x) against sinPi(fl32(x / pi)), ${n}%d arguments per band")
    println("  band            max ulp   mean ulp    max ulp   mean ulp  " +
            "samples")
    println("                    sin        sin       sinPi      sinPi")
    for (k <- -4 to 7) {
      band(k, x => {
        val ref = fl32(sin(x.toDouble))
        if (measurable(ref)) {
          val a = ulp(SinCosSimGeneric(true, cfg.coeffs(Sin), x), ref)
          val b = ulp(SinPiSimGeneric(true, scaledCfg.coeffs(SinPi),
                                      fl32(x.toDouble / Pi)), ref)
          Some((a, b))
        } else {
          None
        }
      })
    }
  }

  it("the relative error of sin(x) and of sinPi(x / pi) around the multiples of pi") {
    println(f"sin(x) against sinPi(fl32(x / pi)) within 2^-10 of k pi, " +
            f"${n}%d arguments each")
    println("    k     max relative error  max relative error")
    println("                    sin               sinPi")
    for (k <- 1 to 4) {
      var maxA = 0.0
      var maxB = 0.0
      for (i <- 1 to n) {
        val x     = fl32(k * Pi + (2.0 * r.nextDouble() - 1.0) * pow(2.0, -10))
        val exact = sin(x.toDouble)
        if (exact != 0.0) {
          val a = SinCosSimGeneric(true, cfg.coeffs(Sin), x)
          val b = SinPiSimGeneric(true, scaledCfg.coeffs(SinPi),
                                  fl32(x.toDouble / Pi))
          maxA = max(maxA, abs(a.toDouble - exact) / abs(exact))
          maxB = max(maxB, abs(b.toDouble - exact) / abs(exact))
        }
      }
      println(f"  ${k}%3d ${maxA}%22.3e ${maxB}%19.3e")
    }
  }
}
