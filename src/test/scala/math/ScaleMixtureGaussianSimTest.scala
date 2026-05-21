package rial.tests

import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.implicits._

import rial.funcs._
import rial.util.ScalaUtil._
import rial.arith._
import rial.approx._

class ScaleMixtureGaussianSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  // Random x in [-p, p] with full random low-bit precision.
  def generateRealWithin(p: Double, spec: RealSpec, r: Random) = {
    val rD: Double = (r.nextDouble() - 0.5) * 2.0 * p
    val x = new RealGeneric(spec, rD)
    new RealGeneric(spec, (x.value & (maskSL(spec.exW + 1) << spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull(spec: RealSpec, r: Random) = {
    new RealGeneric(spec, SafeLong(BigInt(spec.W, r)))
  }

  var counter = 0
  val specialValues = Seq(0.0, -0.0)
  def generateSpecialValues(spec: RealSpec, r: Random) = {
    val idx = counter
    counter += 1
    if (counter >= specialValues.length) counter = 0
    new RealGeneric(spec, specialValues(idx))
  }

  def errorLSB(x: RealGeneric, y: Double): Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm + x.spec.manW)
  }

  def smgTest(t: PolynomialCoefficients, sgmA: Double, sgmB: Double,
              spec: RealSpec, n: Int, r: Random,
              generatorStr: String,
              generator: (RealSpec, Random) => RealGeneric,
              tolerance: Int) = {
    test(s"scaleMixtureGaussian(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError    = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val sgmA2 = sgmA * sgmA
      val sgmB2 = sgmB * sgmB
      val sgmP2 = (sgmA2 * sgmB2) / (sgmA2 - sgmB2)

      println(f"=============================================================")
      println(f"σ_a^2 = ${sgmA2}, σ_b^2 = ${sgmB2}, σ'^2 = ${sgmP2}")

      val refFunc = (w: Double) => {
        val g = exp(-(w * w) / (2 * sgmP2))
        if (g == 0.0) -(w / sgmA2)
        else          -(w / sgmA2) * ((sgmA2 / sgmP2) / (sgmB / (g * sgmA) + 1.0) + 1.0)
      }

      for (i <- 1 to n) {
        val x  = generator(spec, r)
        val x0 = x.toDouble
        val z0 = refFunc(x0)

        val z0r  = new RealGeneric(spec, z0)
        val zi   = ScaleMixtureGaussianSimGeneric(t, x, sgmA, sgmB)
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toLong

        if ((z0r.value - zi.value).abs > tolerance && !z0.isNaN && !z0.isInfinite &&
            (abs(z0) < pow(2.0, spec.exMax + 1))) {
          val xsgn = bit(spec.W - 1, x.value).toInt
          val xexp = slice(spec.manW, spec.exW, x.value)
          val xman = x.value & maskSL(spec.manW)
          val ztestsgn = bit(spec.W - 1, zi.value).toInt
          val ztestexp = slice(spec.manW, spec.exW, zi.value)
          val ztestman = zi.value & maskSL(spec.manW)
          val zrefsgn = bit(spec.W - 1, z0r.value).toInt
          val zrefexp = slice(spec.manW, spec.exW, z0r.value)
          val zrefman = z0r.value & maskSL(spec.manW)
          println(f"x0  = ${x0} (${x.sgn}|${x.ex}|${x.man.toLong.toBinaryString}|)")
          println(f"ref = ${z0}")
          println(f"sim = ${zd}")
          println(f"test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman.toLong.toBinaryString}) != " +
                  f"ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}) = " +
                  f"smg(x = ${xsgn}|${xexp}(${xexp - spec.exBias})|${xman.toLong.toBinaryString})")
        }

        if (x.isInfinite) {
          assert(zi.isInfinite && x.sgn != zi.sgn)
        } else if (x.isNaN) {
          assert(zi.isNaN)
        } else if (x.isZero) {
          assert(zi.isZero)
        } else {
          assert(erri.abs <= tolerance)
          if (maxError < erri.abs) {
            maxError    = erri.abs.toDouble
            xatMaxError = x0
            zatMaxError = zd
          }
        }
      }
      println(f"${generatorStr} Summary")
      if (maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError.toInt}%d ULP where the value is "
              + f"${zatMaxError} != ${refFunc(xatMaxError)}, "
              + f"diff = ${zatMaxError - refFunc(xatMaxError)}, x = ${xatMaxError}")
      }
      println("---------------------------------------------------------------")
    }
  }

  private def mkTable(nOrder: Int, addrW: Int, manW: Int, extraBits: Int,
                      sgmA: Double, sgmB: Double): PolynomialCoefficients = {
    val t = new ScaleMixtureGaussianTable(sgmA, sgmB, nOrder, addrW, manW - addrW, manW + extraBits)
    t.generateTable(t.getBitWidth())
  }

  // --------------------------------------------------------------------------
  // FP32

  val sgmA16 = exp(-1.0)
  val sgmB16 = exp(-6.0)
  val smgMaxBitDigit = ScaleMixtureGaussianTable.tableMaxBitDigit(sgmA16, sgmB16)

  val smgF32 = mkTable(2, 8, RealSpec.Float32Spec.manW, 8, sgmA16, sgmB16)
  val transitionPointF32 =
    pow(2.0, ScaleMixtureGaussianTable.tableDomainDigit(RealSpec.Float32Spec.manW + 8, sgmA16, sgmB16))

  smgTest(smgF32, sgmA16, sgmB16, RealSpec.Float32Spec, n, r,
    f"Test Within (-${transitionPointF32},${transitionPointF32}) where table is used",
    generateRealWithin(transitionPointF32, _, _), pow(2.0, smgMaxBitDigit).toInt)
  smgTest(smgF32, sgmA16, sgmB16, RealSpec.Float32Spec, n, r,
    "Test All range", generateRealFull(_, _), pow(2.0, smgMaxBitDigit).toInt)

  // --------------------------------------------------------------------------
  // BFloat16

  val smgBF16 = mkTable(0, 7, RealSpec.BFloat16Spec.manW, 10, sgmA16, sgmB16)
  val transitionPointBF16 =
    pow(2.0, ScaleMixtureGaussianTable.tableDomainDigit(RealSpec.BFloat16Spec.manW + 10, sgmA16, sgmB16))

  smgTest(smgBF16, sgmA16, sgmB16, RealSpec.BFloat16Spec, n, r,
    f"Test Within (-${transitionPointBF16},${transitionPointBF16}) where table is used",
    generateRealWithin(transitionPointBF16, _, _), pow(2.0, smgMaxBitDigit).toInt)
  smgTest(smgBF16, sgmA16, sgmB16, RealSpec.BFloat16Spec, n, r,
    "Test All range", generateRealFull(_, _), pow(2.0, smgMaxBitDigit).toInt)
}
