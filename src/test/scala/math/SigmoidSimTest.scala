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

class SigmoidSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  // Random x in [-p/2, p/2] but with full random low-bit precision.
  def generateRealWithin(p: Double, spec: RealSpec, r: Random) = {
    val rD: Double = (r.nextDouble() - 0.5) * p
    val x = new RealGeneric(spec, rD)
    new RealGeneric(spec, (x.value & (maskSL(spec.exW + 1) << spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull(spec: RealSpec, r: Random) = {
    new RealGeneric(spec, SafeLong(BigInt(spec.W, r)))
  }

  var counter = 0
  val specialValues = Seq(
    0.0, -0.0,
    1.0, -1.0,
    2.0, -2.0,
    1e10, -1e10,
  )
  def generateSpecialValues(spec: RealSpec, r: Random) = {
    val idx = counter
    counter += 1
    if (counter >= specialValues.length) { counter = 0 }
    new RealGeneric(spec, specialValues(idx))
  }

  def errorLSB(x: RealGeneric, y: Double): Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm + x.spec.manW)
  }

  def sigmoidTest(
    t: PolynomialCoefficients, spec: RealSpec, n: Int, r: Random,
    generatorStr: String,
    generator: (RealSpec, Random) => RealGeneric,
    tolerance: Double
  ) = {
    test(s"sigmoid(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError    = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val sigmoid = (x: Double) => 1.0 / (1.0 + exp(-x))

      for (i <- 1 to n) {
        val x    = generator(spec, r)
        val x0   = x.toDouble
        val z0   = sigmoid(x0)
        val z0r  = new RealGeneric(spec, z0)
        val zi   = SigmoidSimGeneric(t, x)
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble

        if (errf.abs > tolerance) {
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
                  f"sigmoid(x=${xsgn}|${xexp}(${xexp - spec.exBias})|${xman.toLong.toBinaryString})")
        }

        assert(!z0.isInfinity)
        if (z0.isNaN) {
          assert(x.isNaN)
        } else {
          assert(errf.abs <= tolerance)
          if (maxError < errf.abs) {
            maxError    = errf.abs
            xatMaxError = x0
            zatMaxError = zd
          }
        }
      }
      println(f"${generatorStr} Summary")
      if (maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError} where the value is "
              + f"${zatMaxError} != ${sigmoid(xatMaxError)}, "
              + f"diff = ${zatMaxError - sigmoid(xatMaxError)}, x = ${xatMaxError}")
      }
      println("---------------------------------------------------------------")
    }
  }

  private def mkTable(nOrder: Int, addrW: Int, manW: Int, extraBits: Int): PolynomialCoefficients = {
    val t = new SigmoidTable(nOrder, addrW, manW - addrW, manW + extraBits)
    t.generateTable(t.getBitWidth())
  }

  // --------------------------------------------------------------------------
  // FP32

  val sigmoidF32 = mkTable(2, 8, RealSpec.Float32Spec.manW, 2)

  sigmoidTest(sigmoidF32, RealSpec.Float32Spec, n, r,
    "Test Within (-128,128)", generateRealWithin(128.0, _, _), 1e-5)
  sigmoidTest(sigmoidF32, RealSpec.Float32Spec, n, r,
    "Test All range",         generateRealFull(_, _),          1e-5)
  sigmoidTest(sigmoidF32, RealSpec.Float32Spec, n, r,
    "Test special values",    generateSpecialValues(_, _),     1e-5)

  // --------------------------------------------------------------------------
  // BFloat16

  val sigmoidBF16 = mkTable(0, 7, RealSpec.BFloat16Spec.manW, 0)

  val toleranceBF16 = 0.05 // XXX large tolerance

  sigmoidTest(sigmoidBF16, RealSpec.BFloat16Spec, n, r,
    "Test Within (-128,128)", generateRealWithin(128.0, _, _), toleranceBF16)
  sigmoidTest(sigmoidBF16, RealSpec.BFloat16Spec, n, r,
    "Test All range",         generateRealFull(_, _),          toleranceBF16)
  sigmoidTest(sigmoidBF16, RealSpec.BFloat16Spec, n, r,
    "Test special values",    generateSpecialValues(_, _),     toleranceBF16)
}
