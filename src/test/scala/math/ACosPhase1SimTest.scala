package rial.tests

import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.implicits._

import rial.math._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

class ACosPhase1SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRealWithin(from: Double, to: Double, spec: RealSpec, r: Random) = {
    val rD: Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  var counter = 0
  val specialValues = Seq(
       0.0, -0.0,  1.0, -1.0,
       1.1, -1.1,
       Double.PositiveInfinity, Double.NegativeInfinity,
    )
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

  def acosTest(tSqrt: PolynomialCoefficients, spec: RealSpec, n: Int, r: Random,
               generatorStr: String,
               generator: (RealSpec, Random) => RealGeneric,
               tolerance: Int) = {
    test(s"acosPhase1(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0
      var maxError    = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0
      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for (i <- 1 to n) {
        val x  = generator(spec, r)
        val x0 = x.toDouble

        val xref = if (x0 > 1.0) 1.0 else if (x0 < -1.0) -1.0 else x0
        val z0   = sqrt(1.0 - abs(xref))
        val z0r  = new RealGeneric(spec, z0)

        val zs = ACosPhase1SimGeneric(tSqrt, x)

        // expand NaN-boxed special values back into a comparable Real
        val zi = if (zs.ex == maskI(spec.exW)) {
          if      (zs.man == SafeLong(1)) new RealGeneric(spec, 1.0)
          else if (zs.man == SafeLong(2)) new RealGeneric(spec, 0.0)
          else if (zs.man == SafeLong(3)) new RealGeneric(spec, 0.0)
          else                            new RealGeneric(spec, 0, zs.ex, zs.man)
        } else {
          new RealGeneric(spec, 0, zs.ex, zs.man)
        }
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toInt

        // sign should always carry x.sgn
        if (x.sgn == 1) {
          assert(zs.sgn == 1, f"x(${x0}) < 0, but z.sgn = ${zs.sgn}")
        } else {
          assert(zs.sgn == 0, f"x(${x0}) >= 0, but z.sgn = ${zs.sgn}")
        }

        // NaN-boxing checks for special inputs
        if (x.ex == 0) {
          assert(zs.ex == maskI(spec.exW) && zs.man == SafeLong(1))
        } else if (x.ex >= spec.exBias && x.sgn == 0 && !x.isNaN && !x.isInfinite) {
          assert(zs.ex == maskI(spec.exW) && zs.man == SafeLong(2))
        } else if (x.ex >= spec.exBias && x.sgn == 1 && !x.isNaN && !x.isInfinite) {
          assert(zs.ex == maskI(spec.exW) && zs.man == SafeLong(3))
        }

        if (z0r.isNaN) {
          assert(zi.isNaN)
        } else {
          if (erri.abs > tolerance) {
            println(f"x=${x0}, ref=${z0}, sim=${zd}, errLSB=${erri}")
          }
          if (erri != 0) {
            val k = erri.abs
            if (!errs.contains(k)) errs(k) = (0, 0)
            if (erri >= 0) errs(k) = (errs(k)._1 + 1, errs(k)._2)
            else           errs(k) = (errs(k)._1,     errs(k)._2 + 1)
          }
          assert(erri.abs <= tolerance)
          if (maxError < erri.abs) {
            maxError    = erri.abs
            xatMaxError = x0
            zatMaxError = zd
          }
        }
      }
      println(f"${generatorStr} Summary")
      if (maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError.toInt}%d at x=${xatMaxError}")
      }
      for (kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- $k errors positive $errPos%d / negative $errNeg%d")
      }
      println("---------------------------------------------------------------")
    }
  }

  private def mkSqrtTable(nOrder: Int, addrW: Int, manW: Int, extraBits: Int): PolynomialCoefficients = {
    val t = new SqrtTable(nOrder, addrW, manW - addrW, manW + extraBits)
    t.generateTable(t.getBitWidth())
  }

  val sqrtFP32 = mkSqrtTable(2, 8, RealSpec.Float32Spec.manW, 3)
  val sqrtBF16 = mkSqrtTable(0, 7, RealSpec.BFloat16Spec.manW, 1)

  acosTest(sqrtFP32, RealSpec.Float32Spec, n, r,
    "Test 0      < x <= 2^-4",   generateRealWithin(0.0,                pow(2.0, -4),       _, _), 3)
  acosTest(sqrtFP32, RealSpec.Float32Spec, n, r,
    "Test 2^-4   < x <= 1-2^-4", generateRealWithin(pow(2.0, -4),       1.0 - pow(2.0, -4), _, _), 3)
  acosTest(sqrtFP32, RealSpec.Float32Spec, n, r,
    "Test 1-2^-4 < x <= 1",      generateRealWithin(1.0 - pow(2.0, -4), 1.0,                _, _), 3)
  acosTest(sqrtFP32, RealSpec.Float32Spec, n, r,
    "Test special values", generateSpecialValues(_, _), 1)

  acosTest(sqrtBF16, RealSpec.BFloat16Spec, n, r,
    "Test 0      < x <= 2^-4",   generateRealWithin(0.0,                pow(2.0, -4),       _, _), 3)
  acosTest(sqrtBF16, RealSpec.BFloat16Spec, n, r,
    "Test 2^-4   < x <= 1-2^-4", generateRealWithin(pow(2.0, -4),       1.0 - pow(2.0, -4), _, _), 3)
  acosTest(sqrtBF16, RealSpec.BFloat16Spec, n, r,
    "Test 1-2^-4 < x <= 1",      generateRealWithin(1.0 - pow(2.0, -4), 1.0,                _, _), 3)
  acosTest(sqrtBF16, RealSpec.BFloat16Spec, n, r,
    "Test special values", generateSpecialValues(_, _), 1)
}
