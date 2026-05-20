package rial.tests.funcs

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

class ACosPhase2SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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

  def acosTest(tSqrt: PolynomialCoefficients, tACos: PolynomialCoefficients,
               spec: RealSpec, n: Int, r: Random,
               generatorStr: String,
               generator: (RealSpec, Random) => RealGeneric,
               tolerance: Int) = {
    test(s"acos(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0
      var maxError    = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0
      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for (i <- 1 to n) {
        val x  = generator(spec, r)
        val x0 = x.toDouble

        val xref = if (x0 > 1.0 && !x.isNaN && !x.isInfinite) 1.0
                   else if (x0 < -1.0 && !x.isNaN && !x.isInfinite) -1.0
                   else x0

        val z0  = acos(xref)
        val z0r = new RealGeneric(spec, z0)

        val s1 = ACosPhase1SimGeneric(tSqrt, x)
        val zi = ACosPhase2SimGeneric(tACos, s1)
        val zd = zi.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toInt

        if (x.isNaN) {
          assert(zi.isNaN)
        } else {
          if (erri.abs > tolerance) {
            println(f"x=${x0}, ref=${z0}, sim=${zd}, errLSB=${erri}, s1=(${s1.sgn},${s1.ex},${s1.man})")
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
        println(f"N=$n%d : largest errors ${maxError.toInt}%d at x=${xatMaxError} (ref=${acos(xatMaxError)})")
      }
      for (kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- $k errors positive $errPos%d / negative $errNeg%d")
      }
      println("---------------------------------------------------------------")
    }
  }

  private def mkTables(nOrder: Int, addrW: Int, manW: Int, extraBits: Int):
      (PolynomialCoefficients, PolynomialCoefficients) = {
    val sqrtT = new SqrtTable       (nOrder, addrW, manW - addrW, manW + extraBits)
    val acosT = new ACosPhase2Table (nOrder, addrW, manW - addrW, manW + extraBits)
    (sqrtT.generateTable(sqrtT.getBitWidth()),
     acosT.generateTable(acosT.getBitWidth()))
  }

  val (sqrtFP32, acosFP32) = mkTables(2, 8, RealSpec.Float32Spec.manW, 3)
  val (sqrtBF16, acosBF16) = mkTables(0, 7, RealSpec.BFloat16Spec.manW, 1)

  acosTest(sqrtFP32, acosFP32, RealSpec.Float32Spec, n, r,
    "Test 0      < x <= 2^-4",   generateRealWithin(0.0,                pow(2.0, -4),       _, _), 3)
  acosTest(sqrtFP32, acosFP32, RealSpec.Float32Spec, n, r,
    "Test 2^-4   < x <= 1-2^-4", generateRealWithin(pow(2.0, -4),       1.0 - pow(2.0, -4), _, _), 3)
  acosTest(sqrtFP32, acosFP32, RealSpec.Float32Spec, n, r,
    "Test 1-2^-4 < x <= 1",      generateRealWithin(1.0 - pow(2.0, -4), 1.0,                _, _), 3)
  acosTest(sqrtFP32, acosFP32, RealSpec.Float32Spec, n, r,
    "Test 0      >= x > -2^-4",    generateRealWithin(-pow(2.0, -4),     0.0,                _, _), 3)
  acosTest(sqrtFP32, acosFP32, RealSpec.Float32Spec, n, r,
    "Test -2^-4   >= x > -1+2^-4", generateRealWithin(-1.0 + pow(2.0, -4), -pow(2.0, -4),     _, _), 3)
  acosTest(sqrtFP32, acosFP32, RealSpec.Float32Spec, n, r,
    "Test -1+2^-4 >= x > -1",      generateRealWithin(-1.0,              -1.0 + pow(2.0, -4), _, _), 3)
  acosTest(sqrtFP32, acosFP32, RealSpec.Float32Spec, n, r,
    "Test special values", generateSpecialValues(_, _), 1)

  acosTest(sqrtBF16, acosBF16, RealSpec.BFloat16Spec, n, r,
    "Test 0      < x <= 2^-4",   generateRealWithin(0.0,                pow(2.0, -4),       _, _), 3)
  acosTest(sqrtBF16, acosBF16, RealSpec.BFloat16Spec, n, r,
    "Test 2^-4   < x <= 1-2^-4", generateRealWithin(pow(2.0, -4),       1.0 - pow(2.0, -4), _, _), 3)
  acosTest(sqrtBF16, acosBF16, RealSpec.BFloat16Spec, n, r,
    "Test 1-2^-4 < x <= 1",      generateRealWithin(1.0 - pow(2.0, -4), 1.0,                _, _), 3)
  acosTest(sqrtBF16, acosBF16, RealSpec.BFloat16Spec, n, r,
    "Test 0      >= x > -2^-4",    generateRealWithin(-pow(2.0, -4),     0.0,                _, _), 3)
  acosTest(sqrtBF16, acosBF16, RealSpec.BFloat16Spec, n, r,
    "Test -2^-4   >= x > -1+2^-4", generateRealWithin(-1.0 + pow(2.0, -4), -pow(2.0, -4),     _, _), 3)
  acosTest(sqrtBF16, acosBF16, RealSpec.BFloat16Spec, n, r,
    "Test -1+2^-4 >= x > -1",      generateRealWithin(-1.0,              -1.0 + pow(2.0, -4), _, _), 3)
  acosTest(sqrtBF16, acosBF16, RealSpec.BFloat16Spec, n, r,
    "Test special values", generateSpecialValues(_, _), 1)
}
