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

class ATan2Phase1SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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

  def errorLSB(x: RealGeneric, y: Double): Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm + x.spec.manW)
  }

  def atan2Test(tRec: PolynomialCoefficients, spec: RealSpec, n: Int, r: Random,
                generatorStr: String,
                generatorX:      (RealSpec, Random) => RealGeneric,
                generatorYoverX: (RealSpec, Random) => RealGeneric,
                tolerance: Int) = {
    test(s"atan2Phase1(y, x), format ${spec.toStringShort}, ${generatorStr}") {
      var maxError    = 0.0
      var xatMaxError = 0.0
      var yatMaxError = 0.0
      var zatMaxError = 0.0

      var err1lsbPos = 0
      var err1lsbNeg = 0
      var err2lsbPos = 0
      var err2lsbNeg = 0
      var errNlsbPos = 0
      var errNlsbNeg = 0

      for (i <- 1 to n) {
        val yOverX = generatorYoverX(spec, r)
        val x      = generatorX(spec, r)
        val x0     = x.toDouble
        val y      = new RealGeneric(spec, x0 * yOverX.toDouble)
        val y0     = y.toDouble

        val z0  = min(abs(y0), abs(x0)) * (1.0 / max(abs(y0), abs(x0)))
        val z0r = new RealGeneric(spec, z0)
        val zres = ATan2Phase1SimGeneric(tRec, y, x)

        // Decode NaN-boxed special values into comparable Reals
        val zi = if (zres.ex == maskI(spec.exW)) {
          val sgnCoef = if (zres.sgn == 0) 1 else -1
          if      (zres.man == SafeLong(1)) new RealGeneric(spec, 0, 0, 0)            //  0
          else if (zres.man == SafeLong(2)) new RealGeneric(spec, 1, 0, 0)            // -0
          else if (zres.man == SafeLong(3)) new RealGeneric(spec, 0, 0, 0)            //  0
          else if (zres.man == SafeLong(4)) new RealGeneric(spec, sgnCoef *  1.0)     // |y/x| = 1
          else if (zres.man == SafeLong(5)) new RealGeneric(spec, sgnCoef * -1.0)     // |y/x| = 1
          else                              zres
        } else {
          // normal value: strip status bits (exponent MSB, mantissa LSB)
          new RealGeneric(
            spec, 0,
            zres.ex  & maskI(spec.exW - 1),
            zres.man & (maskI(spec.manW - 1) << 1)
          )
        }

        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble)

        if (x0.isInfinity) {
          // |y/x| ~ 0 (encoded), the decoded value is +0
          assert(zi.isZero || zi.isNaN)
        } else if (x0.isNaN) {
          // NaN-boxed
          assert(zres.ex == maskI(spec.exW))
        } else {
          if (erri.abs > tolerance || erri.isNaN) {
            println(f"Error more than ${tolerance} LSB : x = ${x.toDouble}%14.7e, y = ${y.toDouble}%14.7e : refz = $z0%14.7e simz = ${zi.toDouble}%14.7e $errf%14.7e $erri%f")
            val xsgn = bit(spec.W - 1, x.value).toInt
            val xexp = slice(spec.manW, spec.exW, x.value)
            val xman = x.value & maskSL(spec.manW)
            val ztestsgn = bit(spec.W - 1, zi.value).toInt
            val ztestexp = slice(spec.manW, spec.exW, zi.value)
            val ztestman = zi.value & maskSL(spec.manW)
            val zrefsgn = bit(spec.W - 1, z0r.value).toInt
            val zrefexp = slice(spec.manW, spec.exW, z0r.value)
            val zrefman = z0r.value & maskSL(spec.manW)
            println(f"test: x   = ${x0}(${x.sgn}|${x.ex}(${x.ex - x.spec.exBias})|${x.man.toLong.toBinaryString})")
            println(f"test: y   = ${y0}(${y.sgn}|${y.ex}(${y.ex - y.spec.exBias})|${y.man.toLong.toBinaryString})")
            println(f"test: y/x = ${min(abs(x0), abs(y0)) / max(abs(x0), abs(y0))}")
            println(f"test: ref = ${z0}")
            println(f"test: sim = ${zd}(zres = ${zres.sgn}|${zres.ex}|${zres.man.toLong.toBinaryString})")
            println(f"test: test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman.toLong.toBinaryString}) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString})")
            if      (erri >  2.0) errNlsbPos += 1
            else if (erri < -2.0) errNlsbNeg += 1
            else if (erri >  1.0) err2lsbPos += 1
            else                  err2lsbNeg += 1
          } else if (erri >=  1.0) err1lsbPos += 1
          else      if (erri <= -1.0) err1lsbNeg += 1

          assert(erri.abs <= tolerance.toDouble)

          if (maxError < erri.abs) {
            maxError    = erri.abs
            xatMaxError = x0
            yatMaxError = y0
            zatMaxError = zd
          }
        }
      }
      println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
            + f"${zatMaxError} != ${min(abs(xatMaxError), abs(yatMaxError)) / max(abs(xatMaxError), abs(yatMaxError))}, "
            + f"diff = ${zatMaxError - min(abs(xatMaxError), abs(yatMaxError)) / max(abs(xatMaxError), abs(yatMaxError))}")
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
      println(f"N=$n%d : 2LSB errors positive $err2lsbPos%d / negative $err2lsbNeg%d")
      println(f"N=$n%d : 2<   errors positive $errNlsbPos%d / negative $errNlsbNeg%d")
    }
  }

  private def mkRecTable(nOrder: Int, addrW: Int, manW: Int, extraBits: Int): PolynomialCoefficients = {
    val t = new ReciprocalTable(nOrder, addrW, manW - addrW, manW + extraBits)
    t.generateTable(t.getBitWidth())
  }

  // --------------------------------------------------------------------------
  // FP32

  val recFP32 = mkRecTable(2, 8, RealSpec.Float32Spec.manW, 3)

  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within y/x >  2^24",                generateRealWithin(-1.0, 1.0, _, _), generateRealWithin( pow(2.0,  24),  pow(2.0, 128), _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within y/x < -2^24",                generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(-pow(2.0, 128), -pow(2.0,  24), _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within y/x >  2^24, with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin( pow(2.0,  24),  pow(2.0, 128), _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within y/x < -2^24, with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin(-pow(2.0, 128), -pow(2.0,  24), _, _), 3)

  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within y/x >  2^12",                generateRealWithin(-1.0, 1.0, _, _), generateRealWithin( pow(2.0,  12),  pow(2.0,  24), _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within y/x < -2^12",                generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(-pow(2.0,  24), -pow(2.0,  12), _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within y/x >  2^12, with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin( pow(2.0,  12),  pow(2.0,  24), _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within y/x < -2^12, with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin(-pow(2.0,  24), -pow(2.0,  12), _, _), 3)

  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within 1 < y/x <  2^12",            generateRealWithin(-1.0, 1.0, _, _), generateRealWithin( 1.0,             pow(2.0,   8), _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within 1 > y/x > -2^12",            generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(-pow(2.0,   8), -1.0,             _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within 1 < y/x <  2^12, large x",   generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin( 1.0,             pow(2.0,   8), _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within 1 > y/x > -2^12, large x",   generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin(-pow(2.0,   8), -1.0,             _, _), 3)

  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within  2^-12 < y/x <  1",          generateRealWithin(-1.0, 1.0, _, _), generateRealWithin( pow(2.0, -12),  1.0,             _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within -2^-12 > y/x > -1",          generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(-1.0,           -pow(2.0, -12),  _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within  2^-12 < y/x <  1, large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin( pow(2.0, -12),  1.0,             _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within -2^-12 > y/x > -1, large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin(-1.0,           -pow(2.0, -12),  _, _), 3)

  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within 0 < y/x <  2^-12",           generateRealWithin(-1.0, 1.0, _, _), generateRealWithin( 0.0,             pow(2.0, -12),  _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within 0 > y/x > -2^-12",           generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(-pow(2.0, -12),  0.0,             _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within 0 < y/x <  2^-12, large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin( 0.0,             pow(2.0, -12),  _, _), 3)
  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within 0 > y/x > -2^-12, large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100), _, _), generateRealWithin(-pow(2.0, -12),  0.0,             _, _), 3)

  atan2Test(recFP32, RealSpec.Float32Spec, n, r, "Test Within y/x ~ 1",                    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin( 0.999998, 1.000002, _, _), 3)

  // --------------------------------------------------------------------------
  // BFloat16

  val recBF16 = mkRecTable(0, 7, RealSpec.BFloat16Spec.manW, 1)

  atan2Test(recBF16, RealSpec.BFloat16Spec, n, r, "Test Within y/x > 2^24",      generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, 24),  pow(2.0, 128), _, _), 2)
  atan2Test(recBF16, RealSpec.BFloat16Spec, n, r, "Test Within y/x > 2^12",      generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, 12),  pow(2.0,  24), _, _), 2)
  atan2Test(recBF16, RealSpec.BFloat16Spec, n, r, "Test Within 1 < y/x < 2^12",  generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(1.0,           pow(2.0,   8), _, _), 2)
  atan2Test(recBF16, RealSpec.BFloat16Spec, n, r, "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, -12), 1.0,           _, _), 2)
  atan2Test(recBF16, RealSpec.BFloat16Spec, n, r, "Test Within y/x < 2^-12",     generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(0.0,           pow(2.0, -12), _, _), 2)
}
