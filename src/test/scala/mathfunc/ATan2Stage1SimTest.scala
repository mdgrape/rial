//package rial.tests


import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}



//import scopt.OptionParser

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math.ReciprocalSim
import rial.mathfunc._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

class MathFuncATan2Stage1SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRealWithin( from : Double, to : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def atan2Test(t_rec : FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generatorX      : ( (RealSpec, Random) => RealGeneric),
    generatorYoverX : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"atan2(x), format ${spec.toStringShort}, ${generatorStr}") {
      var maxError   = 0.0
      var xatMaxError = 0.0
      var yatMaxError = 0.0
      var zatMaxError = 0.0

      var err1lsbPos = 0
      var err1lsbNeg = 0
      var err2lsbPos = 0
      var err2lsbNeg = 0
      var errNlsbPos = 0
      var errNlsbNeg = 0
      for(i <- 1 to n) {
//         println("------------------")
        val y_over_x = generatorYoverX(spec,r)
        val x  = generatorX(spec,r)
        val x0 = x.toDouble
        val y  = new RealGeneric(spec, x0 * y_over_x.toDouble)
        val y0 = y.toDouble

        val minxy = if(abs(x0) < abs(y0)) {x0} else {y0}
        val maxxy = if(abs(x0) < abs(y0)) {y0} else {x0}

        val z0   = min(abs(y0), abs(x0)) * (1.0 / max(abs(y0), abs(x0)))
        val z0r  = new RealGeneric(spec, z0)
        val zi   = ATan2Stage1Sim.atan2Stage1SimGeneric( t_rec, y, x )._1
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble)

        if (x0.isInfinity) {
          assert(zi.isNaN)
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if (erri.abs>tolerance || erri.isNaN) {
            println(f"Error more than 2 LSB : x = ${x.toDouble}%14.7e, y = ${y.toDouble}%14.7e : refz = $z0%14.7e simz = ${zi.toDouble}%14.7e $errf%14.7e $erri%f")

            val xsgn = bit(spec.W-1, x.value).toInt
            val xexp = slice(spec.manW, spec.exW, x.value)
            val xman = x.value & maskSL(spec.manW)

            val ztestsgn = bit(spec.W-1, zi.value).toInt
            val ztestexp = slice(spec.manW, spec.exW, zi.value)
            val ztestman = zi.value & maskSL(spec.manW)

            val zrefsgn = bit(spec.W-1, z0r.value).toInt
            val zrefexp = slice(spec.manW, spec.exW, z0r.value)
            val zrefman = z0r.value & maskSL(spec.manW)

            println(f"test: x   = ${x0}(${x.sgn}|${x.ex}(${x.ex-x.spec.exBias})|${x.man.toLong.toBinaryString})")
            println(f"test: y   = ${y0}(${y.sgn}|${y.ex}(${y.ex-y.spec.exBias})|${y.man.toLong.toBinaryString})")
            println(f"test: y/x = ${min(abs(x0), abs(y0))/max(abs(x0), abs(y0))}")
            println(f"test: ref = ${z0}")
            println(f"test: sim = ${zd}")
            println(f"test: test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman.toLong.toBinaryString}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")

            if(erri > 2.0) {
              errNlsbPos += 1
            } else if(erri < -2.0) {
              errNlsbNeg += 1
            } else if(erri > 1.0) {
              err2lsbPos += 1
            } else {
              err2lsbNeg += 1
            }
          } else if (erri>=1.0) {
            err1lsbPos+=1
          }
          else if (erri<= -1.0) {
            err1lsbNeg+=1
          }
          assert(erri.abs<=tolerance.toDouble)

          if (maxError < erri.abs) {
            maxError = erri.abs
            xatMaxError = x0
            yatMaxError = y0
            zatMaxError = zd
          }
        }
      }
      println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
            + f"${zatMaxError} != ${min(abs(xatMaxError), abs(yatMaxError))/max(abs(xatMaxError), abs(yatMaxError))}, "
            + f"diff = ${zatMaxError - min(abs(xatMaxError), abs(yatMaxError))/max(abs(xatMaxError), abs(yatMaxError))}")
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
      println(f"N=$n%d : 2LSB errors positive $err2lsbPos%d / negative $err2lsbNeg%d")
      println(f"N=$n%d : 2<   errors positive $errNlsbPos%d / negative $errNlsbNeg%d")
    }
  }

  val nOrderFP32 = 2
  val adrWFP32 = 8
  val extraBitsFP32 = 3
  val atan2FP32ReciprocalTableI = ReciprocalSim.reciprocalTableGeneration(
        nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)

  atan2Test(atan2FP32ReciprocalTableI, RealSpec.Float32Spec, n, r,
    "Test Within y/x > 2^24", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 128),_,_), 2)
  atan2Test(atan2FP32ReciprocalTableI, RealSpec.Float32Spec, n, r,
    "Test Within y/x > 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_), 2)
  atan2Test(atan2FP32ReciprocalTableI, RealSpec.Float32Spec, n, r,
    "Test Within 1 < y/x < 2^12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_), 2)
  atan2Test(atan2FP32ReciprocalTableI, RealSpec.Float32Spec, n, r,
    "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), 1.0,_,_), 2)
  atan2Test(atan2FP32ReciprocalTableI, RealSpec.Float32Spec, n, r,
    "Test Within y/x < 2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_), 2)

  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 1
  val atan2BF16ReciprocalTableI = ReciprocalSim.reciprocalTableGeneration(
        nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)

  atan2Test(atan2BF16ReciprocalTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within y/x > 2^24", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 128),_,_), 2)
  atan2Test(atan2BF16ReciprocalTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within y/x > 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_), 2)
  atan2Test(atan2BF16ReciprocalTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within 1 < y/x < 2^12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_), 2)
  atan2Test(atan2BF16ReciprocalTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), 1.0,_,_), 2)
  atan2Test(atan2BF16ReciprocalTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within y/x < 2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_), 2)
}
