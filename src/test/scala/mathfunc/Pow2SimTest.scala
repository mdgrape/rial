//package rial.tests


import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}



//import scopt.OptionParser

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.mathfunc._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

// import com.sun.jna._
// trait libc extends Library {
//   def powf(x: Float, y:Float):Float
// }

class MathFuncPow2SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 1000000

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

  def pow2Test(t : FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"pow2(x), format ${spec.toStringShort}, ${generatorStr}") {

//       val libc = Native.loadLibrary("c", classOf[libc]).asInstanceOf[libc]

      var maxError    = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      var err1lsbPos = 0
      var err1lsbNeg = 0
      var err2lsbPos = 0
      var err2lsbNeg = 0
      var errNlsbPos = 0
      var errNlsbNeg = 0
      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

//         val z0   = libc.powf(2.0f, x0)
//         val z0b  = java.lang.Float.floatToIntBits(z0)
//         val z0sgn = bit(31, z0b)
//         val z0ex  = slice(spec.manW, spec.exW,  z0b)
//         val z0man = slice(0,         spec.manW, z0b)
//         val z0r  = new RealGeneric(spec, z0sgn, z0ex, z0man)

        val z0   = pow(2.0, x0)
        val z0r  = new RealGeneric(spec, z0)

        val zi   = MathFuncExpSim.expSimGeneric( /*isPow2*/true, t, x )
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r.toDouble)

        if (x0.isInfinity) {
          if(0 < x0) {
            assert(zi.isInfinite)
          } else {
            assert(zi.isZero)
          }
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if (erri.abs>=2.0) {
            val xsgn = bit(spec.W-1, x.value).toInt
            val xexp = slice(spec.manW, spec.exW, x.value)
            val xman = x.value & maskSL(spec.manW)

            val zsimsgn = bit(spec.W-1, zi.value).toInt
            val zsimexp = slice(spec.manW, spec.exW, zi.value)
            val zsimman = zi.value & maskSL(spec.manW)

            val zrefsgn = bit(spec.W-1, z0r.value).toInt
            val zrefexp = slice(spec.manW, spec.exW, z0r.value)
            val zrefman = z0r.value & maskSL(spec.manW)

            if(erri.abs>2.0) {
              println(f"test: x   = ${x0}(${x.sgn}|${x.ex}(${x.ex-x.spec.exBias})|${x.man.toLong.toBinaryString})")
              println(f"test: ref = ${z0}(${zrefsgn}|${zrefexp}(${zrefexp-x.spec.exBias})|${zrefman.toLong.toBinaryString})")
              println(f"test: sim = ${zd}(${zsimsgn}|${zsimexp}(${zsimexp-x.spec.exBias})|${zsimman.toLong.toBinaryString})")
              println(f"test: test(${zsimsgn}|${zsimexp}(${zsimexp - spec.exBias})|${zsimman.toLong.toBinaryString}(${zsimman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
            }

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
          assert(erri.abs<=tolerance)

          if (maxError < erri.abs) {
            maxError    = erri.abs
            xatMaxError = x0
            zatMaxError = zd
          }
        }
      }
      println(f"${generatorStr} Summary")
      if(maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
              + f"${zatMaxError} != ${pow(2.0, xatMaxError)}, "
              + f"diff = ${zatMaxError - pow(2.0, xatMaxError)}, x = ${xatMaxError}")
      }
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
      println(f"N=$n%d : 2LSB errors positive $err2lsbPos%d / negative $err2lsbNeg%d")
      println(f"N=$n%d : 2<   errors positive $errNlsbPos%d / negative $errNlsbNeg%d")
      println( "---------------------------------------------------------------")
    }
  }

  val nOrderFP32 = 2
  val adrWFP32 = 8
  val extraBitsFP32 = 3

  val pow2FP32TableI = MathFuncExpSim.pow2TableGeneration(
    nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)

  pow2Test(pow2FP32TableI, RealSpec.Float32Spec, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 2)
  pow2Test(pow2FP32TableI, RealSpec.Float32Spec, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 2)

  pow2Test(pow2FP32TableI, RealSpec.Float32Spec, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), 2)
  pow2Test(pow2FP32TableI, RealSpec.Float32Spec, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), 2)

  pow2Test(pow2FP32TableI, RealSpec.Float32Spec, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 2)
  pow2Test(pow2FP32TableI, RealSpec.Float32Spec, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 2)

  pow2Test(pow2FP32TableI, RealSpec.Float32Spec, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 2)
  pow2Test(pow2FP32TableI, RealSpec.Float32Spec, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 2)

  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 1

  val pow2BF16TableI = MathFuncExpSim.pow2TableGeneration(
    nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)

  pow2Test(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 1)
  pow2Test(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 1)

  pow2Test(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), 1)
  pow2Test(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), 1)

  pow2Test(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 1)
  pow2Test(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 1)

  pow2Test(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 1)
  pow2Test(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 1)

}
