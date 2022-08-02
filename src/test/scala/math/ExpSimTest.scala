//package rial.tests


import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}



//import scopt.OptionParser

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math._
import rial.util.ScalaUtil._
import rial.testUtil.ScalaTestUtil._ // errorLSB
import rial.arith._
import rial.table._

class ExpSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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

  var counter = 0
  val specialValues = Seq(
      0.0,
      0.0 + (1.0                    ) * pow(2.0, -1022), // for Double
      0.0 + (1.0 + 1 * pow(2.0, -52)) * pow(2.0, -1022), // for Double
      0.0 + (1.0 + 2 * pow(2.0, -52)) * pow(2.0, -1022), // for Double
      0.0 + (1.0 + 3 * pow(2.0, -52)) * pow(2.0, -1022), // for Double
      0.0 + (1.1                    ) * pow(2.0, -1022), // for Double
      0.5,
      1.0,
      2.0,
      log(2.0),
      log(3.0),
      log(4.0),
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }

  def expTest(t : FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"exp(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

//       val libc = Native.loadLibrary("c", classOf[libc]).asInstanceOf[libc]

      var maxError    = SafeLong(0)
      var xatMaxError = 0.0
      var zatMaxError = 0.0
      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val z0   = exp(x0) // .toFloat
        val z0r  = new RealGeneric(spec, z0)

        val zi   = ExpSim.expSimGeneric(/*isPow2*/ false, t, x )
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r).toLong

        if (x0.isInfinity) {
          if(0 < x0) {
            assert(zi.isInfinite)
          } else {
            assert(zi.isZero)
          }
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if (erri.abs > tolerance) {
            val xsgn = bit(spec.W-1, x.value).toInt
            val xexp = slice(spec.manW, spec.exW, x.value)
            val xman = x.value & maskSL(spec.manW)

            val zsimsgn = bit(spec.W-1, zi.value).toInt
            val zsimexp = slice(spec.manW, spec.exW, zi.value)
            val zsimman = zi.value & maskSL(spec.manW)

            val zrefsgn = bit(spec.W-1, z0r.value).toInt
            val zrefexp = slice(spec.manW, spec.exW, z0r.value)
            val zrefman = z0r.value & maskSL(spec.manW)

            println(f"z0  = ${z0}")
            println(f"z0r = ${z0r.sgn}| ${z0r.ex}| ${z0r.man.toLong.toBinaryString}")
            println(f"zi  = ${zi.sgn }| ${zi .ex}| ${zi .man.toLong.toBinaryString}")
            println(f"zd  = ${zd}")
            println(f"test: x   = ${x0}(${x.sgn}|${x.ex}(${x.ex-x.spec.exBias})|${x.man.toLong.toBinaryString})")
            println(f"test: ref = ${z0}(${zrefsgn}|${zrefexp}(${zrefexp-x.spec.exBias})|${zrefman.toLong.toBinaryString})")
            println(f"test: sim = ${zd}(${zsimsgn}|${zsimexp}(${zsimexp-x.spec.exBias})|${zsimman.toLong.toBinaryString})")
            println(f"test: test(${zsimsgn}|${zsimexp}(${zsimexp - spec.exBias})|${zsimman.toLong.toBinaryString}(${zsimman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
          }

          if(erri != 0) {
            val errkey = erri.abs
            if( ! errs.contains(errkey)) {
              errs(errkey) = (0, 0)
            }
            if (erri >= 0) {
              errs(errkey) = (errs(errkey)._1 + 1, errs(errkey)._2)
            } else {
              errs(errkey) = (errs(errkey)._1, errs(errkey)._2 + 1)
            }
          }

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
              + f"${zatMaxError} != ${exp(xatMaxError)}, "
              + f"diff = ${zatMaxError - exp(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- ${k}%4d errors (${log2DownL(k)+1}%2d ULPs) positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")

    }
  }
// 
//   val nOrderFP32    = 2
//   val adrWFP32      = 8
//   val extraBitsFP32 = 3
// 
//   val pow2FP32TableI = ExpSim.pow2TableGeneration(
//     nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)
// 
//   expTest(pow2FP32TableI, RealSpec.Float32Spec, n, r,
//     "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 2)
//   expTest(pow2FP32TableI, RealSpec.Float32Spec, n, r,
//     "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 2)
// 
//   expTest(pow2FP32TableI, RealSpec.Float32Spec, n, r,
//     "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), 2)
//   expTest(pow2FP32TableI, RealSpec.Float32Spec, n, r,
//     "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), 2)
// 
//   expTest(pow2FP32TableI, RealSpec.Float32Spec, n, r,
//     "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 2)
//   expTest(pow2FP32TableI, RealSpec.Float32Spec, n, r,
//     "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 2)
// 
//   expTest(pow2FP32TableI, RealSpec.Float32Spec, n, r,
//     "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 2)
//   expTest(pow2FP32TableI, RealSpec.Float32Spec, n, r,
//     "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 2)
// 
//   expTest(pow2FP32TableI, RealSpec.Float32Spec, n, r,
//     "Test Special Values", generateSpecialValues(_,_), 1)
// 
//   val nOrderBF16    = 0
//   val adrWBF16      = 7
//   val extraBitsBF16 = 1
// 
//   val pow2BF16TableI = ExpSim.pow2TableGeneration(
//     nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)
// 
//   expTest(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 1)
//   expTest(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 1)
// 
//   expTest(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), 1)
//   expTest(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), 1)
// 
//   expTest(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 1)
//   expTest(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 1)
// 
//   expTest(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 1)
//   expTest(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 1)
// 
//   expTest(pow2BF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Special Values", generateSpecialValues(_,_), 2)

  val nOrderFP64    = 3
  val adrWFP64      = 12
  val extraBitsFP64 = 4

  val pow2FP64TableI = ExpSim.pow2TableGeneration(
    nOrderFP64, adrWFP64, RealSpec.Float64Spec.manW, RealSpec.Float64Spec.manW+extraBitsFP64)

//   expTest(pow2FP64TableI, RealSpec.Float64Spec, n, r,
//     "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 7)
//   expTest(pow2FP64TableI, RealSpec.Float64Spec, n, r,
//     "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 7)
// 
//   expTest(pow2FP64TableI, RealSpec.Float64Spec, n, r,
//     "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), 7)
//   expTest(pow2FP64TableI, RealSpec.Float64Spec, n, r,
//     "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), 7)
// 
//   expTest(pow2FP64TableI, RealSpec.Float64Spec, n, r,
//     "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 7)
//   expTest(pow2FP64TableI, RealSpec.Float64Spec, n, r,
//     "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 7)
// 
//   expTest(pow2FP64TableI, RealSpec.Float64Spec, n, r,
//     "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 7)
//   expTest(pow2FP64TableI, RealSpec.Float64Spec, n, r,
//     "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 7)

  expTest(pow2FP64TableI, RealSpec.Float64Spec, n, r,
    "Test Special Values", generateSpecialValues(_,_), 2)

}
