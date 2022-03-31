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

class MathFuncATan2Stage2SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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

  def atan2Test(t_rec : FuncTableInt, ts : Seq[FuncTableInt], spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generatorX      : ( (RealSpec, Random) => RealGeneric),
    generatorYoverX : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int,
    toleranceAtStage1 : Int = 1
    ) = {
    test(s"atan2(x), format ${spec.toStringShort}, ${generatorStr}") {
      var maxError   = 0.0
      var xatMaxError = 0.0
      var yatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for(i <- 1 to n) {
        val y_over_x = generatorYoverX(spec,r)
        val ysgn = if(r.nextBoolean()) {1.0} else {-1.0}
        val x  = generatorX(spec,r)
        val x0 = x.toDouble
        val y  = new RealGeneric(spec, x0 * y_over_x.toDouble * ysgn)
        val y0 = y.toDouble

        val minxy = if(abs(x0) < abs(y0)) {x0} else {y0}
        val maxxy = if(abs(x0) < abs(y0)) {y0} else {x0}
        val stage1z = new RealGeneric(spec, abs(minxy / maxxy))

        val z0   = math.atan2(y0, x0)
        val z0r  = new RealGeneric(spec, z0)
        val s1   = ATan2Stage1Sim.atan2Stage1SimGeneric( t_rec, y, x )
        val erriS1 = errorLSB(s1._1, stage1z.toDouble)
        if(abs(erriS1.toInt) > toleranceAtStage1) {
          println(f"WARN: ${toleranceAtStage1} < error in Stage1: sim(${s1._1.toDouble}) != ref(${stage1z.toDouble})")
        }

        val zi   = ATan2Stage2Sim.atan2Stage2SimGeneric( ts, s1._1, s1._2, s1._3, s1._4 )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toInt

        if (x0.isInfinity) {
          assert(zi.isNaN)
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if (erri.abs>tolerance) {
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

            println(f"test: x   = ${x0}")
            println(f"test: y   = ${y0}")
            println(f"test: y/x = ${min(abs(x0), abs(y0))/max(abs(x0),abs(y0))}")
            println(f"test: ref = ${z0}")
            println(f"test: sim = ${zd}")
            println(f"test: s1  = (${s1._1.sgn}|${s1._1.ex}|${s1._1.man.toLong.toBinaryString}), status = ${s1._2}, special = ${s1._3}, ysgn = ${s1._4}")
            println(f"test: test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman.toLong.toBinaryString}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
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
          assert(erri.abs<=tolerance)

          if (maxError < erri.abs) {
            maxError = erri.abs
            xatMaxError = x0
            yatMaxError = y0
            zatMaxError = zd
          }
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      println(f"${generatorStr} Summary")
      if(maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
              + f"${zatMaxError} != ${sin(xatMaxError)}, "
              + f"diff = ${zatMaxError - sin(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- $k errors positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val nOrderFP32 = 2
  val adrWFP32 = 8
  val extraBitsFP32 = 3
  val atan2FP32ReciprocalTableI = ReciprocalSim.reciprocalTableGeneration(
        nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)
  val atan2FP32ATanTableI       = ATan2Stage2Sim.atanTableGeneration(
        nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)

  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r,
    "Test Within y/x > 2^24", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 128),_,_), 1)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r,
    "Test Within y/x > 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_), 1)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r,
    "Test Within 1 < y/x < 2^12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_), 2)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r,
    "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), 1.0,_,_), 2)  // XXX
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r,
    "Test Within y/x < 2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_), 1)


  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 1
  val atan2BF16ReciprocalTableI = ReciprocalSim.reciprocalTableGeneration(
        nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)
  val atan2BF16ATanTableI       = ATan2Stage2Sim.atanTableGeneration(
        nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)

  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within y/x > 2^24", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 128),_,_), 1)
  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within y/x > 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_), 1)
  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within 1 < y/x < 2^12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_), 2)
  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), 1.0,_,_), 2)  // XXX
  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r,
    "Test Within y/x < 2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_), 1)

  val float48Spec = new RealSpec(10, 511, 37)

  val nOrderFP48 = 3
  val adrWFP48 = 10
  val extraBitsFP48 = 4
  val atan2FP48ReciprocalTableI = ReciprocalSim.reciprocalTableGeneration(
        nOrderFP48, adrWFP48, float48Spec.manW, float48Spec.manW+extraBitsFP48)
  val atan2FP48ATanTableI       = ATan2Stage2Sim.atanTableGeneration(
        nOrderFP48, adrWFP48, float48Spec.manW, float48Spec.manW+extraBitsFP48)

  atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r,
    "Test Within y/x > 2^24", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 128),_,_), 1)
  atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r,
    "Test Within y/x > 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_), 1)
  atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r,
    "Test Within 1 < y/x < 2^12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_), 2)
  atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r,
    "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), 1.0,_,_), 2)  // XXX
  atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r,
    "Test Within y/x < 2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_), 1)

  val nOrderFP64 = 3
  val adrWFP64 = 12
  val extraBitsFP64 = 4
  val atan2FP64ReciprocalTableI = ReciprocalSim.reciprocalTableGeneration(
        nOrderFP64, adrWFP64, RealSpec.Float64Spec.manW, RealSpec.Float64Spec.manW+extraBitsFP64)
  val atan2FP64ATanTableI       = ATan2Stage2Sim.atanTableGeneration(
        nOrderFP64, adrWFP64, RealSpec.Float64Spec.manW, RealSpec.Float64Spec.manW+extraBitsFP64)

  atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r,
    "Test Within y/x > 2^24", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 128),_,_), 1, 3)
  atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r,
    "Test Within y/x > 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_), 1, 3)
  atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r,
    "Test Within 1 < y/x < 2^12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_), 3, 3)
  atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r,
    "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), 1.0,_,_), 7, 3)  // XXX
  atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r,
    "Test Within y/x < 2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_), 7, 3)



}
