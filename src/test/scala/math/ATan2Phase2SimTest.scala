//package rial.tests

import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

class ATan2Phase2SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 1000000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRealWithin( from : Double, to : Double, spec: RealSpec, idx: Int, r : Random ) = {
    val rD : Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  val specialValues = Seq(
      ( 1.0,  1.0),
      ( 1.0, -1.0),
      (-1.0,  1.0),
      (-1.0, -1.0),

      ( 1.99999988,  1.99999988),
      ( 1.99999988, -1.99999988),
      (-1.99999988,  1.99999988),
      (-1.99999988, -1.99999988),

      ( 2.0,  2.0),
      ( 2.0, -2.0),
      (-2.0,  2.0),
      (-2.0, -2.0),

      ( 1.41421356,  1.0),
      ( 1.41421356, -1.0),
      (-1.41421356,  1.0),
      (-1.41421356, -1.0),

      ( 1.0,  1.41421356),
      ( 1.0, -1.41421356),
      (-1.0,  1.41421356),
      (-1.0, -1.41421356),

      ( 1.41421356,  1.41421356),
      ( 1.41421356, -1.41421356),
      (-1.41421356,  1.41421356),
      (-1.41421356, -1.41421356),

      ( 3.14159265,  3.14159265),
      ( 3.14159265, -3.14159265),
      (-3.14159265,  3.14159265),
      (-3.14159265, -3.14159265),

      ( 100.0,  100.0),
      ( 100.0, -100.0),
      (-100.0,  100.0),
      (-100.0, -100.0),

      (Double.PositiveInfinity,  3.1415926),
      (Double.PositiveInfinity, -3.1415926),
      (Double.NegativeInfinity,  3.1415926),
      (Double.NegativeInfinity, -3.1415926),

      ( 3.1415926, Double.PositiveInfinity),
      (-3.1415926, Double.PositiveInfinity),
      ( 3.1415926, Double.NegativeInfinity),
      (-3.1415926, Double.NegativeInfinity),

      (Double.PositiveInfinity, Double.PositiveInfinity),
      (Double.PositiveInfinity, Double.NegativeInfinity),
      (Double.NegativeInfinity, Double.PositiveInfinity),
      (Double.NegativeInfinity, Double.NegativeInfinity),
    )
  def generateSpecialValuesX( spec: RealSpec, idx: Int, r: Random ) = {

    val i = idx % specialValues.length
    val x = specialValues(i)._1

    val eps = math.pow(2.0, -spec.manW)

    val rnd = r.nextInt(5) - 2 // [-2, 2]
    new RealGeneric(spec, x + (rnd * eps))
  }
  def generateSpecialValuesY( spec: RealSpec, idx: Int, r: Random ) = {

    val i = idx % specialValues.length
    val y = specialValues(i)._2

    val eps = math.pow(2.0, -spec.manW)

    val rnd = r.nextInt(5) - 2 // [-2, 2]
    new RealGeneric(spec, y + (rnd * eps))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def atan2Test(t_rec : FuncTableInt, t : FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generatorX      : ( (RealSpec, Int, Random) => RealGeneric),
    generatorY      : ( (RealSpec, Int, Random) => RealGeneric),
    tolerance       : Int,
    toleranceAtPhase1 : Int = 3
    ) = {
    test(s"atan2(x), format ${spec.toStringShort}, ${generatorStr}") {
      var maxError   = 0.0
      var xatMaxError = 0.0
      var yatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for(i <- 1 to n) {
        val x  = generatorX(spec,i,r)
        val y  = generatorY(spec,i,r)
        val x0 = x.toDouble
        val y0 = y.toDouble

        val minxy = if(abs(x0) < abs(y0)) {x0} else {y0}
        val maxxy = if(abs(x0) < abs(y0)) {y0} else {x0}
        val stage1z = new RealGeneric(spec, abs(minxy / maxxy))

        val z0   = math.atan2(y0, x0)
        val z0r  = new RealGeneric(spec, z0)
        val s1   = ATan2Phase1Sim.atan2Phase1SimGeneric( t_rec, y, x )
        val erriS1 = errorLSB(s1, stage1z.toDouble)
        // if(abs(erriS1.toInt) > toleranceAtPhase1) {
        //   println(f"WARN: ${toleranceAtPhase1} < error in Phase1: sim(${s1.toDouble}) != ref(${stage1z.toDouble})")
        // }

        val zi   = ATan2Phase2Sim.atan2Phase2SimGeneric( t, s1 )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toInt

        if (z0r.isInfinite && z0r.toDouble > 0) {
          assert(zi.isInfinite, f"x = ${x0}, y = ${y0}, zref = Inf, zsim = not inf")
        } else if (z0r.isInfinite && z0r.toDouble < 0) {
          assert(zi.isInfinite, f"x = ${x0}, y = ${y0}, zref = -Inf, zsim = not inf")
        } else if (z0r.isNaN) {
          assert(zi.isNaN, f"x = ${x0}, y = ${y0}, zref = NaN, zsim = not nan")
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
            println(f"test: s1  = (${s1.sgn}|${s1.ex}|${s1.man.toLong.toBinaryString})")
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
              + f"${zatMaxError} != ${atan2(yatMaxError, xatMaxError)}, "
              + f"diff = ${zatMaxError - atan2(yatMaxError, xatMaxError)}, x=${xatMaxError}, y=${yatMaxError}")
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
  val atan2FP32ATanTableI       = ATan2Phase2Sim.atanTableGeneration(
        nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)

  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within  2^24  < y/x <  inf",               generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within -2^24  > y/x > -inf",               generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within  2^24  < y/x <  inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within -2^24  > y/x > -inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_,_), 4)

  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within  2^12  < y/x <  2^24",              generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within -2^12  > y/x > -2^24",              generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within  2^12  < y/x <  2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within -2^12  > y/x > -2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_,_), 4)

  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within  1     < y/x <  2^12",              generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within -1     > y/x > -2^12",              generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within  1     < y/x <  2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within -1     > y/x > -2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_,_), 4)

  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within  2^-12 < y/x <  1",                 generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within -2^-12 > y/x > -1",                 generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within  2^-12 < y/x <  1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within -2^-12 > y/x > -1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_,_), 4)

  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within 0     < y/x <  2^-12",              generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin(           0.0, pow(2.0, -12),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within 0     > y/x > -2^-12",              generateRealWithin(-1.0,           1.0,          _,_,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin(           0.0, pow(2.0, -12),_,_,_), 4)
  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_,_), 4)

  atan2Test(atan2FP32ReciprocalTableI, atan2FP32ATanTableI, RealSpec.Float32Spec, n, r, "Test Special Values", generateSpecialValuesX(_,_,_), generateSpecialValuesY(_,_,_), 4)

  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 1
  val atan2BF16ReciprocalTableI = ReciprocalSim.reciprocalTableGeneration(
        nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)
  val atan2BF16ATanTableI       = ATan2Phase2Sim.atanTableGeneration(
        nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)

  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r, "Test Within y/x > 2^24",      generateRealWithin(-1.0, 1.0,_,_,_), generateRealWithin(pow(2.0, 24),  pow(2.0, 128),_,_,_), 3)
  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r, "Test Within y/x > 2^12",      generateRealWithin(-1.0, 1.0,_,_,_), generateRealWithin(pow(2.0, 12),  pow(2.0, 24), _,_,_), 3)
  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r, "Test Within 1 < y/x < 2^12",  generateRealWithin(-1.0, 1.0,_,_,_), generateRealWithin(1.0,           pow(2.0, 8),  _,_,_), 3)
  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r, "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_,_), generateRealWithin(pow(2.0, -12), 1.0,          _,_,_), 7) // XXX
  atan2Test(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, RealSpec.BFloat16Spec, n, r, "Test Within y/x < 2^-12",     generateRealWithin(-1.0, 1.0,_,_,_), generateRealWithin(0.0,           pow(2.0, -12),_,_,_), 7) // XXX

  // val float48Spec = new RealSpec(10, 511, 37)
  //
  // val nOrderFP48 = 3
  // val adrWFP48 = 10
  // val extraBitsFP48 = 4
  // val atan2FP48ReciprocalTableI = ReciprocalSim.reciprocalTableGeneration(
  //       nOrderFP48, adrWFP48, float48Spec.manW, float48Spec.manW+extraBitsFP48)
  // val atan2FP48ATanTableI       = ATan2Phase2Sim.atanTableGeneration(
  //       nOrderFP48, adrWFP48, float48Spec.manW, float48Spec.manW+extraBitsFP48)
  //
  // atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r, "Test Within y/x > 2^24",      generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 128),_,_), 4)
  // atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r, "Test Within y/x > 2^12",      generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_),  4)
  // atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r, "Test Within 1 < y/x < 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_),            4)
  // atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r, "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), 1.0,_,_),          4)  // XXX
  // atan2Test(atan2FP48ReciprocalTableI, atan2FP48ATanTableI, float48Spec, n, r, "Test Within y/x < 2^-12",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_),          4)
  //
  // val nOrderFP64 = 3
  // val adrWFP64 = 12
  // val extraBitsFP64 = 4
  // val atan2FP64ReciprocalTableI = ReciprocalSim.reciprocalTableGeneration(
  //       nOrderFP64, adrWFP64, RealSpec.Float64Spec.manW, RealSpec.Float64Spec.manW+extraBitsFP64)
  // val atan2FP64ATanTableI       = ATan2Phase2Sim.atanTableGeneration(
  //       nOrderFP64, adrWFP64, RealSpec.Float64Spec.manW, RealSpec.Float64Spec.manW+extraBitsFP64)
  //
  // atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r, "Test Within y/x > 2^24",      generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 128),_,_), 7, 3)
  // atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r, "Test Within y/x > 2^12",      generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_),  7, 3)
  // atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r, "Test Within 1 < y/x < 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_),            7, 3)
  // atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r, "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), 1.0,_,_),          7, 3)  // XXX
  // atan2Test(atan2FP64ReciprocalTableI, atan2FP64ATanTableI, RealSpec.Float64Spec, n, r, "Test Within y/x < 2^-12",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_),          7, 3)
}
