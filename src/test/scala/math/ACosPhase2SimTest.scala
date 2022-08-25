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

class ACosPhase2SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 1000000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000000").toInt
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
      -0.0,
       1.0,
      -1.0,
       1.1, // should be acos(1.0)
      -1.1, // should be acos(-1.0)
      Double.PositiveInfinity,
      Double.NegativeInfinity,
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def acosTest(tSqrt : FuncTableInt, tACos : FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int,
    toleranceAtPhase1 : Int = 1
    ) = {
    test(s"acos(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0;
      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val xref = if (x.toDouble > 1.0 && !x.isNaN && !x.isInfinite) {
          1.0
        } else if (x.toDouble < -1.0 && !x.isNaN && !x.isInfinite) {
          -1.0
        } else {
          x.toDouble
        }
        val stage1z = new RealGeneric(spec, sqrt(1.0 - abs(xref)))

        val z0   = math.acos(xref)
        val z0r  = new RealGeneric(spec, z0)

        val s1   = ACosPhase1Sim.acosPhase1SimGeneric( tSqrt, x )
        val erriS1 = errorLSB(s1._1, stage1z.toDouble)

        if(abs(erriS1.toInt) > toleranceAtPhase1) {
          println(f"WARN: ${toleranceAtPhase1} < error in Phase1: sim(${s1._1.toDouble}) != ref(${stage1z.toDouble})")
        }
        assert(abs(erriS1.toInt) <= toleranceAtPhase1)

        val zi   = ACosPhase2Sim.acosPhase2SimGeneric( tACos, s1._1, s1._2, s1._3 )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toInt

        if (z0r.isNaN) {
          assert(zi.isNaN, f"x = ${x0}, zref = NaN, zsim = not nan")
        } else {
          if (erri.abs>tolerance) {
            println(f"Error more than 2 LSB : x = ${x.toDouble}%14.7e, : refz = $z0%14.7e simz = ${zi.toDouble}%14.7e $errf%14.7e $erri%f")

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
            println(f"test: ref = ${z0}")
            println(f"test: sim = ${zd}")
            println(f"test: s1  = (${s1._1.sgn}|${s1._1.ex}|${s1._1.man.toLong.toBinaryString}), xneg = ${s1._2}, special = ${s1._3}")
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
            zatMaxError = zd
          }
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      println(f"${generatorStr} Summary")
      if(maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
              + f"${zatMaxError} != ${acos(xatMaxError)}, "
              + f"diff = ${zatMaxError - acos(xatMaxError)}, x = ${xatMaxError}")
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
  val manWFP32 = RealSpec.Float32Spec.manW
  val fracWFP32 = RealSpec.Float32Spec.manW + extraBitsFP32
  val acosFP32SqrtTableI = ACosPhase1Sim.sqrtTableGeneration(
        nOrderFP32, adrWFP32, manWFP32, fracWFP32)
  val acosFP32ACosTableI = ACosPhase2Sim.acosTableGeneration(
        nOrderFP32, adrWFP32, manWFP32, fracWFP32)

  acosTest(acosFP32SqrtTableI, acosFP32ACosTableI, RealSpec.Float32Spec, n, r, "Test Within 0      < x <= 2^-4",   generateRealWithin(0.0,                pow(2.0, -4)      ,_,_), 3, 1)
  acosTest(acosFP32SqrtTableI, acosFP32ACosTableI, RealSpec.Float32Spec, n, r, "Test Within 2^-4   < x <= 1-2^-4", generateRealWithin(pow(2.0, -4),       1.0 - pow(2.0, -4),_,_), 3, 1)
  acosTest(acosFP32SqrtTableI, acosFP32ACosTableI, RealSpec.Float32Spec, n, r, "Test Within 1-2^-4 < x <= 1",      generateRealWithin(1.0 - pow(2.0, -4), 1.0               ,_,_), 3, 1)

  acosTest(acosFP32SqrtTableI, acosFP32ACosTableI, RealSpec.Float32Spec, n, r, "Test Within  0      >= x > -2^-4",   generateRealWithin(-pow(2.0,-4),     0.0,             _,_), 3, 1)
  acosTest(acosFP32SqrtTableI, acosFP32ACosTableI, RealSpec.Float32Spec, n, r, "Test Within -2^-4   >= x > -1-2^-4", generateRealWithin(-1.0+pow(2.0,-4), -pow(2.0,-4),    _,_), 3, 1)
  acosTest(acosFP32SqrtTableI, acosFP32ACosTableI, RealSpec.Float32Spec, n, r, "Test Within -1-2^-4 >= x > -1",      generateRealWithin(-1.0,             -1.0+pow(2.0,-4),_,_), 3, 1)

  acosTest(acosFP32SqrtTableI, acosFP32ACosTableI, RealSpec.Float32Spec, n, r, "Test Special Values", generateSpecialValues(_,_), 1)



  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 1
  val acosBF16SqrtTableI = ACosPhase1Sim.sqrtTableGeneration(
        nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)
  val acosBF16ACosTableI       = ACosPhase2Sim.acosTableGeneration(
        nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)

  acosTest(acosBF16SqrtTableI, acosBF16ACosTableI, RealSpec.BFloat16Spec, n, r, "Test Within 0      < x <= 2^-4",   generateRealWithin(0.0,                pow(2.0, -4)      ,_,_), 3, 1)
  acosTest(acosBF16SqrtTableI, acosBF16ACosTableI, RealSpec.BFloat16Spec, n, r, "Test Within 2^-4   < x <= 1-2^-4", generateRealWithin(pow(2.0, -4),       1.0 - pow(2.0, -4),_,_), 3, 1)
  acosTest(acosBF16SqrtTableI, acosBF16ACosTableI, RealSpec.BFloat16Spec, n, r, "Test Within 1-2^-4 < x <= 1",      generateRealWithin(1.0 - pow(2.0, -4), 1.0               ,_,_), 3, 1)

  acosTest(acosBF16SqrtTableI, acosBF16ACosTableI, RealSpec.BFloat16Spec, n, r, "Test Within  0      >= x > -2^-4",   generateRealWithin(-pow(2.0,-4),     0.0,             _,_), 3, 1)
  acosTest(acosBF16SqrtTableI, acosBF16ACosTableI, RealSpec.BFloat16Spec, n, r, "Test Within -2^-4   >= x > -1-2^-4", generateRealWithin(-1.0+pow(2.0,-4), -pow(2.0,-4),    _,_), 3, 1)
  acosTest(acosBF16SqrtTableI, acosBF16ACosTableI, RealSpec.BFloat16Spec, n, r, "Test Within -1-2^-4 >= x > -1",      generateRealWithin(-1.0,             -1.0+pow(2.0,-4),_,_), 3, 1)

  acosTest(acosBF16SqrtTableI, acosBF16ACosTableI, RealSpec.BFloat16Spec, n, r, "Test Special Values", generateSpecialValues(_,_), 1)


//   val float48Spec = new RealSpec(10, 511, 37)
//
//   val nOrderFP48 = 3
//   val adrWFP48 = 10
//   val extraBitsFP48 = 4
//   val acosFP48SqrtTableI = ACosPhase1Sim.sqrtTableGeneration(
//         nOrderFP48, adrWFP48, float48Spec.manW, float48Spec.manW+extraBitsFP48)
//   val acosFP48ACosTableI       = ACosPhase2Sim.acosTableGeneration(
//         nOrderFP48, adrWFP48, float48Spec.manW, float48Spec.manW+extraBitsFP48)
//
//   acosTest(acosFP48SqrtTableI, acosFP48ACosTableI, float48Spec, n, r, "Test Within 0      < x <= 2^-4",   generateRealWithin(0.0,                pow(2.0, -4)      ,_,_), 3, 1)
//   acosTest(acosFP48SqrtTableI, acosFP48ACosTableI, float48Spec, n, r, "Test Within 2^-4   < x <= 1-2^-4", generateRealWithin(pow(2.0, -4),       1.0 - pow(2.0, -4),_,_), 3, 1)
//   acosTest(acosFP48SqrtTableI, acosFP48ACosTableI, float48Spec, n, r, "Test Within 1-2^-4 < x <= 1",      generateRealWithin(1.0 - pow(2.0, -4), 1.0               ,_,_), 3, 1) // XXX
//
//   acosTest(acosFP48SqrtTableI, acosFP48ACosTableI, float48Spec, n, r, "Test Within  0      >= x > -2^-4",   generateRealWithin(-pow(2.0,-4),     0.0,             _,_), 3, 1)
//   acosTest(acosFP48SqrtTableI, acosFP48ACosTableI, float48Spec, n, r, "Test Within -2^-4   >= x > -1-2^-4", generateRealWithin(-1.0+pow(2.0,-4), -pow(2.0,-4),    _,_), 3, 1)
//   acosTest(acosFP48SqrtTableI, acosFP48ACosTableI, float48Spec, n, r, "Test Within -1-2^-4 >= x > -1",      generateRealWithin(-1.0,             -1.0+pow(2.0,-4),_,_), 3, 1)
//
//   acosTest(acosFP48SqrtTableI, acosFP48ACosTableI, float48Spec, n, r, "Test Special Values", generateSpecialValues(_,_), 1)
//
//
//
//   val nOrderFP64 = 3
//   val adrWFP64 = 12
//   val extraBitsFP64 = 4
//   val acosFP64SqrtTableI = ACosPhase1Sim.sqrtTableGeneration(
//         nOrderFP64, adrWFP64, RealSpec.Float64Spec.manW, RealSpec.Float64Spec.manW+extraBitsFP64)
//   val acosFP64ACosTableI       = ACosPhase2Sim.acosTableGeneration(
//         nOrderFP64, adrWFP64, RealSpec.Float64Spec.manW, RealSpec.Float64Spec.manW+extraBitsFP64)
//
//   acosTest(acosFP64SqrtTableI, acosFP64ACosTableI, RealSpec.Float64Spec, n, r, "Test Within 0      < x <= 2^-4",   generateRealWithin(0.0,                pow(2.0, -4)      ,_,_), 15, 7)
//   acosTest(acosFP64SqrtTableI, acosFP64ACosTableI, RealSpec.Float64Spec, n, r, "Test Within 2^-4   < x <= 1-2^-4", generateRealWithin(pow(2.0, -4),       1.0 - pow(2.0, -4),_,_), 15, 7)
//   acosTest(acosFP64SqrtTableI, acosFP64ACosTableI, RealSpec.Float64Spec, n, r, "Test Within 1-2^-4 < x <= 1",      generateRealWithin(1.0 - pow(2.0, -4), 1.0               ,_,_), 15, 7)
//
//   acosTest(acosFP64SqrtTableI, acosFP64ACosTableI, RealSpec.Float64Spec, n, r, "Test Within  0      >= x > -2^-4",   generateRealWithin(-pow(2.0,-4),     0.0,             _,_), 15, 7)
//   acosTest(acosFP64SqrtTableI, acosFP64ACosTableI, RealSpec.Float64Spec, n, r, "Test Within -2^-4   >= x > -1-2^-4", generateRealWithin(-1.0+pow(2.0,-4), -pow(2.0,-4),    _,_), 15, 7)
//   acosTest(acosFP64SqrtTableI, acosFP64ACosTableI, RealSpec.Float64Spec, n, r, "Test Within -1-2^-4 >= x > -1",      generateRealWithin(-1.0,             -1.0+pow(2.0,-4),_,_), 15, 7)
//
//   acosTest(acosFP64SqrtTableI, acosFP64ACosTableI, RealSpec.Float64Spec, n, r, "Test Special Values", generateSpecialValues(_,_), 1)

}
