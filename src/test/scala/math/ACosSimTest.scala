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

class ACosSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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
      0.0,                -0.0,
      1.0,                -1.0,
      1.0 / 2.0,          -1.0 / 2.0,
      sqrt(2.0) / 2.0,    -sqrt(2.0) / 2.0,
      sqrt(3.0) / 2.0,    -sqrt(3.0) / 2.0,
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

  def acosTest(t : Seq[FuncTableInt], tSqrt: FuncTableInt, exTable: Option[Seq[FuncTableInt]],
    spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"acos(x), format ${spec.toStringShort}, ${generatorStr}") {

      var maxError    = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val z0   = acos(x0)
        val z0r  = new RealGeneric(spec, z0)

        val zi   = ACosSim.acosSimGeneric(t, tSqrt, x, exTable )
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toLong
//         val errf = zi.toDouble - z0r.toDouble

        if (x0.isInfinity) {
          if(0 < x0) {
            assert(zi.isInfinite)
          } else {
            assert(zi.isZero)
          }
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else if (z0r.isNaN) {
          assert(zi.isNaN)
        } else {
          if(erri.abs>tolerance) {
            val xsgn = bit(spec.W-1, x.value).toInt
            val xexp = slice(spec.manW, spec.exW, x.value)
            val xman = x.value & maskSL(spec.manW)

            val zsimsgn = bit(spec.W-1, zi.value).toInt
            val zsimexp = slice(spec.manW, spec.exW, zi.value)
            val zsimman = zi.value & maskSL(spec.manW)

            val zrefsgn = bit(spec.W-1, z0r.value).toInt
            val zrefexp = slice(spec.manW, spec.exW, z0r.value)
            val zrefman = z0r.value & maskSL(spec.manW)

            println(f"test: x   = ${x0}%20.15g(${x.sgn}|${x.ex}(${x.ex-x.spec.exBias})|${x.man.toLong.toBinaryString})")
            println(f"test: ref = ${z0}%20.15g(${zrefsgn}|${zrefexp}(${zrefexp-x.spec.exBias})|${zrefman.toLong.toBinaryString})")
            println(f"test: sim = ${zd}%20.15g(${zsimsgn}|${zsimexp}(${zsimexp-x.spec.exBias})|${zsimman.toLong.toBinaryString})")
            println(f"test: test(${zsimsgn}|${zsimexp}(${zsimexp - spec.exBias})|${zsimman.toLong.toBinaryString}(${zsimman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
          }
          assert(erri.abs<=tolerance)

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
            maxError    = erri.abs.toDouble
            xatMaxError = x0
            zatMaxError = zd
          }
        }
      }
      println(f"${generatorStr} Summary")
      if(maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
              + f"${zatMaxError} != ${acos(xatMaxError)}, "
              + f"diff = ${zatMaxError - acos(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- ${k}%4d errors (${log2DownL(k)+1}%2d ULPs) positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3

  val acosFP32Table = ACosSim.acosTableGeneration(RealSpec.Float32Spec,
    nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)
  val sqrtFP32Table = ACosSim.sqrtTableGeneration(
    nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)

  val smallValueFP32 = -4

  // acos(-|x|) is in [pi/2, pi]. this does not require super high resolution.
  acosTest(acosFP32Table, sqrtFP32Table, None, RealSpec.Float32Spec, n, r,
    "Test close to 0.0:  [0.0, -2^-4]", generateRealWithin(-pow(2.0, smallValueFP32), 0.0,_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, None, RealSpec.Float32Spec, n, r,
    "Test normal table x < 0.5:  [-2^-4, -0.5]", generateRealWithin(-0.5, -pow(2.0, smallValueFP32),_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, None, RealSpec.Float32Spec, n, r,
    "Test normal table x > 0.5:  [-0.5, -1+2^-4]", generateRealWithin(-1.0+pow(2.0, smallValueFP32)-pow(2.0, -23), -0.5-pow(2.0, -23),_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, None, RealSpec.Float32Spec, n, r,
    "Test Puiseux: close to 1.0: [-1-2^+4, -1.0]", generateRealWithin(-1.0, -1.0+pow(2.0, smallValueFP32)-pow(2.0, -23),_,_), 3) // 2ULPs

  acosTest(acosFP32Table, sqrtFP32Table, None, RealSpec.Float32Spec, n, r,
    "Test close to 0.0:  [0.0, 2^-4]", generateRealWithin(0.0, pow(2.0, smallValueFP32),_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, None, RealSpec.Float32Spec, n, r,
    "Test normal table x < 0.5:  [2^-4, 0.5]", generateRealWithin(pow(2.0, smallValueFP32), 0.5,_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, None, RealSpec.Float32Spec, n, r,
    "Test normal table x > 0.5:  [0.5, 1-2^-4]", generateRealWithin(0.5+pow(2.0, -23), 1.0-pow(2.0, smallValueFP32)+pow(2.0, -23),_,_), 3) // 2ULPs
  acosTest(acosFP32Table, sqrtFP32Table, None, RealSpec.Float32Spec, n, r,
    "Test Puiseux: close to 1.0: [1-2^-4, 1.0]", generateRealWithin(1.0-pow(2.0, smallValueFP32)+pow(2.0, -23), 1.0,_,_), 3) // 2ULPs

  acosTest(acosFP32Table, sqrtFP32Table, None, RealSpec.Float32Spec, n, r,
    "Test special value", generateSpecialValues(_,_), 1)

  val nOrderBF16    = 0
  val adrWBF16      = 7
  val extraBitsBF16 = 1

  val acosBF16Table = ACosSim.acosTableGeneration(RealSpec.BFloat16Spec,
    nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)
  val sqrtBF16Table = ACosSim.sqrtTableGeneration(
    nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)

  val acosExBF16Table = ACosSim.acosExTableGeneration(RealSpec.BFloat16Spec,
    nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)

  // since acos[BF16] does not use series expansion, we actually don't need to check the edge cases basically.

  val smallValueBF16 = -4

  acosTest(acosBF16Table, sqrtBF16Table, Some(acosExBF16Table), RealSpec.BFloat16Spec, n, r,
    f"Test close to -1.0: [-1, -1+2^${smallValueBF16}]", generateRealWithin(-1.0, -1.0+pow(2.0, smallValueBF16),_,_), 1)
  acosTest(acosBF16Table, sqrtBF16Table, Some(acosExBF16Table), RealSpec.BFloat16Spec, n, r,
    f"Test normal table -1.0 < x < -0.5:  [-1+2^${smallValueBF16}, -0.5]", generateRealWithin(-1.0, -0.5,_,_), 1)
  acosTest(acosBF16Table, sqrtBF16Table, Some(acosExBF16Table), RealSpec.BFloat16Spec, n, r,
    f"Test normal table -0.5 < x < 0:  [-0.5, 0.0]", generateRealWithin(-0.5, 0.0,_,_), 1)
  acosTest(acosBF16Table, sqrtBF16Table, Some(acosExBF16Table), RealSpec.BFloat16Spec, n, r,
    f"Test close to 0.0: [-2^${smallValueBF16},  0.0]", generateRealWithin(-pow(2.0, smallValueBF16), 0.0,_,_), 1)

  acosTest(acosBF16Table, sqrtBF16Table, Some(acosExBF16Table), RealSpec.BFloat16Spec, n, r,
    f"Test close to 0.0:  [0.0, 2^${smallValueBF16}]", generateRealWithin(0.0, pow(2.0, smallValueBF16),_,_), 1)
  acosTest(acosBF16Table, sqrtBF16Table, Some(acosExBF16Table), RealSpec.BFloat16Spec, n, r,
    f"Test normal table 0.0 < x < 0.5:  [0.0, 0.5]", generateRealWithin(0.0, 0.5,_,_), 1)
  acosTest(acosBF16Table, sqrtBF16Table, Some(acosExBF16Table), RealSpec.BFloat16Spec, n, r,
    f"Test normal table 0.5 < x < 1.0:  [0.5, 1.0]", generateRealWithin(0.5, 1.0,_,_), 1)
  acosTest(acosBF16Table, sqrtBF16Table, Some(acosExBF16Table), RealSpec.BFloat16Spec, n, r,
    f"Test close to 1.0: [1-2^${smallValueBF16}, 1.0]", generateRealWithin(1.0-pow(2.0, smallValueBF16)+pow(2.0, -RealSpec.BFloat16Spec.manW), 1.0,_,_), 1)
  acosTest(acosBF16Table, sqrtBF16Table, Some(acosExBF16Table), RealSpec.BFloat16Spec, n, r,
    "Test special value", generateSpecialValues(_,_), 1)


  val nOrderFP48 = 3
  val adrWFP48 = 10
  val extraBitsFP48 = 4
  val float48Spec = new RealSpec(12, 2047, 35)

  val acosFP48Table = ACosSim.acosTableGeneration(float48Spec,
    nOrderFP48, adrWFP48, float48Spec.manW, float48Spec.manW+extraBitsFP48)
  val sqrtFP48Table = ACosSim.sqrtTableGeneration(
    nOrderFP48, adrWFP48, float48Spec.manW, float48Spec.manW+extraBitsFP48)

  val smallValueFP48 = -7

  // acos(-|x|) is in [pi/2, pi]. this does not require super high resolution.
  acosTest(acosFP48Table, sqrtFP48Table, None, float48Spec, n, r,
    "Test close to 0.0:  [0.0, -2^-4]", generateRealWithin(-pow(2.0, smallValueFP48), 0.0,_,_), 3)
  acosTest(acosFP48Table, sqrtFP48Table, None, float48Spec, n, r,
    "Test normal table x < 0.5:  [-2^-4, -0.5]", generateRealWithin(-0.5, -pow(2.0, smallValueFP48),_,_), 3)
  acosTest(acosFP48Table, sqrtFP48Table, None, float48Spec, n, r,
    "Test normal table x > 0.5:  [-0.5, -1+2^-4]", generateRealWithin(-1.0+pow(2.0, smallValueFP48)-pow(2.0, -35), -0.5-pow(2.0, -35),_,_), 3)
  acosTest(acosFP48Table, sqrtFP48Table, None, float48Spec, n, r,
    "Test Puiseux: close to 1.0: [-1-2^+4, -1.0]", generateRealWithin(-1.0, -1.0+pow(2.0, smallValueFP48)-pow(2.0, -35),_,_), 3) // 2ULPs

  acosTest(acosFP48Table, sqrtFP48Table, None, float48Spec, n, r,
    "Test close to 0.0:  [0.0, 2^-4]", generateRealWithin(0.0, pow(2.0, smallValueFP48),_,_), 3)
  acosTest(acosFP48Table, sqrtFP48Table, None, float48Spec, n, r,
    "Test normal table x < 0.5:  [2^-4, 0.5]", generateRealWithin(pow(2.0, smallValueFP48), 0.5,_,_), 3)
  acosTest(acosFP48Table, sqrtFP48Table, None, float48Spec, n, r,
    "Test normal table x > 0.5:  [0.5, 1-2^-4]", generateRealWithin(0.5+pow(2.0, -35), 1.0-pow(2.0, smallValueFP48)+pow(2.0, -35),_,_), 31) // XXX
  acosTest(acosFP48Table, sqrtFP48Table, None, float48Spec, n, r,
    "Test Puiseux: close to 1.0: [1-2^-4, 1.0]", generateRealWithin(1.0-pow(2.0, smallValueFP48)+pow(2.0, -35), 1.0,_,_), 3) // 2ULPs

  acosTest(acosFP48Table, sqrtFP48Table, None, float48Spec, n, r,
    "Test special value", generateSpecialValues(_,_), 1)
}
