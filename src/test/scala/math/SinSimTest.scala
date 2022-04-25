import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

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

class SinSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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

  def sinTest(ts : Seq[FuncTableInt], taylorOrder: Int, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int
    ) = {
    test(s"sin(x), format ${spec.toStringShort}, ${generatorStr}") {

      var maxError    = SafeLong(0)
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val z0   = sin(x0)
        val z0r  = new RealGeneric(spec, z0)

        val zi   = SinCosSim.sincosSimGeneric(true, ts, x, taylorOrder )
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r).toInt

        if (x0.isInfinity) {
          if(0 < x0) {
            assert(zi.isInfinite)
          } else {
            assert(zi.isZero)
          }
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if(erri.abs > tolerance) {
            val xsgn = bit(spec.W-1, x.value).toInt
            val xexp = slice(spec.manW, spec.exW, x.value)
            val xman = x.value & maskSL(spec.manW)

            val zsimsgn = bit(spec.W-1, zi.value).toInt
            val zsimexp = slice(spec.manW, spec.exW, zi.value)
            val zsimman = zi.value & maskSL(spec.manW)

            val zrefsgn = bit(spec.W-1, z0r.value).toInt
            val zrefexp = slice(spec.manW, spec.exW, z0r.value)
            val zrefman = z0r.value & maskSL(spec.manW)

            println(f"test: x   = ${x0}%16g(${x.sgn}|${x.ex}(${x.ex-x.spec.exBias})|${x.man.toLong.toBinaryString})")
            println(f"test: ref = ${z0}%16g(${zrefsgn}|${zrefexp}(${zrefexp-x.spec.exBias})|${zrefman.toLong.toBinaryString})")
            println(f"test: sim = ${zd}%16g(${zsimsgn}|${zsimexp}(${zsimexp-x.spec.exBias})|${zsimman.toLong.toBinaryString})")
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

          assert(erri.abs <= tolerance)

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

  val taylorOrder5th = 5
  val taylorOrder3rd = 3

  val nOrderFP32     = 2
  val adrWFP32       = 8
  val extraBitsFP32  = 3

  val sinFP32TableItaylor5th = SinCosSim.sincosTableGeneration(
    nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32,
    None, None, taylorOrder5th)

  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [-10pi, -2pi]", generateRealWithin(-10 * Pi, -2 * Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [2pi, 10pi]", generateRealWithin(2.0 * Pi, 10.0 * Pi,_,_), 3)

  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [16pi, 32pi]", generateRealWithin(16.0 * Pi, 32.0 * Pi,_,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [32pi, 64pi]", generateRealWithin(32.0 * Pi, 64.0 * Pi,_,_), 3)

  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [-32pi, -16pi]", generateRealWithin( -32.0 * Pi,-16.0 * Pi, _,_), 3)
  sinTest(sinFP32TableItaylor5th, taylorOrder5th, RealSpec.Float32Spec, n, r,
    "Test 5th-order Within [-64pi, -32pi]", generateRealWithin( -64.0 * Pi,-32.0 * Pi, _,_), 3)

  val sinFP32TableItaylor3rd = SinCosSim.sincosTableGeneration(
    nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32,
    None, None, taylorOrder3rd)

  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [-10pi, -2pi]", generateRealWithin(-10 * Pi, -2 * Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [2pi, 10pi]", generateRealWithin(2.0 * Pi, 10.0 * Pi,_,_), 3)

  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [16pi, 32pi]", generateRealWithin(16.0 * Pi, 32.0 * Pi,_,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [32pi, 64pi]", generateRealWithin(32.0 * Pi, 64.0 * Pi,_,_), 3)

  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [-32pi, -16pi]", generateRealWithin( -32.0 * Pi,-16.0 * Pi, _,_), 3)
  sinTest(sinFP32TableItaylor3rd, taylorOrder3rd, RealSpec.Float32Spec, n, r,
    "Test 3rd-order Within [-64pi, -32pi]", generateRealWithin( -64.0 * Pi,-32.0 * Pi, _,_), 3)

  val nOrderBF16     = 0
  val adrWBF16       = 7
  val extraBitsBF16  = 1

  val sinBF16TableItaylor3rd = SinCosSim.sincosTableGeneration(
    nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16,
    None, None, taylorOrder3rd)

  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [-10pi, -2pi]", generateRealWithin(-10 * Pi, -2 * Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [2pi, 10pi]", generateRealWithin(2.0 * Pi, 10.0 * Pi,_,_), 3)

  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [16pi, 32pi]", generateRealWithin(16.0 * Pi, 32.0 * Pi,_,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [32pi, 64pi]", generateRealWithin(32.0 * Pi, 64.0 * Pi,_,_), 3)

  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [-32pi, -16pi]", generateRealWithin( -32.0 * Pi,-16.0 * Pi, _,_), 3)
  sinTest(sinBF16TableItaylor3rd, taylorOrder3rd, RealSpec.BFloat16Spec, n, r,
    "Test 3rd-order Within [-64pi, -32pi]", generateRealWithin( -64.0 * Pi,-32.0 * Pi, _,_), 3)

  val float48Spec = new RealSpec(10, 511, 37)

  val nOrderFP48     = 3
  val adrWFP48       = 10
  val extraBitsFP48  = 4

  val sinFP48TableItaylor3rd = SinCosSim.sincosTableGeneration(
    nOrderFP48, adrWFP48, float48Spec.manW, float48Spec.manW+extraBitsFP48,
    None, None, taylorOrder3rd)

  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [-10pi, -2pi]", generateRealWithin(-10 * Pi, -2 * Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 3)
  sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
    "Test 3rd-order Within [2pi, 10pi]", generateRealWithin(2.0 * Pi, 10.0 * Pi,_,_), 3)

//   sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
//     "Test 3rd-order Within [16pi, 32pi]", generateRealWithin(16.0 * Pi, 32.0 * Pi,_,_), 3)
//   sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
//     "Test 3rd-order Within [32pi, 64pi]", generateRealWithin(32.0 * Pi, 64.0 * Pi,_,_), 3)

//   sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
//     "Test 3rd-order Within [-32pi, -16pi]", generateRealWithin( -32.0 * Pi,-16.0 * Pi, _,_), 3)
//   sinTest(sinFP48TableItaylor3rd, taylorOrder3rd, float48Spec, n, r,
//     "Test 3rd-order Within [-64pi, -32pi]", generateRealWithin( -64.0 * Pi,-32.0 * Pi, _,_), 3)

  val nOrderFP64     = 3
  val adrWFP64       = 12
  val extraBitsFP64  = 4

  val sinFP64TableItaylor3rd = SinCosSim.sincosTableGeneration(
    nOrderFP64, adrWFP64, RealSpec.Float64Spec.manW, RealSpec.Float64Spec.manW+extraBitsFP64,
    None, None, taylorOrder3rd)

  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [-10pi, -2pi]", generateRealWithin(-10 * Pi, -2 * Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 7)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test 3rd-order Within [2pi, 10pi]", generateRealWithin(2.0 * Pi, 10.0 * Pi,_,_), 7)

//   sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
//     "Test 3rd-order Within [16pi, 32pi]", generateRealWithin(16.0 * Pi, 32.0 * Pi,_,_), 7)
//   sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
//     "Test 3rd-order Within [32pi, 64pi]", generateRealWithin(32.0 * Pi, 64.0 * Pi,_,_), 7)
//
//   sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
//     "Test 3rd-order Within [-32pi, -16pi]", generateRealWithin( -32.0 * Pi,-16.0 * Pi, _,_), 7)
//   sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
//     "Test 3rd-order Within [-64pi, -32pi]", generateRealWithin( -64.0 * Pi,-32.0 * Pi, _,_), 7)

  val delta = pow(2.0, -32)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test sin around 0", generateRealWithin(-delta*Pi, delta*Pi,_,_), 8)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test sin around pi/2", generateRealWithin((0.5-delta)*Pi, (0.5+delta)*Pi,_,_), 8)
// XXX since sin(x) is evaluated by converting x into the range [0, pi/2],
//     sin(pi - eps) is converted to sin(eps). So, if x is close to pi,
//     cancellation error extremely reduces the precision of the result.
//     So the following test case cannot gain as much precision as others.
//     If eps is enough small, the result from taylor path will be selected,
//     but sin(x) = x - x^3/pi cannot gain precision if x^2/pi < 2^-manW
//     because of numerical error.
//   sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
//     "Test sin around pi", generateRealWithin((1-delta)*Pi, (1+delta)*Pi,_,_), 8)
  sinTest(sinFP64TableItaylor3rd, taylorOrder3rd, RealSpec.Float64Spec, n, r,
    "Test sin around 3pi/2", generateRealWithin((1.5-delta)*Pi, (1.5+delta)*Pi,_,_), 8)
}
