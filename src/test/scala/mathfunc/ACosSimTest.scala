import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.mathfunc._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

class MathFuncACosSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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

  def acosTest(t : Seq[FuncTableInt], tSqrt : FuncTableInt,
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

        val zi   = MathFuncACosSim.acosSimGeneric(0, t, tSqrt, x )
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

  val acosFP32Table = MathFuncACosSim.acosTableGeneration(RealSpec.Float32Spec,
    nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)
  val sqrtFP32Table = MathFuncACosSim.sqrtTableGeneration(
    nOrderFP32, adrWFP32, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBitsFP32)

  // acos(-|x|) is in [pi/2, pi]. this does not require super high resolution.
  acosTest(acosFP32Table, sqrtFP32Table, RealSpec.Float32Spec, n, r,
    "Test Taylor: close to 0.0:  [0.0, -2^-4]", generateRealWithin(-pow(2.0, -4), 0.0,_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, RealSpec.Float32Spec, n, r,
    "Test normal table x < 0.5:  [-2^-4, -0.5]", generateRealWithin(-0.5, -pow(2.0, -4),_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, RealSpec.Float32Spec, n, r,
    "Test normal table x > 0.5:  [-0.5, -1+2^-4]", generateRealWithin(-1.0+pow(2.0, -4)-pow(2.0, -23), -0.5-pow(2.0, -23),_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, RealSpec.Float32Spec, n, r,
    "Test Puiseux: close to 1.0: [-1-2^+4, -1.0]", generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_), 3) // 2ULPs

  acosTest(acosFP32Table, sqrtFP32Table, RealSpec.Float32Spec, n, r,
    "Test Taylor: close to 0.0:  [0.0, 2^-4]", generateRealWithin(0.0, pow(2.0, -4),_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, RealSpec.Float32Spec, n, r,
    "Test normal table x < 0.5:  [2^-4, 0.5]", generateRealWithin(pow(2.0, -4), 0.5,_,_), 1)
  acosTest(acosFP32Table, sqrtFP32Table, RealSpec.Float32Spec, n, r,
    "Test normal table x > 0.5:  [0.5, 1-2^-4]", generateRealWithin(0.5+pow(2.0, -23), 1.0-pow(2.0, -4)+pow(2.0, -23),_,_), 3) // 2ULPs
  acosTest(acosFP32Table, sqrtFP32Table, RealSpec.Float32Spec, n, r,
    "Test Puiseux: close to 1.0: [1-2^-4, 1.0]", generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_), 3) // 2ULPs


//   val nOrderBF16    = 0
//   val adrWBF16      = 7
//   val extraBitsBF16 = 1
// 
//   val acosBF16Table = MathFuncACosSim.acosTableGeneration(
//     nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)
//   val sqrtBF16Table = MathFuncACosSim.sqrtTableGeneration(
//     nOrderBF16, adrWBF16, RealSpec.BFloat16Spec.manW, RealSpec.BFloat16Spec.manW+extraBitsBF16)
// 
//   val taylorThresholdBF16 = MathFuncACosSim.calcTaylorThreshold(RealSpec.BFloat16Spec.manW, 3)
// 
//   acosTest(acosBF16Table, sqrtBF16Table, RealSpec.BFloat16Spec, n, r,
//     f"Test Taylor: close to 0.0:  [0.0, -2^${taylorThresholdBF16}]", generateRealWithin(-pow(2.0, taylorThresholdBF16), 0.0,_,_), 1)
//   acosTest(acosBF16Table, sqrtBF16Table, RealSpec.BFloat16Spec, n, r,
//     f"Test normal table x < 0.5:  [-2^${taylorThresholdBF16}, -0.5]", generateRealWithin(-0.5, -pow(2.0, taylorThresholdBF16),_,_), 1)
//   acosTest(acosBF16Table, sqrtBF16Table, RealSpec.BFloat16Spec, n, r,
//     f"Test normal table x > 0.5:  [-0.5, -1+2^${taylorThresholdBF16}]", generateRealWithin(-1.0+pow(2.0, taylorThresholdBF16)-pow(2.0, -RealSpec.BFloat16Spec.manW), -0.5-pow(2.0, -RealSpec.BFloat16Spec.manW),_,_), 1)
//   acosTest(acosBF16Table, sqrtBF16Table, RealSpec.BFloat16Spec, n, r,
//     f"Test Puiseux: close to 1.0: [-1-2^${taylorThresholdBF16}, -1.0]", generateRealWithin(-1.0, -1.0+pow(2.0, taylorThresholdBF16)-pow(2.0, -RealSpec.BFloat16Spec.manW),_,_), 1)
// 
//   acosTest(acosBF16Table, sqrtBF16Table, RealSpec.BFloat16Spec, n, r,
//     f"Test Taylor: close to 0.0:  [0.0, 2^${taylorThresholdBF16}]", generateRealWithin(0.0, pow(2.0, taylorThresholdBF16),_,_), 1)
//   acosTest(acosBF16Table, sqrtBF16Table, RealSpec.BFloat16Spec, n, r,
//     f"Test normal table x < 0.5:  [2^${taylorThresholdBF16}, 0.5]", generateRealWithin(pow(2.0, taylorThresholdBF16), 0.5,_,_), 1)
//   acosTest(acosBF16Table, sqrtBF16Table, RealSpec.BFloat16Spec, n, r,
//     f"Test normal table x > 0.5:  [0.5, 1-2^${taylorThresholdBF16}]", generateRealWithin(0.5+pow(2.0, -RealSpec.BFloat16Spec.manW), 1.0-pow(2.0, taylorThresholdBF16)+pow(2.0, -RealSpec.BFloat16Spec.manW),_,_), 1)
//   acosTest(acosBF16Table, sqrtBF16Table, RealSpec.BFloat16Spec, n, r,
//     f"Test Puiseux: close to 1.0: [1-2^${taylorThresholdBF16}, 1.0]", generateRealWithin(1.0-pow(2.0, taylorThresholdBF16)+pow(2.0, -RealSpec.BFloat16Spec.manW), 1.0,_,_), 1)

}
