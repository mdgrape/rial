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

class MathFuncLogSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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

  def logTest(t : FuncTableInt, tSmallPos : FuncTableInt, tSmallNeg : FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"log(x), format ${spec.toStringShort}, ${generatorStr}") {

      var maxError    = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val z0   = log(x0)
        val z0r  = new RealGeneric(spec, z0)

        val zi   = MathFuncLogSim.logSimGeneric( /*islog2*/false, t, tSmallPos, tSmallNeg, x )
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toInt
        val errf = zi.toDouble - z0r.toDouble

        if (x0.isInfinity) {
          if(0 < x0) {
            assert(zi.isInfinite)
          } else {
            assert(zi.isNaN)
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
              + f"${zatMaxError} != ${log(xatMaxError)}, "
              + f"diff = ${zatMaxError - log(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- $k errors positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3
  val log2FP32TableI              = MathFuncLogSim.logNormalTableGeneration       (RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)
  val log2FP32SmallPositiveTableI = MathFuncLogSim.logSmallPositiveTableGeneration(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)
  val log2FP32SmallNegativeTableI = MathFuncLogSim.logSmallNegativeTableGeneration(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)

  logTest(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Large More Than 1 [2, inf]", generateRealWithin(2.0, pow(2.0, 128.0),_,_), 3)
  logTest(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Small More Than 1 [1, 1+2^-8]",   generateRealWithin(1.0, 1.0+pow(2.0, -8) - pow(2.0,-23),_,_), 3)
  logTest(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Small More Than 1 [1-2^-8, 1]",   generateRealWithin(1.0-pow(2.0, -8) + pow(2.0,-23), 1.0,_,_), 3)
  logTest(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Small More Than 1 [1, 2]",   generateRealWithin(1.0+pow(2.0, -8), 2.0,_,_), 4) // XXX not in 2 ULPs

//   val smallPositive = (-1 until -23 by -1).map( ex => {
//     val xmax = 1.0 + pow(2.0, ex)
//     val xmin = 1.0 + pow(2.0, ex-1)
//     logTest(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
//       f"Test Small More Than 1 [1+2^${ex-1}%3d, 1+2^${ex}%3d]", generateRealWithin(xmin, xmax,_,_), 2)
//   })
//
//   val smallNegative = (-1 until -23 by -1).map( ex => {
//     val xmax = 1.0 - pow(2.0, ex-1)
//     val xmin = 1.0 - pow(2.0, ex)
//     logTest(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
//       f"Test Small More Than 1 [1-2^${ex}%3d, 1-2^${ex-1}%3d]", generateRealWithin(xmin, xmax,_,_), 2)
//   })

  logTest(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Large Less Than 1 [0.5, 1]", generateRealWithin(0.5,1.0,_,_), 3)
  logTest(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Small Less Than 1 [0, 0.5]", generateRealWithin(0.0,0.5,_,_), 3)

  logTest(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Any Negative [-inf, 0]", generateRealWithin(-pow(2.0, 128), 0.0,_,_), 1)

  val nOrderBF16    = 0
  val adrWBF16      = 7
  val extraBitsBF16 = 1
  val log2BF16TableI              = MathFuncLogSim.logNormalTableGeneration       (RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16)
  val log2BF16SmallPositiveTableI = MathFuncLogSim.logSmallPositiveTableGeneration(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16)
  val log2BF16SmallNegativeTableI = MathFuncLogSim.logSmallNegativeTableGeneration(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16)

  logTest(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Large More Than 1 [2, inf]", generateRealWithin(2.0, pow(2.0, 128.0),_,_), 1)
  logTest(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Small More Than 1 [1, 1+2^-8]",   generateRealWithin(1.0, 1.0+pow(2.0, -8) - pow(2.0,-23),_,_), 1)
  logTest(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Small More Than 1 [1-2^-8, 1]",   generateRealWithin(1.0-pow(2.0, -8) + pow(2.0,-23), 1.0,_,_), 1)
  logTest(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Small More Than 1 [1, 2]",   generateRealWithin(1.0+pow(2.0, -8), 2.0,_,_), 3)

//   val smallPositive = (-1 until -23 by -1).map( ex => {
//     val xmax = 1.0 + pow(2.0, ex)
//     val xmin = 1.0 + pow(2.0, ex-1)
//     logTest(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
//       f"Test Small More Than 1 [1+2^${ex-1}%3d, 1+2^${ex}%3d]", generateRealWithin(xmin, xmax,_,_), 2)
//   })
//
//   val smallNegative = (-1 until -23 by -1).map( ex => {
//     val xmax = 1.0 - pow(2.0, ex-1)
//     val xmin = 1.0 - pow(2.0, ex)
//     logTest(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
//       f"Test Small More Than 1 [1-2^${ex}%3d, 1-2^${ex-1}%3d]", generateRealWithin(xmin, xmax,_,_), 2)
//   })

  logTest(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Large Less Than 1 [0.5, 1]", generateRealWithin(0.5,1.0,_,_), 3)
  logTest(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Small Less Than 1 [0, 0.5]", generateRealWithin(0.0,0.5,_,_), 1)

  logTest(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Any Negative [-inf, 0]", generateRealWithin(-pow(2.0, 128), 0.0,_,_), 1)

}
