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

class MathFuncLog2SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def log2Test(t : FuncTableInt, tSmallPos : FuncTableInt, tSmallNeg : FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    toleranceUlps   : Int ) = {
    test(s"log2(x), format ${spec.toStringShort}, ${generatorStr}") {

      val log2 = (a:Double) => {log(a) / log(2.0)}
      val tolerance = pow(2, toleranceUlps) - 1

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

        val z0   = log2(x0)
        val z0r  = new RealGeneric(spec, z0)

        val zi   = MathFuncLogSim.logSimGeneric(/*islog2*/true, t, tSmallPos, tSmallNeg, x )
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r.toDouble)
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

            if(erri.abs>tolerance) {
              println(f"test: x   = ${x0}(${x.sgn}|${x.ex}(${x.ex-x.spec.exBias})|${x.man.toLong.toBinaryString})")
              println(f"test: ref = ${z0}%16g(${zrefsgn}|${zrefexp}(${(zrefexp-x.spec.exBias).toInt}%4d)|${zrefman.toLong.toBinaryString})")
              println(f"test: sim = ${zd}%16g(${zsimsgn}|${zsimexp}(${(zsimexp-x.spec.exBias).toInt}%4d)|${zsimman.toLong.toBinaryString})")
//               println(f"test: test(${zsimsgn}|${zsimexp}(${zsimexp - spec.exBias})|${zsimman.toLong.toBinaryString}(${zsimman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
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
          assert(erri.abs<=tolerance)// || errf.abs <= pow(2.0, -spec.manW)) // Fixed point precision

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
              + f"${zatMaxError} != ${log2(xatMaxError)}, "
              + f"diff = ${zatMaxError - log2(xatMaxError)}, x = ${xatMaxError}")
      }
      println(f"N=$n%d : +/- 1 errors(the last 1 bit) positive $err1lsbPos%d / negative $err1lsbNeg%d")
      println(f"N=$n%d : +/- 2 errors(the last 2 bit) positive $err2lsbPos%d / negative $err2lsbNeg%d")
      println(f"N=$n%d : +/- 3 errors(the last 2 bit) positive $errNlsbPos%d / negative $errNlsbNeg%d")
      println( "---------------------------------------------------------------")
    }
  }

  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3
  val log2FP32TableI              = MathFuncLogSim.logNormalTableGeneration       (RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)
  val log2FP32SmallPositiveTableI = MathFuncLogSim.logSmallPositiveTableGeneration(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)
  val log2FP32SmallNegativeTableI = MathFuncLogSim.logSmallNegativeTableGeneration(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)

  val taylorThresholdFP32 = MathFuncLogSim.calcTaylorThreshold(RealSpec.Float32Spec)

  log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Large More Than 1 [2, inf]", generateRealWithin(2.0, pow(2.0, 128.0),_,_), 1)
  log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    f"Test Small More Than 1 [1+2^-${taylorThresholdFP32}, 2]",   generateRealWithin(1.0+pow(2.0, -taylorThresholdFP32), 2.0,_,_), 2)

  log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    f"Test Taylor More Than 1 [1, 1+2^-${taylorThresholdFP32}]",   generateRealWithin(1.0, 1.0+pow(2.0, -taylorThresholdFP32) - pow(2.0,-RealSpec.Float32Spec.manW),_,_), 1)
  log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    f"Test Taylor More Than 1 [1-2^-${taylorThresholdFP32}, 1]",   generateRealWithin(1.0-pow(2.0, -taylorThresholdFP32) + pow(2.0,-RealSpec.Float32Spec.manW), 1.0,_,_), 1)

//   log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
//     "Test x = 2.0",   generateRealWithin(2.0, 2.0,_,_), 2)
//   log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
//     "Test x = 1.0",   generateRealWithin(1.0, 1.0,_,_), 2)
//   log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
//     "Test x = 0.5",   generateRealWithin(0.5, 0.5,_,_), 2)

//   val smallPositive = (-1 until -23 by -1).map( ex => {
//     val xmax = 1.0 + pow(2.0, ex)
//     val xmin = 1.0 + pow(2.0, ex-1)
//     log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
//       f"Test Small More Than 1 [1+2^${ex-1}%3d, 1+2^${ex}%3d]", generateRealWithin(xmin, xmax,_,_), 2)
//   })
//
//   val smallNegative = (-7 until -23 by -1).map( ex => {
//     val xmax = 1.0 - pow(2.0, ex-1)
//     val xmin = 1.0 - pow(2.0, ex)
//     log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
//       f"Test Small More Than 1 [1-2^${ex}%3d, 1-2^${ex-1}%3d]", generateRealWithin(xmin, xmax,_,_), 2)
//   })

  log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    f"Test Small Less Than 1 [0.5, 1-2^-${taylorThresholdFP32}]", generateRealWithin(0.5,1.0-pow(2.0, -taylorThresholdFP32),_,_), 2)
  log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Large Less Than 1 [0, 0.5]", generateRealWithin(0.0,0.5,_,_), 1)

  log2Test(log2FP32TableI, log2FP32SmallPositiveTableI, log2FP32SmallNegativeTableI, RealSpec.Float32Spec, n, r,
    "Test Any Negative [-inf, 0]", generateRealWithin(-pow(2.0, 128), 0.0,_,_), 1)


  val nOrderBF16    = 0
  val adrWBF16      = 7
  val extraBitsBF16 = 1
  val log2BF16TableI              = MathFuncLogSim.logNormalTableGeneration       (RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16)
  val log2BF16SmallPositiveTableI = MathFuncLogSim.logSmallPositiveTableGeneration(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16)
  val log2BF16SmallNegativeTableI = MathFuncLogSim.logSmallNegativeTableGeneration(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16)

  val taylorThresholdBF16 = MathFuncLogSim.calcTaylorThreshold(RealSpec.BFloat16Spec)

  log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Large More Than 1 [2, inf]", generateRealWithin(2.0, pow(2.0, 128.0),_,_), 1)
  log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    f"Test Small More Than 1 [1+2^-${taylorThresholdBF16}, 2]",   generateRealWithin(1.0+pow(2.0, -taylorThresholdBF16), 2.0,_,_), 2)

  log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    f"Test Taylor More Than 1 [1, 1+2^-${taylorThresholdBF16}]",   generateRealWithin(1.0, 1.0+pow(2.0, -taylorThresholdBF16) - pow(2.0,-RealSpec.BFloat16Spec.manW),_,_), 1)
  log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    f"Test Taylor More Than 1 [1-2^-${taylorThresholdBF16}, 1]",   generateRealWithin(1.0-pow(2.0, -taylorThresholdBF16) + pow(2.0,-RealSpec.BFloat16Spec.manW), 1.0,_,_), 2) // XXX if we set extraBits = 2, we can reduce the error

//   log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
//     "Test x = 2.0",   generateRealWithin(2.0, 2.0,_,_), 2)
//   log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
//     "Test x = 1.0",   generateRealWithin(1.0, 1.0,_,_), 2)
//   log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
//     "Test x = 0.5",   generateRealWithin(0.5, 0.5,_,_), 2)

//   val smallPositive = (-1 until -23 by -1).map( ex => {
//     val xmax = 1.0 + pow(2.0, ex)
//     val xmin = 1.0 + pow(2.0, ex-1)
//     log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
//       f"Test Small More Than 1 [1+2^${ex-1}%3d, 1+2^${ex}%3d]", generateRealWithin(xmin, xmax,_,_), 2)
//   })
//
//   val smallNegative = (-7 until -23 by -1).map( ex => {
//     val xmax = 1.0 - pow(2.0, ex-1)
//     val xmin = 1.0 - pow(2.0, ex)
//     log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
//       f"Test Small More Than 1 [1-2^${ex}%3d, 1-2^${ex-1}%3d]", generateRealWithin(xmin, xmax,_,_), 2)
//   })

  log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    f"Test Small Less Than 1 [0.5, 1-2^-${taylorThresholdBF16}]", generateRealWithin(0.5,1.0-pow(2.0, -taylorThresholdBF16),_,_), 2)
  log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Large Less Than 1 [0, 0.5]", generateRealWithin(0.0,0.5,_,_), 1)

  log2Test(log2BF16TableI, log2BF16SmallPositiveTableI, log2BF16SmallNegativeTableI, RealSpec.BFloat16Spec, n, r,
    "Test Any Negative [-inf, 0]", generateRealWithin(-pow(2.0, 128), 0.0,_,_), 1)

}
