import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.mathfunc._
import rial.util.ScalaUtil._
import rial.testUtil.ScalaTestUtil._ // errorLSB
import rial.arith._
import rial.table._

import com.sun.jna._
trait libc extends Library {
  def cosf(x: Float):Float
}

class MathFuncCosSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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

  def cosTest(ts : Seq[FuncTableInt], spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"cos(x), format ${spec.toStringShort}, ${generatorStr}") {

//       val libc = Native.loadLibrary("c", classOf[libc]).asInstanceOf[libc]

      var maxError    = SafeLong(0)
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val z0  = cos(x0)
        val z0r  = new RealGeneric(spec, z0)
//         val z0f  = libc.cosf(x.toFloat)
//         assert(z0f.toDouble - pow(2.0, -22) < z0 && z0 < z0f.toDouble + pow(2.0, -22))

        val zi   = MathFuncCosSim.cosSimGeneric( ts, x )
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
              + f"${zatMaxError} != ${cos(xatMaxError)}, "
              + f"diff = ${zatMaxError - cos(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- $k errors positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val nOrder = 2
  val adrW = 8
  val extraBits = 3
  val sinF32TableI = MathFuncSinSim.sinTableGeneration(
    nOrder, adrW, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+extraBits)

  //XXX allowing error in 2ULPs

  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-10pi, -2pi]", generateRealWithin(-10 * Pi, -2 * Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 3)
  cosTest(sinF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [2pi, 10pi]", generateRealWithin(2.0 * Pi, 10.0 * Pi,_,_), 3)
}
