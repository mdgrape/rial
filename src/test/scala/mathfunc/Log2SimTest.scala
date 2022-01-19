//package rial.tests

//import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.FunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

//import org.scalatest.flatspec.AnyFlatSpec
//import org.scalatest.matchers.should.Matchers
//import scopt.OptionParser

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math.Log2Sim
import rial.mathfunc._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

class MathFuncLog2SimTest extends FunSuite with BeforeAndAfterAllConfigMap {
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

  def log2Test(t : FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"log2(x), format ${spec.toStringShort}, ${generatorStr}") {

      val log2 = (a:Double) => {log(a) / log(2.0)}

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

        val zi   = MathFuncLog2Sim.log2SimGeneric( t, x )
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r.toDouble)
        val errf = zi.toDouble - z0r.toDouble

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

            if(erri.abs>2.0) {
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
          assert(erri.abs<=tolerance || errf.abs <= pow(2.0, -spec.manW)) // Fixed point precision

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
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
      println(f"N=$n%d : 2LSB errors positive $err2lsbPos%d / negative $err2lsbNeg%d")
      println(f"N=$n%d : 2<   errors positive $errNlsbPos%d / negative $errNlsbNeg%d")
      println( "---------------------------------------------------------------")
    }
  }

  val log2F32TableI = Log2Sim.log2TableGeneration(
    2, 8, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+2,
    Some(Seq(27, 21, 20)), Some(Seq(27, 22, 20)))

  log2Test(log2F32TableI, RealSpec.Float32Spec, n, r,
    "Test Large More Than 1 [2, inf]", generateRealWithin(2.0, pow(2.0, 128.0),_,_), 2)
  log2Test(log2F32TableI, RealSpec.Float32Spec, n, r,
    "Test Small More Than 1 [1, 1+2^-12]",   generateRealWithin(1.0, 1.0+pow(2.0, -12),_,_), 2)
  log2Test(log2F32TableI, RealSpec.Float32Spec, n, r,
    "Test Small More Than 1 [1-2^-12, 1]",   generateRealWithin(1.0-pow(2.0, -12), 1.0,_,_), 2)
  log2Test(log2F32TableI, RealSpec.Float32Spec, n, r,
    "Test Small More Than 1 [1, 2]",   generateRealWithin(1.0, 2.0,_,_), 2)

  log2Test(log2F32TableI, RealSpec.Float32Spec, n, r,
    "Test Large Less Than 1 [0.5, 1]", generateRealWithin(0.5,1.0,_,_), 2)
  log2Test(log2F32TableI, RealSpec.Float32Spec, n, r,
    "Test Small Less Than 1 [0, 0.5]", generateRealWithin(0.0,0.5,_,_), 2)

  log2Test(log2F32TableI, RealSpec.Float32Spec, n, r,
    "Test Any Negative [-inf, 0]", generateRealWithin(-pow(2.0, 128), 0.0,_,_), 2)
}
