
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

import rial.math._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

// This test checks (mantissas exactly match) **OR** (the error is less than 2^-manW).
// it has an error in the last 5 bits around acos(pi) and acos(2pi).
class ACosSimTest extends FunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
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

  def acosTest(t: Seq[FuncTableInt], spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric),
    tolerance : Int ) = {
    test(s"acos(x), format ${spec.toStringShort}, ${generatorStr}") {
      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      var err1lsbPos = 0
      var err1lsbNeg = 0
      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = math.acos(x0)
        val z0r  = new RealGeneric(spec, z0)
        val zi   = ACosSim.acosSimGeneric( t, x )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble)
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if (z0r.value != zi.value) {
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
          println(f"test: test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman.toLong.toBinaryString}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
        }
        if (x0.isInfinity) {
          assert(zi.isNaN)
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if (erri.abs>=2.0) {
            println(f"Error more than 2 LSB : x = ${x.toDouble}%14.7e : refz = $z0%14.7e simz = ${zi.toDouble}%14.7e $errf%14.7e $erri%f")
          } else if (erri>=1.0) {
            err1lsbPos+=1
          }
          else if (erri<= -1.0) {
            err1lsbNeg+=1
          }
          assert(erri.abs<=tolerance.toDouble || (z0 - zd).abs < tolerance * pow(2.0, -spec.manW)) // XXX

          if (maxError < erri.abs) {
            maxError = erri.abs
            xatMaxError = x0
            zatMaxError = zd
          }
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is ${zatMaxError} != ${math.acos(xatMaxError)} (acos(${xatMaxError})), diff = ${zatMaxError - math.acos(xatMaxError)}")
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
    }
  }

  val acosBF16TableI = ACosSim.acosTableGeneration( 0, 7, 7, 7+2)

  acosTest(acosBF16TableI, RealSpec.BFloat16Spec, n, r,
     "Test Within [-1, -2^-8]", generateRealWithin(-0.96, -pow(2.0, -8),_,_), 1)
  acosTest(acosBF16TableI, RealSpec.BFloat16Spec, n, r,
     "Test Within [-2^-8, 2^-23]",  generateRealWithin(-pow(2.0, -8),-pow(2.0, -23),_,_), 1)
  acosTest(acosBF16TableI, RealSpec.BFloat16Spec, n, r,
     "Test Within [-2^-23, 0]",  generateRealWithin(-pow(2.0, -23), 0.0,_,_), 1)
  acosTest(acosBF16TableI, RealSpec.BFloat16Spec, n, r,
     "Test Within [ 0,  2^-23]", generateRealWithin(0.0,pow(2.0, -23),_,_), 1)
  acosTest(acosBF16TableI, RealSpec.BFloat16Spec, n, r,
     "Test Within [ 2^-23,  2^-8]", generateRealWithin(pow(2.0,-23),pow(2.0, -8),_,_), 1)
  acosTest(acosBF16TableI, RealSpec.BFloat16Spec, n, r,
     "Test Within [ 2^-8,  0.96]",    generateRealWithin(pow(2.0, -8), 0.96,_,_), 1)

  val acosF32TableI = ACosSim.acosTableGeneration( 2, 8, 23, 23+2 )

  acosTest(acosF32TableI, RealSpec.Float32Spec, n, r,
     "Test Within [-0.96, -2^-8]", generateRealWithin(-0.96, -pow(2.0, -8),_,_), 2)
  acosTest(acosF32TableI, RealSpec.Float32Spec, n, r,
     "Test Within [-2^-8, 2^-23]",  generateRealWithin(-pow(2.0, -8),-pow(2.0, -23),_,_), 1)
  acosTest(acosF32TableI, RealSpec.Float32Spec, n, r,
     "Test Within [-2^-23, 0]",  generateRealWithin(-pow(2.0, -23), 0.0,_,_), 1)
  acosTest(acosF32TableI, RealSpec.Float32Spec, n, r,
     "Test Within [ 0,  2^-23]", generateRealWithin(0.0,pow(2.0, -23),_,_), 1)
  acosTest(acosF32TableI, RealSpec.Float32Spec, n, r,
     "Test Within [ 2^-23,  2^-8]", generateRealWithin(pow(2.0,-23),pow(2.0, -8),_,_), 1)
  acosTest(acosF32TableI, RealSpec.Float32Spec, n, r,
     "Test Within [ 2^-8,  0.96]",    generateRealWithin(pow(2.0, -8), 0.96,_,_), 2)
  // XXX currently the errors in the range in 0.96 ~ 1.0 are too big
}
