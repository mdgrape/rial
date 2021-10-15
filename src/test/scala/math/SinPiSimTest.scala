
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
// it has an error in the last 5 bits around sin(pi) and sin(2pi).
class SinPiSimTest extends FunSuite with BeforeAndAfterAllConfigMap {
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

  def sinPiTest(t: Seq[FuncTableInt], spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    test(s"sinPi(x), format ${spec.toStringShort}, ${generatorStr}") {
      var err1lsbPos = 0
      var err1lsbNeg = 0
      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = math.sin(math.Pi * x0)
        val z0r  = new RealGeneric(spec, z0)
        val zi   = SinPiSim.sinPiSimGeneric( t, x )
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

//           println(f"x   = ${x0})
//           println(f"ref = ${z0}")
//           println(f"sim = ${zd}")
//           println(f"test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman.toLong.toBinaryString}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
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
          assert(erri.abs<=16.0 || (z0 - zd).abs < pow(2.0, -spec.manW)) // XXX
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
    }
  }

  val sinPiF32TableI = SinPiSim.sinPiTableGeneration( 2, 8, 23, 23+6 ) // +6 to suppress error in [2^-12, 2^-1)

  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
     "Test Within [-1, 0]", generateRealWithin(-1.0, 0.0,_,_))

  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [0, 2^-12]", generateRealWithin(0.0, pow(2.0, -12)-pow(2.0, -35),_,_))
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [2^-12, 2^-1]", generateRealWithin(pow(2.0, -12), 0.5,_,_))
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [2^-1, 1]", generateRealWithin(0.5, 1.0,_,_))
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [1, 2]", generateRealWithin(1.0, 2.0,_,_))

//   val sinPiBF16TableI = SinPiSim.sinPiTableGeneration(0, 7, 7, 7 )
//
//   sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Within (-128,128)",generateRealWithin(128.0,_,_))
//   sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test All range",generateRealFull(_,_) )
}
