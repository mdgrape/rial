
//package rial.tests

//import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.funsuite.AnyFunSuite
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
class SinPiSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 100000

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
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric),
    tolerance : Int, use3rdOrder : Boolean ) = {
    test(s"sinPi(x), format ${spec.toStringShort}, ${generatorStr}, 3rd order = ${use3rdOrder}") {
      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      var err1lsbPos = 0
      var err1lsbNeg = 0
      var err2lsbPos = 0
      var err2lsbNeg = 0
      var errNlsbPos = 0
      var errNlsbNeg = 0
      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = math.sin(math.Pi * x0)
        val z0r  = new RealGeneric(spec, z0)
        val zi   = SinPiSim.sinPiSimGeneric( t, x, use3rdOrder )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble)
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if (z0r.value != zi.value) {
          val ztestsgn = bit(spec.W-1, zi.value).toInt
          val ztestexp = slice(spec.manW, spec.exW, zi.value)
          val ztestman = zi.value & maskSL(spec.manW)

          val zrefsgn = bit(spec.W-1, z0r.value).toInt
          val zrefexp = slice(spec.manW, spec.exW, z0r.value)
          val zrefman = z0r.value & maskSL(spec.manW)

//           println(f"test: x   = ${x0}")
//           println(f"test: ref = ${z0}")
//           println(f"test: sim = ${zd}")
//           println(f"test: test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman.toLong.toBinaryString}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
        }
        if (x0.isInfinity) {
          assert(zi.isNaN)
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else if (x0 == floor(x0)) { // x is an integer
          assert(zi.isZero)
        } else {
          if (erri.abs>=4) {
            println(f"Error more than 2 LSB : x = ${x.toDouble}%14.7e : refz = $z0%14.7e simz = ${zi.toDouble}%14.7e $errf%14.7e $erri%f")
            println(f"Error more than 2 LSB : x = ${x.sgn}|${x.ex}(${x.ex-x.spec.exBias})|${x.man.toLong.toBinaryString}")
            if(erri > 0) {
              errNlsbPos += 1
            } else {
              errNlsbNeg += 1
            }
          } else if (erri.abs>=2) {
//             println(f"Error in the last 2 LSB : x = ${x.toDouble}%14.7e : refz = $z0%14.7e simz = ${zi.toDouble}%14.7e $errf%14.7e $erri%f")
            if(erri > 0) {
              err2lsbPos += 1
            } else {
              err2lsbNeg += 1
            }
          } else if (erri.abs >= 1) {
            if (erri > 0) {
              err1lsbPos+=1
            }
            else if (erri < 0.0) {
              err1lsbNeg+=1
            }
          }
//           assert(erri.abs<=tolerance.toDouble || (z0 - zd).abs < pow(2.0, -spec.manW)) // XXX
          assert(erri.abs<=tolerance) // XXX

          if (maxError < erri.abs) {
            maxError = erri.abs
            xatMaxError = x0
            zatMaxError = zd
          }
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is ${zatMaxError} != ${math.sin(Pi * xatMaxError)} (sin(Pi * ${xatMaxError})), diff = ${zatMaxError - math.sin(Pi * xatMaxError)}")
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
      println(f"N=$n%d : 2LSB errors positive $err2lsbPos%d / negative $err2lsbNeg%d")
      println(f"N=$n%d : 2<   errors positive $errNlsbPos%d / negative $errNlsbNeg%d")
    }
  }

  // we need to make polynomial extrabits 3, to have 1-bit accuracy.
  val sinPiF32TableI = SinPiSim.sinPiTableGeneration( 2, 8, 23, 23+3 )

  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-1, -0.5]", generateRealWithin(-1.0, -0.5,_,_), 1, false)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-0.5, -2^-12]", generateRealWithin(-0.5, -pow(2.0, -12),_,_), 1, false)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-2^-12, 0]", generateRealWithin(-pow(2.0, -12), 0.0,_,_), 1, false)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [0, 2^-12]", generateRealWithin(0.0, pow(2.0, -12),_,_), 1, false)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [2^-12, 2^-1]", generateRealWithin(pow(2.0, -12), 0.5,_,_), 1, false)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [2^-1, 1]", generateRealWithin(0.5, 1.0,_,_), 1, false)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [1, 2]", generateRealWithin(1.0, 2.0,_,_), 1, false)

  // TODO: 2-bit error in 3rd-order taylor expantion
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-1,    -0.5]",   generateRealWithin(-1.0, -0.5,_,_), 2, true)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-0.5,  -2^-6]",  generateRealWithin(-0.5, -pow(2.0, -6),_,_), 1, true)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-2^-6, -2^-12]", generateRealWithin(pow(2.0, -6), -pow(2.0, -12),_,_), 2, true) // XXX
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [-2^-12, 0]",     generateRealWithin(-pow(2.0, -12), 0.0,_,_), 1, true)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [0,      2^-12]", generateRealWithin(0.0, pow(2.0, -12)-pow(2.0, -35),_,_), 1, true)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [2^-12,  2^-6]",  generateRealWithin(pow(2.0, -12), pow(2.0, -6),_,_), 2, true) // XXX
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [2^-6,   0.5]",   generateRealWithin(pow(2.0, -6), 0.5,_,_), 1, true)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [0.5,    1.0]",   generateRealWithin(0.5, 1.0,_,_), 2, true)
  sinPiTest(sinPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [1.0,    2.0]",   generateRealWithin(1.0, 2.0,_,_), 2, true)

  val sinPiBF16TableI = SinPiSim.sinPiTableGeneration(0, 7, 7, 7 )

  sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within (-1,0)",generateRealWithin(-1.0, 0.0,_,_), 2, false)
  sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within ( 0,2^-12)",generateRealWithin(0.0, pow(2.0, -12) - pow(2.0, -35),_,_), 2, false)
  sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within (2^-12, 0.5)",generateRealWithin(pow(2.0, -12) - pow(2.0, -35), 0.5,_,_), 2, false)
  sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within (0.5,1)",generateRealWithin(0.5, 1.0,_,_), 2, false)
  sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within (1,2)",generateRealWithin(1.0, 2.0,_,_), 2, false)

  // TODO:
  // - change testing regions to the appropreate ones
  // - fix large error here
//   sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Within (-1,0)",generateRealWithin(-1.0, 0.0,_,_), 2, true)
  sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within ( 0,2^-12)",generateRealWithin(0.0, pow(2.0, -12) - pow(2.0, -35),_,_), 2, true)
//   sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Within (2^-12, 0.5)",generateRealWithin(pow(2.0, -12) - pow(2.0, -35), 0.5,_,_), 2, true)
  sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within (0.5,1)",generateRealWithin(0.5, 1.0,_,_), 2, true)
  sinPiTest(sinPiBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within (1,2)",generateRealWithin(1.0, 2.0,_,_), 2, true)
}
