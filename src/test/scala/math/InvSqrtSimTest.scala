
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

class InvSqrtSimTest extends FunSuite with BeforeAndAfterAllConfigMap {
//class InvSqrtSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def invsqrtTest(t_even: FuncTableInt, t_odd: FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    test(s"invsqrt(x), format ${spec.toStringShort}, ${generatorStr}") {
      var err1lsbPos = 0
      var err1lsbNeg = 0
      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = 1.0 / math.sqrt(x0)
        val z0r  = new RealGeneric(spec, z0)
        val zi   = InvSqrtSim.invsqrtSimGeneric( t_even, t_odd, x )
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

//           val zref = new RealGeneric(spec, zrefsgn, zrefexp.toInt, zrefman.toInt)
//           val zrefd = zref.toDouble

          println(f"gen = ${x0}")
          println(f"ref = ${z0}, 1/ref^2 = ${1.0 / (z0*z0)}")
          println(f"sim = ${zd}, 1/sim^2 = ${1.0 / (zd*zd)}")
          println(f"test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman}(${zrefman.toLong}%x)) = invsqrt(x = ${xsgn}|${xexp}(${xexp - spec.exBias})|${xman}(${xman.toLong}%x))")
        }
        if (z0.isInfinity) {
          assert(zd.isInfinity)
        } else if (zd.isInfinite) {
          assert(z0r.isInfinite)
        } else if (x.isNaN) {
          assert(zi.isNaN)
        } else if (x.sgn == 1 && !x.isZero) {
          assert(zi.isNaN)
        } else {
          if (erri.abs>=2.0) {
            println(f"Error more than 2 LSB : ${x.toDouble}%14.7e : $z0%14.7e ${zi.toDouble}%14.7e $errf%14.7e $erri%f")
          } else if (erri>=1.0) err1lsbPos+=1
          else if (erri<= -1.0) err1lsbNeg+=1
          assert(erri.abs<=1)
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
    }
  }

  val invsqrtF32TableIEven = InvSqrtSim.invsqrtTableGenerationEven( 2, 8, 23, 23+2 )
  val invsqrtF32TableIOdd  = InvSqrtSim.invsqrtTableGenerationOdd ( 2, 8, 23, 23+2 )

  invsqrtTest(invsqrtF32TableIEven, invsqrtF32TableIOdd, RealSpec.Float32Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  invsqrtTest(invsqrtF32TableIEven, invsqrtF32TableIOdd, RealSpec.Float32Spec, n, r,
    "Test All range",generateRealFull(_,_) )

  val invsqrtBF16TableIEven = InvSqrtSim.invsqrtTableGenerationEven( 0, 7, 7, 7 )
  val invsqrtBF16TableIOdd  = InvSqrtSim.invsqrtTableGenerationOdd ( 0, 7, 7, 7 )

  invsqrtTest(invsqrtBF16TableIEven, invsqrtBF16TableIOdd, RealSpec.BFloat16Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  invsqrtTest(invsqrtBF16TableIEven, invsqrtBF16TableIOdd, RealSpec.BFloat16Spec, n, r,
    "Test All range",generateRealFull(_,_) )
}
