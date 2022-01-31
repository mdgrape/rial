
//package rial.tests

//import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

//import org.scalatest.flatspec.AnyFlatSpec
//import org.scalatest.matchers.should.Matchers
//import scopt.OptionParser

import scala.collection.mutable.Map
import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

class InvSqrtSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
//class InvSqrtSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 100000

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

  def invsqrtTest(t : FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric), tolerance : Int ) = {
    test(s"invsqrt(x), format ${spec.toStringShort}, ${generatorStr}") {
      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0
      var errorCount: Map[Int, Int] = Map()

      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = 1.0 / math.sqrt(x0)
        val z0r  = new RealGeneric(spec, z0)
        val zi   = InvSqrtSim.invsqrtSimGeneric( t, x )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble)
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

//           println(f"gen = ${x0}")
//           println(f"ref = ${z0}, 1/ref^2 = ${1.0 / (z0*z0)}")
//           println(f"sim = ${zd}, 1/sim^2 = ${1.0 / (zd*zd)}")
//           println(f"test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman}(${zrefman.toLong}%x)) = invsqrt(x = ${xsgn}|${xexp}(${xexp - spec.exBias})|${xman}(${xman.toLong}%x))")
        }

        if (x.isNaN) {
          assert(zi.isNaN)
        } else if (x.isInfinite) {
          assert(zi.isZero)
        } else if (x.sgn == 1 && !x.isZero) {
          assert(zi.isInfinite)
        } else if (x.isZero) {
          assert(zi.isInfinite)
        } else {
          assert(erri.abs <= tolerance)

          if (erri.toInt != 0) {
            if (!errorCount.contains(erri.toInt)) {
              errorCount += (erri.toInt -> 0)
            }
            errorCount(erri.toInt) += 1
          }
          if (maxError < erri.abs) {
            maxError = erri.abs
            xatMaxError = x0
            zatMaxError = zd
          }
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      for (kv <- errorCount.toSeq.sortBy(_._1)) {
        val (level, count) = kv
        println(f"${level}%2d errors happen ${count}%6d / ${n} times (${count.toDouble / n})")
      }
      println(f"N=$n%d : the largest error is ${maxError.toInt}%d where the value is ${zatMaxError} != ${1.0 / math.sqrt(xatMaxError)} diff = ${zatMaxError - 1.0 / math.sqrt(xatMaxError)}")
    }
  }

  val invsqrtF32TableI  = InvSqrtSim.invsqrtTableGeneration( 2, 8, 23, 23+2 )

  invsqrtTest(invsqrtF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  invsqrtTest(invsqrtF32TableI, RealSpec.Float32Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )

  val invsqrtBF16TableI  = InvSqrtSim.invsqrtTableGeneration( 0, 7, 7, 7 )

  invsqrtTest(invsqrtBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  invsqrtTest(invsqrtBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )
}
