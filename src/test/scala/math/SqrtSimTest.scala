
//package rial.tests


import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}



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

class SqrtSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
//class SqrtSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateReal1to4( spec: RealSpec, r : Random) = {
    val rD : Double = r.nextDouble()*3.0+1.0
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

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

  def sqrtTest(t: FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric),
    tolerance: Int) = {
    test(s"sqrt(x), format ${spec.toStringShort}, ${generatorStr}") {

      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = if(x0 < 0) {0.0} else {math.sqrt(x0)}
        val z0r  = new RealGeneric(spec, z0)
        val zi   = SqrtSim.sqrtSimGeneric( t, x )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toLong
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if ((z0r.value - zi.value).abs > tolerance) {
          val xsgn = bit(spec.W-1, x.value).toInt
          val xexp = slice(spec.manW, spec.exW, x.value)
          val xman = x.value & maskSL(spec.manW)

          val ztestsgn = bit(spec.W-1, zi.value).toInt
          val ztestexp = slice(spec.manW, spec.exW, zi.value)
          val ztestman = zi.value & maskSL(spec.manW)

          val zrefsgn = bit(spec.W-1, z0r.value).toInt
          val zrefexp = slice(spec.manW, spec.exW, z0r.value)
          val zrefman = z0r.value & maskSL(spec.manW)

          println(f"x0  = ${x0} (${x.sgn}|${x.ex}|${x.man.toLong.toBinaryString}|)")
          println(f"ref = ${z0}, ref^2 = ${z0*z0}")
          println(f"sim = ${zd}, sim^2 = ${zd*zd}")
          println(f"test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman}(${zrefman.toLong}%x)) = sqrt(x = ${xsgn}|${xexp}(${xexp - spec.exBias})|${xman}(${xman.toLong}%x))")
        }
        if (z0.isInfinity) {
          assert(zd.isInfinity)
        } else if (zd.isInfinite) {
          assert(z0r.isInfinite)
        } else if (x.isNaN) {
          assert(zi.isNaN)
        } else if (x.sgn == 1 && !x.isZero) {
          assert(zi.isZero)
        } else {
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
              + f"${zatMaxError} != ${sqrt(xatMaxError)}, "
              + f"diff = ${zatMaxError - sqrt(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- ${k}%4d errors (${log2DownL(k)+1}%2d ULPs) positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val sqrtF32TableI = SqrtSim.sqrtTableGeneration( 2, 8, 23, 23+2 )

  sqrtTest(sqrtF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [0, 4)",generateReal1to4(_,_), 1)
  sqrtTest(sqrtF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  sqrtTest(sqrtF32TableI, RealSpec.Float32Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )

  val sqrtBF16TableI = SqrtSim.sqrtTableGeneration(0, 7, 7, 7 ) // [1,2) + [2,4) + 1.0

  sqrtTest(sqrtBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  sqrtTest(sqrtBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )

  val float48Spec = new RealSpec(10, 511, 37)
  val sqrtFP48TableI = SqrtSim.sqrtTableGeneration(3, 10, 37, 37+4 )

  sqrtTest(sqrtFP48TableI, float48Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  sqrtTest(sqrtFP48TableI, float48Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )

  val sqrtFP64TableI = SqrtSim.sqrtTableGeneration(3, 12, 52, 52+4 )

  sqrtTest(sqrtFP64TableI, RealSpec.Float64Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 3) // XXX 3!
  sqrtTest(sqrtFP64TableI, RealSpec.Float64Spec, n, r,
    "Test All range",generateRealFull(_,_), 3 )
}
