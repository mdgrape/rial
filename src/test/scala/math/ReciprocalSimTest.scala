//package rial.tests

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


class ReciprocalSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
//class ReciprocalSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ): RealGeneric  = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ): RealGeneric  = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  var counter = 0
  val specialValues = Seq(
      0.0,
     -0.0,
      1.0,
      2.0,
      4.0,
      8.0,
      16.0,
    )
  def generateSpecialValues( spec: RealSpec, r: Random ): RealGeneric = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }

    val x = specialValues(idx)

    val eps = new RealGeneric(spec, 0, spec.exBias - spec.manW, 0).toDouble

    val rnd = r.nextInt(9) - 4 // [-4, 4]
    new RealGeneric(spec, x + (rnd * eps))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def reciprocalTest(t: FuncTableInt, spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric),
    tolerance: Int) = {
    test(s"1/x, format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      var errored = false
      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble
        val z0 = 1.0/x0
        val z0r = new RealGeneric(spec, z0)
        val zi = ReciprocalSim.reciprocalSimGeneric( t, x )
        val zd = zi.toDouble
        val errf = zd-z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toLong
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")

        if (z0r.value != zi.value) {
          // println(f"error: x=${x.toDouble}(${x.sgn}|${x.ex}|${x.man.toLong.toBinaryString}), "+
          //   f"${z0r.toDouble}(${z0r.sgn}|${z0r.ex}|${z0r.man.toLong.toBinaryString}) != "+
          //   f"${zi.toDouble}(${zi.sgn}|${zi.ex}|${zi.man.toLong.toBinaryString})")
        }
        if (z0.isInfinity) {
          assert(zd.isInfinity)
        } else if (zd.isInfinite) {
          assert(z0r.isInfinite)
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {

          // assert(erri.abs<=tolerance)
          if(erri.abs > tolerance) {
            errored = true
            println(f"error: x=${x.toDouble}(${x.sgn}|${x.ex}|${x.man.toLong.toBinaryString}), "+
              f"${z0r.toDouble}(${z0r.sgn}|${z0r.ex}|${z0r.man.toLong.toBinaryString}) != "+
              f"${zi.toDouble}(${zi.sgn}|${zi.ex}|${zi.man.toLong.toBinaryString})")
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
      assert(!errored)
      println( "---------------------------------------------------------------")

    }
  }

  val reciprocalF32TableI = ReciprocalSim.reciprocalTableGeneration( 2, 8, 23, 23+2 )

  reciprocalTest(reciprocalF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  reciprocalTest(reciprocalF32TableI, RealSpec.Float32Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )
  reciprocalTest(reciprocalF32TableI, RealSpec.Float32Spec, n, r,
    "Test Special Values",generateSpecialValues(_,_), 1 )

  val reciprocalBF16TableI = ReciprocalSim.reciprocalTableGeneration( 0, 7, 7, 7 )

  reciprocalTest(reciprocalBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  reciprocalTest(reciprocalBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )
  reciprocalTest(reciprocalBF16TableI, RealSpec.BFloat16Spec, n, r,
    "Test Special Values",generateSpecialValues(_,_), 1 )

  val float48Spec = new RealSpec(10, 511, 37)
  val reciprocalFP48TableI = ReciprocalSim.reciprocalTableGeneration(3, 10, 37, 37+2 )

  reciprocalTest(reciprocalFP48TableI, float48Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  reciprocalTest(reciprocalFP48TableI, float48Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )
  reciprocalTest(reciprocalFP48TableI, float48Spec, n, r,
    "Test Special Values",generateSpecialValues(_,_), 1 )

  val reciprocalFP64TableI = ReciprocalSim.reciprocalTableGeneration(3, 12, 52, 52+2 )

  reciprocalTest(reciprocalFP64TableI, RealSpec.Float64Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 3) // XXX 3!
  reciprocalTest(reciprocalFP64TableI, RealSpec.Float64Spec, n, r,
    "Test All range",generateRealFull(_,_), 3 )
  reciprocalTest(reciprocalFP64TableI, RealSpec.Float64Spec, n, r,
    "Test Special Values",generateSpecialValues(_,_), 3 )
}
