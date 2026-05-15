package rial.tests.funcs

import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.funcs._
import rial.util.ScalaUtil._
import rial.arith._
import rial.approx._

class ReciprocalSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
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

  var counter = 0
  val specialValues = Seq(
      0.0,
     -0.0,
      1.0,
      2.0,
      4.0,
     -1.0,
     -2.0,
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def reciprocalTest(t: PolynomialCoefficients, spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric),
    tolerance: Int) = {
    test(s"reciprocal(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = if(x0 == 0.0) {Double.PositiveInfinity} else {1.0/x0}
        val z0r  = new RealGeneric(spec, z0)
        val zi   = ReciprocalSimGeneric( t, x )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toLong

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
          println(f"ref = ${z0}, 1/ref = ${1.0/z0}")
          println(f"sim = ${zd}, 1/sim = ${1.0/zd}")
          println(f"test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman}(${zrefman.toLong}%x)) = 1/(x = ${xsgn}|${xexp}(${xexp - spec.exBias})|${xman}(${xman.toLong}%x))")
        }
        if (z0.isInfinity) {
          assert(zd.isInfinity)
        } else if (zd.isInfinite) {
          assert(z0r.isInfinite)
        } else if (x.isNaN) {
          assert(zi.isNaN)
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
              + f"${zatMaxError} != ${1.0/xatMaxError}, "
              + f"diff = ${zatMaxError - 1.0/xatMaxError}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- ${k}%4d errors (${log2DownL(k)+1}%2d ULPs) positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val reciprocalF32Table = ReciprocalTable.fp32()
  val reciprocalF32Coeff = reciprocalF32Table.generateTable(reciprocalF32Table.getBitWidth())

  reciprocalTest(reciprocalF32Coeff, RealSpec.Float32Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  reciprocalTest(reciprocalF32Coeff, RealSpec.Float32Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )
  reciprocalTest(reciprocalF32Coeff, RealSpec.Float32Spec, n, r,
    "Test special value",generateSpecialValues(_,_), 1 )

  val reciprocalBF16Table = ReciprocalTable.bf16()
  val reciprocalBF16Coeff = reciprocalBF16Table.generateTable(reciprocalBF16Table.getBitWidth())

  reciprocalTest(reciprocalBF16Coeff, RealSpec.BFloat16Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
  reciprocalTest(reciprocalBF16Coeff, RealSpec.BFloat16Spec, n, r,
    "Test All range",generateRealFull(_,_), 1 )
  reciprocalTest(reciprocalBF16Coeff, RealSpec.BFloat16Spec, n, r,
    "Test special value",generateSpecialValues(_,_), 1 )
}
