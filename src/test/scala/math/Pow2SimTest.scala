package rial.tests

import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.implicits._

import rial.math._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

class Pow2SimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRealWithin( from : Double, to : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  var counter = 0
  val specialValues = Seq(
      0.0,
      0.5,
      1.0,
      2.0,
      log(2.0),
      log(3.0),
      log(4.0),
      -0.5,
      -1.0,
      -2.0,
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) counter = 0
    new RealGeneric(spec, specialValues(idx))
  }

  def errorLSB( x : RealGeneric, y : RealGeneric ) : SafeLong = x.value - y.value

  def pow2Test(t : PolynomialCoefficients, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"pow2(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError    = SafeLong(0)
      var xatMaxError = 0.0
      var zatMaxError = 0.0
      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val z0   = pow(2.0, x0)
        val z0r  = new RealGeneric(spec, z0)

        val zi   = Pow2SimGeneric(t, x)
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r).toLong

        if (x0.isInfinity) {
          if(0 < x0) {
            assert(zi.isInfinite)
          } else {
            assert(zi.isZero)
          }
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
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

          assert(erri.abs <= tolerance,
            f"x = ${x0}, sim = ${zd}, ref = ${z0}, error = ${erri}%d LSB")

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
              + f"${zatMaxError} != ${pow(2.0, xatMaxError)}, "
              + f"diff = ${zatMaxError - pow(2.0, xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- ${k}%4d errors (${log2DownL(k)+1}%2d ULPs) positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  private def mkTable(nOrder: Int, addrW: Int, manW: Int, extraBits: Int): PolynomialCoefficients = {
    val t = new ExpTable(nOrder, addrW, manW - addrW, manW + extraBits)
    t.generateTable(t.getBitWidth())
  }

  // --------------------------------------------------------------------------
  // FP32

  val pow2FP32Table = mkTable(nOrder = 2, addrW = 8,
                              manW = RealSpec.Float32Spec.manW, extraBits = 3)

  pow2Test(pow2FP32Table, RealSpec.Float32Spec, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 2)
  pow2Test(pow2FP32Table, RealSpec.Float32Spec, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 2)

  pow2Test(pow2FP32Table, RealSpec.Float32Spec, n, r,
    "Test Large Positive [128, inf]", generateRealWithin(128.0, Double.PositiveInfinity,_,_), 2)
  pow2Test(pow2FP32Table, RealSpec.Float32Spec, n, r,
    "Test Large Negative [-inf, -127]", generateRealWithin(Double.NegativeInfinity, -127.0,_,_), 2)

  pow2Test(pow2FP32Table, RealSpec.Float32Spec, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 2)
  pow2Test(pow2FP32Table, RealSpec.Float32Spec, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 2)

  pow2Test(pow2FP32Table, RealSpec.Float32Spec, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 2)
  pow2Test(pow2FP32Table, RealSpec.Float32Spec, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 2)

  pow2Test(pow2FP32Table, RealSpec.Float32Spec, n, r,
    "Test Special Values", generateSpecialValues(_,_), 1)

  // --------------------------------------------------------------------------
  // BFloat16

  val pow2BF16Table = mkTable(nOrder = 0, addrW = 7,
                              manW = RealSpec.BFloat16Spec.manW, extraBits = 1)

  pow2Test(pow2BF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 1)
  pow2Test(pow2BF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 1)

  pow2Test(pow2BF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Large Positive [128, inf]", generateRealWithin(128.0, Double.PositiveInfinity,_,_), 1)
  pow2Test(pow2BF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Large Negative [-inf, -127]", generateRealWithin(Double.NegativeInfinity, -127.0,_,_), 1)

  pow2Test(pow2BF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 1)
  pow2Test(pow2BF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 1)

  pow2Test(pow2BF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 1)
  pow2Test(pow2BF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 1)

  pow2Test(pow2BF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Special Values", generateSpecialValues(_,_), 1)
}
