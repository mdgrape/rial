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

class SinPiSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  /** sin(pi x), with the argument reduced so that an integral x gives zero. */
  def sinPi(x: Double): Double = {
    val n = rint(x)
    val s = if (rint(n / 2.0) * 2.0 == n) { 1.0 } else { -1.0 }
    s * sin(Pi * (x - n))
  }

  def generateRealWithin( from : Double, to : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  var counter = 0
  val specialValues = Seq(
      +0.0,      -0.0,
      0.25,      -0.25,
      1.0 / 3.0, -1.0 / 3.0,
      0.5,       -0.5,
      2.0 / 3.0, -2.0 / 3.0,
      0.75,      -0.75,
      1.0,       -1.0,
      1.25,      -1.25,
      4.0 / 3.0, -4.0 / 3.0,
      1.5,       -1.5,
      5.0 / 3.0, -5.0 / 3.0,
      1.75,      -1.75,
      2.0,       -2.0
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }

  // integer bit-pattern difference, in LSB units of the result format.
  def errorLSB( x : RealGeneric, y : RealGeneric ) : SafeLong = x.value - y.value

  def sinPiTest( t: PolynomialCoefficients, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int
    ) = {
    test(s"sinPi(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError    = SafeLong(0)
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val z0   = sinPi(x0)
        val z0r  = new RealGeneric(spec, z0)

        val zi   = SinPiSimGeneric(true, t, x)
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r).toInt

        if (x.ex == 0) {
          assert(zi.ex == 0 && x.sgn == zi.sgn)
        } else if (zi.isZero && z0r.isZero) {
          // an integral x gives an exact zero, whose sign the fold of the
          // quadrant does not carry
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
              + f"${zatMaxError} != ${sinPi(xatMaxError)}, "
              + f"diff = ${zatMaxError - sinPi(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- $k errors positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  private def mkTable(nOrder: Int, addrW: Int, manW: Int, extraBits: Int): PolynomialCoefficients = {
    val t = new SinCosTable(nOrder, addrW, manW - addrW, manW + extraBits)
    t.generateTable(t.getBitWidth())
  }

  // --------------------------------------------------------------------------
  // FP32

  val sinPiFP32Table = mkTable(
    nOrder    = 2,
    addrW     = 8,
    manW      = RealSpec.Float32Spec.manW,
    extraBits = 3
  )

  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [-2, -1.5]", generateRealWithin(-2.0, -1.5,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [-1.5, -1]", generateRealWithin(-1.5, -1.0,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [-1, -0.5]", generateRealWithin(-1.0, -0.5,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [-0.5, -2^-12]", generateRealWithin(-0.5, -pow(2.0, -12),_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [-2^-12, 0]", generateRealWithin(-pow(2.0, -12), 0.0,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [0, 2^-12]", generateRealWithin(0.0, pow(2.0, -12),_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [2^-12, 0.5]", generateRealWithin(pow(2.0, -12), 0.5,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [0.5, 1]", generateRealWithin(0.5, 1.0,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [1, 1.5]", generateRealWithin(1.0, 1.5,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [1.5, 2]", generateRealWithin(1.5, 2.0,_,_), 3)

  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi with special values", generateSpecialValues(_,_), 1)

  // beyond +-2
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [-10, -2]", generateRealWithin(-10.0, -2.0,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [2, 10]", generateRealWithin(2.0, 10.0,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [16, 32]", generateRealWithin(16.0, 32.0,_,_), 3)
  sinPiTest(sinPiFP32Table, RealSpec.Float32Spec, n, r,
    "Test sinPi Within [-64, -32]", generateRealWithin(-64.0, -32.0,_,_), 3)

  // --------------------------------------------------------------------------
  // BFloat16

  val sinPiBF16Table = mkTable(
    nOrder    = 0,
    addrW     = 7,
    manW      = RealSpec.BFloat16Spec.manW,
    extraBits = 1
  )

  sinPiTest(sinPiBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test sinPi Within [-2, -1.5]", generateRealWithin(-2.0, -1.5,_,_), 3)
  sinPiTest(sinPiBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test sinPi Within [-1, -0.5]", generateRealWithin(-1.0, -0.5,_,_), 3)
  sinPiTest(sinPiBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test sinPi Within [-0.5, -2^-12]", generateRealWithin(-0.5, -pow(2.0, -12),_,_), 3)
  sinPiTest(sinPiBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test sinPi Within [2^-12, 0.5]", generateRealWithin(pow(2.0, -12), 0.5,_,_), 3)
  sinPiTest(sinPiBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test sinPi Within [0.5, 1]", generateRealWithin(0.5, 1.0,_,_), 3)
  sinPiTest(sinPiBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test sinPi Within [1.5, 2]", generateRealWithin(1.5, 2.0,_,_), 3)

  sinPiTest(sinPiBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test sinPi with special values", generateSpecialValues(_,_), 1)

  sinPiTest(sinPiBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test sinPi Within [2, 10]", generateRealWithin(2.0, 10.0,_,_), 3)
  sinPiTest(sinPiBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test sinPi Within [16, 32]", generateRealWithin(16.0, 32.0,_,_), 3)
}
