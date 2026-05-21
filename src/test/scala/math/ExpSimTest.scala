package rial.tests

import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.implicits._

import rial.funcs._
import rial.util.ScalaUtil._
import rial.arith._
import rial.approx._

class ExpSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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
      0.0 + (1.0                    ) * pow(2.0, -1022), // for Double
      0.0 + (1.0 + 1 * pow(2.0, -52)) * pow(2.0, -1022), // for Double
      0.0 + (1.0 + 2 * pow(2.0, -52)) * pow(2.0, -1022), // for Double
      0.0 + (1.0 + 3 * pow(2.0, -52)) * pow(2.0, -1022), // for Double
      0.0 + (1.1                    ) * pow(2.0, -1022), // for Double
      0.5,
      1.0,
      2.0,
      log(2.0),
      log(3.0),
      log(4.0),
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) counter = 0
    new RealGeneric(spec, specialValues(idx))
  }

  def errorLSB( x : RealGeneric, y : RealGeneric ) : SafeLong = x.value - y.value

  def expTest(t : PolynomialCoefficients, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"exp(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError    = SafeLong(0)
      var xatMaxError = 0.0
      var zatMaxError = 0.0
      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val z0   = exp(x0)
        val z0r  = new RealGeneric(spec, z0)

        val zi   = ExpSimGeneric(t, x)
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
          if (erri.abs > tolerance) {
            val xsgn = bit(spec.W-1, x.value).toInt
            val xexp = slice(spec.manW, spec.exW, x.value)
            val xman = x.value & maskSL(spec.manW)

            val zsimsgn = bit(spec.W-1, zi.value).toInt
            val zsimexp = slice(spec.manW, spec.exW, zi.value)
            val zsimman = zi.value & maskSL(spec.manW)

            val zrefsgn = bit(spec.W-1, z0r.value).toInt
            val zrefexp = slice(spec.manW, spec.exW, z0r.value)
            val zrefman = z0r.value & maskSL(spec.manW)

            println(f"z0  = ${z0}")
            println(f"test: x   = ${x0}(${x.sgn}|${x.ex}(${x.ex-x.spec.exBias})|${x.man.toLong.toBinaryString})")
            println(f"test: ref = ${z0}(${zrefsgn}|${zrefexp}(${zrefexp-x.spec.exBias})|${zrefman.toLong.toBinaryString})")
            println(f"test: sim = ${zd}(${zsimsgn}|${zsimexp}(${zsimexp-x.spec.exBias})|${zsimman.toLong.toBinaryString})")
            println(f"test: test(${zsimsgn}|${zsimexp}(${zsimexp - spec.exBias})|${zsimman.toLong.toBinaryString}(${zsimman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
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

          assert(erri.abs <= tolerance)

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
              + f"${zatMaxError} != ${exp(xatMaxError)}, "
              + f"diff = ${zatMaxError - exp(xatMaxError)}, x = ${xatMaxError}")
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

  val expFP32Table = mkTable(nOrder = 2, addrW = 8,
                             manW = RealSpec.Float32Spec.manW, extraBits = 3)

  expTest(expFP32Table, RealSpec.Float32Spec, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 2)
  expTest(expFP32Table, RealSpec.Float32Spec, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 2)

  expTest(expFP32Table, RealSpec.Float32Spec, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), 2)
  expTest(expFP32Table, RealSpec.Float32Spec, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), 2)

  expTest(expFP32Table, RealSpec.Float32Spec, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 2)
  expTest(expFP32Table, RealSpec.Float32Spec, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 2)

  expTest(expFP32Table, RealSpec.Float32Spec, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 2)
  expTest(expFP32Table, RealSpec.Float32Spec, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 2)

  expTest(expFP32Table, RealSpec.Float32Spec, n, r,
    "Test Special Values", generateSpecialValues(_,_), 1)

  // --------------------------------------------------------------------------
  // BFloat16

  val expBF16Table = mkTable(nOrder = 0, addrW = 7,
                             manW = RealSpec.BFloat16Spec.manW, extraBits = 1)

  expTest(expBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 1)
  expTest(expBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 1)

  expTest(expBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), 1)
  expTest(expBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), 1)

  expTest(expBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 1)
  expTest(expBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 1)

  expTest(expBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 1)
  expTest(expBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 1)

  expTest(expBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test Special Values", generateSpecialValues(_,_), 2)

  // --------------------------------------------------------------------------
  // FP64

  val expFP64Table = mkTable(nOrder = 3, addrW = 12,
                             manW = RealSpec.Float64Spec.manW, extraBits = 4)

  expTest(expFP64Table, RealSpec.Float64Spec, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), 7)
  expTest(expFP64Table, RealSpec.Float64Spec, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), 7)

  expTest(expFP64Table, RealSpec.Float64Spec, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), 7)
  expTest(expFP64Table, RealSpec.Float64Spec, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), 7)

  expTest(expFP64Table, RealSpec.Float64Spec, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), 7)
  expTest(expFP64Table, RealSpec.Float64Spec, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), 7)

  expTest(expFP64Table, RealSpec.Float64Spec, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), 7)
  expTest(expFP64Table, RealSpec.Float64Spec, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), 7)

  expTest(expFP64Table, RealSpec.Float64Spec, n, r,
    "Test Special Values", generateSpecialValues(_,_), 2)
}
