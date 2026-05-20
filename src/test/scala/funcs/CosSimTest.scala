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

class CosSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
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
      +0.0,                -0.0,
      math.Pi / 4.0,       -math.Pi / 4.0,
      math.Pi / 3.0,       -math.Pi / 3.0,
      math.Pi / 2.0,       -math.Pi / 2.0,
      math.Pi * 2.0 / 3.0, -math.Pi * 2.0 / 3.0,
      math.Pi * 3.0 / 4.0, -math.Pi * 3.0 / 4.0,
      math.Pi,             -math.Pi,
      math.Pi * 5.0 / 4.0, -math.Pi * 5.0 / 4.0,
      math.Pi * 4.0 / 3.0, -math.Pi * 4.0 / 3.0,
      math.Pi * 3.0 / 2.0, -math.Pi * 3.0 / 2.0,
      math.Pi * 5.0 / 3.0, -math.Pi * 5.0 / 3.0,
      math.Pi * 7.0 / 4.0, -math.Pi * 7.0 / 4.0,
      math.Pi * 2.0,       -math.Pi * 2.0
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }

  def errorLSB( x : RealGeneric, y : RealGeneric ) : SafeLong = x.value - y.value

  def cosTest( t: PolynomialCoefficients, spec : RealSpec, n : Int, r : Random,
    generatorStr    : String,
    generator       : ( (RealSpec, Random) => RealGeneric),
    tolerance       : Int ) = {
    test(s"cos(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError    = SafeLong(0)
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map[Int, (Int, Int)]()

      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble

        val z0   = cos(x0)
        val z0r  = new RealGeneric(spec, z0)

        val zi   = SinCosSimGeneric(false, t, x)
        val zd   = zi.toDouble
        val erri = errorLSB(zi, z0r).toInt

        if (x0.isInfinity) {
          if(0 < x0) {
            assert(zi.isInfinite)
          } else {
            assert(zi.isZero)
          }
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if(erri.abs > tolerance) {
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
            println(f"z0r = ${z0r.sgn}| ${z0r.ex}| ${z0r.man.toLong.toBinaryString}")
            println(f"zi  = ${zi.sgn }| ${zi .ex}| ${zi .man.toLong.toBinaryString}")
            println(f"zd  = ${zd}")
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
              + f"${zatMaxError} != ${cos(xatMaxError)}, "
              + f"diff = ${zatMaxError - cos(xatMaxError)}, x = ${xatMaxError}")
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

  val cosFP32Table = mkTable(
    nOrder    = 2,
    addrW     = 8,
    manW      = RealSpec.Float32Spec.manW,
    extraBits = 3
  )

  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 3)

  // beyond ±2pi
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [-10pi, -2pi]", generateRealWithin(-10 * Pi, -2 * Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [2pi, 10pi]", generateRealWithin(2.0 * Pi, 10.0 * Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [16pi, 32pi]", generateRealWithin(16.0 * Pi, 32.0 * Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [32pi, 64pi]", generateRealWithin(32.0 * Pi, 64.0 * Pi,_,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [-32pi, -16pi]", generateRealWithin( -32.0 * Pi,-16.0 * Pi, _,_), 3)
  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos Within [-64pi, -32pi]", generateRealWithin( -64.0 * Pi,-32.0 * Pi, _,_), 3)

  cosTest(cosFP32Table, RealSpec.Float32Spec, n, r,
    "Test cos with special values", generateSpecialValues(_,_), 2)

  // --------------------------------------------------------------------------
  // BFloat16

  val cosBF16Table = mkTable(
    nOrder    = 0,
    addrW     = 7,
    manW      = RealSpec.BFloat16Spec.manW,
    extraBits = 1
  )

  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 3)

  // beyond ±2pi
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [-10pi, -2pi]", generateRealWithin(-10 * Pi, -2 * Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [2pi, 10pi]", generateRealWithin(2.0 * Pi, 10.0 * Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [16pi, 32pi]", generateRealWithin(16.0 * Pi, 32.0 * Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [32pi, 64pi]", generateRealWithin(32.0 * Pi, 64.0 * Pi,_,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [-32pi, -16pi]", generateRealWithin( -32.0 * Pi,-16.0 * Pi, _,_), 3)
  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos Within [-64pi, -32pi]", generateRealWithin( -64.0 * Pi,-32.0 * Pi, _,_), 3)

  cosTest(cosBF16Table, RealSpec.BFloat16Spec, n, r,
    "Test cos with special values", generateSpecialValues(_,_), 1)

  // --------------------------------------------------------------------------
  // a hypothetical Float48 (10-bit exponent, 37-bit mantissa)

  val float48Spec = new RealSpec(10, 511, 37)
  val cosFP48Table = mkTable(
    nOrder    = 3,
    addrW     = 10,
    manW      = float48Spec.manW,
    extraBits = 4
  )

  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 4)
  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 4)
  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 4)
  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 4)
  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 4)
  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 4)
  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 4)
  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 4)
  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 4)
  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 4)

  cosTest(cosFP48Table, float48Spec, n, r,
    "Test cos with special values", generateSpecialValues(_,_), 1)

  // --------------------------------------------------------------------------
  // FP64

  val cosFP64Table = mkTable(
    nOrder    = 3,
    addrW     = 12,
    manW      = RealSpec.Float64Spec.manW,
    extraBits = 4
  )

  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [-2pi, -1.5pi]", generateRealWithin(-2 * Pi, -1.5 * Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [-1.5pi, -pi]", generateRealWithin(-1.5 * Pi, -Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [-pi, -pi/2]", generateRealWithin(-Pi, -0.5 * Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [-pi/2, -2^-12pi]", generateRealWithin(-0.5*Pi, -pow(2.0, -12)*Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [-2^-12pi, 0]", generateRealWithin(-pow(2.0, -12)*Pi, 0.0,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [0, 2^-12pi]", generateRealWithin(0.0, pow(2.0, -12)*Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [2^-12pi, pi/2]", generateRealWithin(pow(2.0, -12)*Pi, 0.5*Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [pi/2, pi]", generateRealWithin(0.5*Pi, Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [pi, 3/2pi]", generateRealWithin(Pi, 1.5*Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos Within [3/2pi, 2pi]", generateRealWithin(1.5*Pi, 2.0*Pi,_,_), 7)

  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos with special values", generateSpecialValues(_,_), 2)

  val delta = pow(2.0, -32)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos around 0", generateRealWithin(-delta*Pi, delta*Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos around pi/2", generateRealWithin((0.5-delta)*Pi, (0.5+delta)*Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos around pi", generateRealWithin((1-delta)*Pi, (1+delta)*Pi,_,_), 7)
  cosTest(cosFP64Table, RealSpec.Float64Spec, n, r,
    "Test cos around 3pi/2", generateRealWithin((1.5-delta)*Pi, (1.5+delta)*Pi,_,_), 7)
}
