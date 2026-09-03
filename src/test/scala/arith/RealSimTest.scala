package rial.tests


import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}



//import scopt.OptionParser

import scala.util.Random
import scala.math._
import rial.arith._
import rial.util.ScalaUtil._

import spire.math.SafeLong

class RealSimTest extends AnyFunSuite with Matchers with BeforeAndAfterAllConfigMap {
//class RealSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

  test("Double set->get test") {
    for(i <- 1 to n) {
      //val x0 = (r.nextDouble()-0.5)*128.0
      val x0 = 1.0
      val xr = RealGeneric.fromDouble(RealSpec.Float64Spec, x0)
      val xd = xr.toDouble
      xd should be (x0)
    }
  }

  test("Double mul test") {
    for(i <- 1 to n) {
      val x = (r.nextDouble()-0.5)*128.0
      val y = (r.nextDouble()-0.5)*128.0
      val z0 = x * y
      val xr = RealGeneric.fromDouble(RealSpec.Float64Spec, x)
      val yr = RealGeneric.fromDouble(RealSpec.Float64Spec, y)
      val zr = xr.multiply(RealSpec.Float64Spec, RoundSpec.roundToEven,yr)
      val zd = zr.toDouble
      zd should be (z0)
    }
  }

  test("Double add test") {
    for(i <- 1 to n) {
      val x = (r.nextDouble()-0.5)*128.0
      val y = (r.nextDouble()-0.5)*128.0
      val z0 = x + y
      val xr = RealGeneric.fromDouble(RealSpec.Float64Spec, x)
      val yr = RealGeneric.fromDouble(RealSpec.Float64Spec, y)
      val zr = xr.add(RealSpec.Float64Spec, RoundSpec.roundToEven,yr)
      val zd = zr.toDouble
      zd should be (z0)
    }
  }

  test("Double negate test") {
    for(i <- 1 to n) {
      val x = (r.nextDouble()-0.5)*128.0
      val z0 = -x
      val xr = RealGeneric.fromDouble(RealSpec.Float64Spec, x)
      val zr = xr.negate()
      val zd = zr.toDouble
      zd should be (z0)
    }
  }

  test("Float fmadd test") {
    val tolerance = 1
    for(i <- 1 to n) {
      val x = (r.nextFloat()-0.5f)*128.0f
      val y = (r.nextFloat()-0.5f)*128.0f
      val z = (r.nextFloat()-0.5f)*128.0f
      val w0 = RealGeneric.fromFloat(RealSpec.Float32Spec, java.lang.Math.fma(x, y, z))
      val xr = RealGeneric.fromFloat(RealSpec.Float32Spec, x)
      val yr = RealGeneric.fromFloat(RealSpec.Float32Spec, y)
      val zr = RealGeneric.fromFloat(RealSpec.Float32Spec, z)
      val wr = xr.fmadd(RealSpec.Float32Spec, RoundSpec.roundToEven, yr, zr)

      val diff = (w0.value - wr.value).toInt.abs
      if(diff > tolerance) {
        println(f"x  = ${xr.sgn}|${xr.ex}|${xr.man.toLong.toBinaryString}(${xr.toDouble})")
        println(f"y  = ${yr.sgn}|${yr.ex}|${yr.man.toLong.toBinaryString}(${yr.toDouble})")
        println(f"z  = ${zr.sgn}|${zr.ex}|${zr.man.toLong.toBinaryString}(${zr.toDouble})")
        println(f"w0 = ${w0.sgn}|${w0.ex}|${w0.man.toLong.toBinaryString}(${w0.toDouble})")
        println(f"wr = ${wr.sgn}|${wr.ex}|${wr.man.toLong.toBinaryString}(${wr.toDouble})")
      }

      diff should be <= tolerance
    }
  }
  test("Double fmadd test") {
    val tolerance = 1
    for(i <- 1 to n) {
      val x = (r.nextDouble()-0.5)*128.0
      val y = (r.nextDouble()-0.5)*128.0
      val z = (r.nextDouble()-0.5)*128.0
      val w0 = RealGeneric.fromDouble(RealSpec.Float64Spec, java.lang.Math.fma(x, y, z))
      val xr = RealGeneric.fromDouble(RealSpec.Float64Spec, x)
      val yr = RealGeneric.fromDouble(RealSpec.Float64Spec, y)
      val zr = RealGeneric.fromDouble(RealSpec.Float64Spec, z)
      val wr = xr.fmadd(RealSpec.Float64Spec, RoundSpec.roundToEven, yr, zr)

      val diff = (w0.value - wr.value).toInt.abs
      if(diff > tolerance) {
        println(f"x  = ${xr.sgn}|${xr.ex}|${xr.man.toLong.toBinaryString}(${xr.toDouble})")
        println(f"y  = ${yr.sgn}|${yr.ex}|${yr.man.toLong.toBinaryString}(${yr.toDouble})")
        println(f"z  = ${zr.sgn}|${zr.ex}|${zr.man.toLong.toBinaryString}(${zr.toDouble})")
        println(f"w0 = ${w0.sgn}|${w0.ex}|${w0.man.toLong.toBinaryString}(${w0.toDouble})")
        println(f"wr = ${wr.sgn}|${wr.ex}|${wr.man.toLong.toBinaryString}(${wr.toDouble})")
      }

      diff should be <= tolerance
    }
  }

  test("Float 3-op add test") {
    val tolerance = 3
    for(i <- 1 to n) {
      val x = (r.nextFloat()-0.5f)*128.0f
      val y = (r.nextFloat()-0.5f)*128.0f
      val z = (r.nextFloat()-0.5f)*128.0f
      val w0 = RealGeneric.fromFloat(RealSpec.Float32Spec, x + y + z)
      val xr = RealGeneric.fromFloat(RealSpec.Float32Spec, x)
      val yr = RealGeneric.fromFloat(RealSpec.Float32Spec, y)
      val zr = RealGeneric.fromFloat(RealSpec.Float32Spec, z)
      val wr = xr.add3op(RealSpec.Float32Spec, RoundSpec.roundToEven, yr, zr)

      val diff = (w0.value - wr.value).toInt.abs
      if(diff > tolerance) {
        println(f"x  = ${xr.sgn}|${xr.ex}|${xr.man.toLong.toBinaryString}(${xr.toDouble})")
        println(f"y  = ${yr.sgn}|${yr.ex}|${yr.man.toLong.toBinaryString}(${yr.toDouble})")
        println(f"z  = ${zr.sgn}|${zr.ex}|${zr.man.toLong.toBinaryString}(${zr.toDouble})")
        println(f"w0 = ${w0.sgn}|${w0.ex}|${w0.man.toLong.toBinaryString}(${w0.toDouble})")
        println(f"wr = ${wr.sgn}|${wr.ex}|${wr.man.toLong.toBinaryString}(${wr.toDouble})")
      }
      diff should be <= tolerance
    }
  }
  test("Double 3-op add test") {
    val tolerance = 3
    for(i <- 1 to n) {
      val x = (r.nextDouble()-0.5)*128.0
      val y = (r.nextDouble()-0.5)*128.0
      val z = (r.nextDouble()-0.5)*128.0
      val w0 = RealGeneric.fromDouble(RealSpec.Float64Spec, x + y + z)
      val xr = RealGeneric.fromDouble(RealSpec.Float64Spec, x)
      val yr = RealGeneric.fromDouble(RealSpec.Float64Spec, y)
      val zr = RealGeneric.fromDouble(RealSpec.Float64Spec, z)
      val wr = xr.add3op(RealSpec.Float64Spec, RoundSpec.roundToEven, yr, zr)

      val diff = (w0.value - wr.value).toInt.abs
      if(diff > tolerance) {
        println(f"x  = ${xr.sgn}|${xr.ex}|${xr.man.toLong.toBinaryString}(${xr.toDouble})")
        println(f"y  = ${yr.sgn}|${yr.ex}|${yr.man.toLong.toBinaryString}(${yr.toDouble})")
        println(f"z  = ${zr.sgn}|${zr.ex}|${zr.man.toLong.toBinaryString}(${zr.toDouble})")
        println(f"w0 = ${w0.sgn}|${w0.ex}|${w0.man.toLong.toBinaryString}(${w0.toDouble})")
        println(f"wr = ${wr.sgn}|${wr.ex}|${wr.man.toLong.toBinaryString}(${wr.toDouble})")
      }
      diff should be <= tolerance
    }
  }

  //--------------------------------------------------------------------------
  // Sign of zero and conversion into the output spec.
  //
  // IEEE 754 (2019) 6.3 : an exact zero result is +0 (unless both addends are
  // -0), while a non-zero result that underflows to zero keeps the sign of the
  // exact result. `x + 0` is `x` rounded into the output spec, so a narrower
  // output has to round the mantissa and to check the output exponent range.
  // The expected values below are the ones AddFPGeneric produces.

  val allRoundSpecs = Seq(RoundSpec.round, RoundSpec.roundToEven, RoundSpec.truncate,
    RoundSpec.truncateToZero, RoundSpec.ceil, RoundSpec.ceilToInfinite)

  val f64  = RealSpec.Float64Spec
  val f32  = RealSpec.Float32Spec
  val f16  = RealSpec.Float16Spec
  val bf16 = RealSpec.BFloat16Spec

  test("add: an exact zero sum is +0 in every rounding mode") {
    for(roundSpec <- allRoundSpecs) {
      val x = new RealGeneric(f64, 0, f64.exBias, SafeLong(0x123456789ABCL))
      val z = x.add(f64, roundSpec, x.negate())
      withClue(s"roundSpec=$roundSpec: ") {
        z.value should be (SafeLong(0))
      }
    }
  }

  test("add: a non-zero sum that underflows keeps its sign") {
    // -(1+2^-52) * 2^-1022 + 2^-1022 = -2^-1074, which flushes to -0
    val xNeg = new RealGeneric(f64, 1, 1, SafeLong(1))
    val yPos = new RealGeneric(f64, 0, 1, SafeLong(0))
    for(roundSpec <- allRoundSpecs) {
      withClue(s"roundSpec=$roundSpec: ") {
        val zNeg = xNeg.add(f64, roundSpec, yPos)
        zNeg.isZero should be (true)
        zNeg.sgn   should be (1)

        val zPos = xNeg.negate().add(f64, roundSpec, yPos.negate())
        zPos.isZero should be (true)
        zPos.sgn    should be (0)
      }
    }
  }

  test("fmadd: an exact zero result is +0, an underflowed one keeps its sign") {
    val one = new RealGeneric(f64, 0, f64.exBias,   SafeLong(0))
    val two = new RealGeneric(f64, 0, f64.exBias+1, SafeLong(0))
    one.fmadd(f64, RoundSpec.roundToEven, two, two.negate()).value should be (SafeLong(0))

    // (1+2^-52) * -2^-1022 + 2^-1022 = -2^-1074, which flushes to -0
    val x = new RealGeneric(f64, 0, f64.exBias, SafeLong(1))
    val y = new RealGeneric(f64, 1, 1, SafeLong(0))
    val z = new RealGeneric(f64, 0, 1, SafeLong(0))

    val wNeg = x.fmadd(f64, RoundSpec.roundToEven, y, z)
    wNeg.isZero should be (true)
    wNeg.sgn    should be (1)

    val wPos = x.fmadd(f64, RoundSpec.roundToEven, y.negate(), z.negate())
    wPos.isZero should be (true)
    wPos.sgn    should be (0)
  }

  test("fmadd: a zero product returns the addend in the output spec") {
    val zero32 = RealGeneric.zero(f32)
    val one32  = new RealGeneric(f32, 0, f32.exBias, SafeLong(0))
    // 0 * 1 + 1 = 1
    zero32.fmadd(f32, RoundSpec.roundToEven, one32, one32).value should be (SafeLong(0x3F800000L))
    // the addend is rounded into the narrower output spec
    val zRound = new RealGeneric(f32, 0, f32.exBias, (SafeLong(1) << 15) + 1)
    zero32.fmadd(bf16, RoundSpec.roundToEven, one32, zRound).value should be (SafeLong(0x3F81))
    // and it underflows to a signed zero
    val zTiny = RealGeneric.fromDouble(f64, -1.0e-300)
    zero32.fmadd(f32, RoundSpec.roundToEven, one32, zTiny).value should be (SafeLong(0x80000000L))
    // 0 * 0 + 0 : -0 only when both are -0
    val nzero32 = RealGeneric.zero(f32, 1)
    zero32.fmadd(f32, RoundSpec.roundToEven, zero32, nzero32).value should be (SafeLong(0))
    nzero32.fmadd(f32, RoundSpec.roundToEven, zero32, nzero32).value should be (SafeLong(0x80000000L))
  }

  test("add3op: an exact zero sum is +0, an underflowed one keeps its sign") {
    val one = new RealGeneric(f64, 0, f64.exBias,   SafeLong(0))
    val two = new RealGeneric(f64, 0, f64.exBias+1, SafeLong(0))
    one.add3op(f64, RoundSpec.roundToEven, one, two.negate()).value should be (SafeLong(0))

    // -(1+2^-52) * 2^-1022 + 2^-1022 + 0 = -2^-1074, which flushes to -0
    val xNeg  = new RealGeneric(f64, 1, 1, SafeLong(1))
    val yPos  = new RealGeneric(f64, 0, 1, SafeLong(0))
    val zZero = RealGeneric.zero(f64)

    val wNeg = xNeg.add3op(f64, RoundSpec.roundToEven, yPos, zZero)
    wNeg.isZero should be (true)
    wNeg.sgn    should be (1)

    val wPos = xNeg.negate().add3op(f64, RoundSpec.roundToEven, yPos.negate(), zZero)
    wPos.isZero should be (true)
    wPos.sgn    should be (0)
  }

  test("add: x + 0 is rounded into the output spec (FP32 -> BFloat16)") {
    val zero32 = RealGeneric.zero(f32)
    val nzero32 = RealGeneric.zero(f32, 1)

    // 1 + 2^-8 : the dropped bits are exactly half an lsb of BFloat16, so
    // round-to-even keeps the even mantissa while `round` rounds up.
    val xTie = new RealGeneric(f32, 0, f32.exBias, SafeLong(1) << 15)
    xTie.add(bf16, RoundSpec.roundToEven, zero32).value  should be (SafeLong(0x3F80))
    xTie.add(bf16, RoundSpec.roundToEven, nzero32).value should be (SafeLong(0x3F80))
    zero32.add(bf16, RoundSpec.roundToEven, xTie).value  should be (SafeLong(0x3F80))
    xTie.add(bf16, RoundSpec.round, zero32).value        should be (SafeLong(0x3F81))
    xTie.add(bf16, RoundSpec.truncate, zero32).value     should be (SafeLong(0x3F80))

    // one ulp above the tie always rounds up
    val xUp = new RealGeneric(f32, 0, f32.exBias, (SafeLong(1) << 15) + 1)
    xUp.add(bf16, RoundSpec.roundToEven, zero32).value should be (SafeLong(0x3F81))
    zero32.add(bf16, RoundSpec.roundToEven, xUp).value should be (SafeLong(0x3F81))
  }

  test("add: x + 0 overflows to Inf when it does not fit the output spec") {
    val zero32 = RealGeneric.zero(f32)
    // the largest finite FP32 rounds up to Inf in BFloat16
    val xMax = new RealGeneric(f32, SafeLong(0x7F7FFFFFL))
    xMax.add(bf16, RoundSpec.roundToEven, zero32).value should be (SafeLong(0x7F80))
    zero32.add(bf16, RoundSpec.roundToEven, xMax).value should be (SafeLong(0x7F80))
    xMax.negate().add(bf16, RoundSpec.roundToEven, zero32).value should be (SafeLong(0xFF80))

    // the exponent is out of the FP16 range
    val big32 = new RealGeneric(f32, 0, f32.exBias+20, SafeLong(0))
    big32.add(f16, RoundSpec.roundToEven, zero32).value should be (SafeLong(0x7C00))
  }

  test("add: x + 0 underflows to a signed zero when it does not fit the output spec") {
    val zero64 = RealGeneric.zero(f64)
    val xTiny  = RealGeneric.fromDouble(f64, 1.0e-300)
    xTiny.add(f32, RoundSpec.roundToEven, zero64).value should be (SafeLong(0))
    zero64.add(f32, RoundSpec.roundToEven, xTiny).value should be (SafeLong(0))
    xTiny.negate().add(f32, RoundSpec.roundToEven, zero64).value should be (SafeLong(0x80000000L))
    zero64.add(f32, RoundSpec.roundToEven, xTiny.negate()).value should be (SafeLong(0x80000000L))

    // BFloat16 zero (disableSubnormal, so the mantissa bits are ignored) plus a
    // BFloat16 value whose exponent is below the FP16 range : AddFPGeneric
    // returns -0 here.
    val x = new RealGeneric(bf16, SafeLong(0x8032))
    val y = new RealGeneric(bf16, SafeLong(0x9FB5))
    x.add(f16, RoundSpec.roundToEven, y).value should be (SafeLong(0x8000))
  }
}
