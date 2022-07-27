

import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}



//import scopt.OptionParser

import scala.util.Random
import scala.math._
import rial.arith._
import rial.util.ScalaUtil._

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
}
