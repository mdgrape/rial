package rial.tests

import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._
//import org.scalatest.flatspec.AnyFlatSpec
//import org.scalatest.matchers.should.Matchers
//import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}
import org.scalatest.FlatSpec
import org.scalatest.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}
import rial.math._
import rial.util.ScalaUtil._

import scala.util.Random
import scala.math._
import scala.language.reflectiveCalls
//
// Testing Pow2F32 using ChiselTest
//

class Pow2F32Test extends FlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test 2^x, float32"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

  it should "test 2^x in float" in {
    test(new Pow2F32) { c =>
      {
        for(i <- 1 to n) {
          val x  = ((r.nextDouble()-0.5)*128.0).toFloat
          val xi = java.lang.Float.floatToRawIntBits(x)
          val xl = xi.toLong & mask(32)
          c.io.x.poke(xl.U(32.W))
          val zi = c.io.z.peek.litValue.toInt
          val z  = java.lang.Float.intBitsToFloat(zi)
          val z0i= ExponentialSim.pow2F32Sim(xi).toInt
          val z0 = java.lang.Float.intBitsToFloat(zi)
          val errf = z-z0
          val erri = zi-z0i
          c.io.z.expect(z0i.U(32.W))
          c.clock.step(1)
        }
      }
    }
  }
//    iotesters.Driver.execute(arg0, () => new InvSqrt(25, 6, 0) ) {
//      c => new ExponentialUnitTester(c, ncycle)
//    }
}

