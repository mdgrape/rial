
//import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.FunSuite
import org.scalatest.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

//import org.scalatest.flatspec.AnyFlatSpec
//import org.scalatest.matchers.should.Matchers
//import scopt.OptionParser

import scala.util.Random
import scala.math._
import rial.arith._
import rial.util.ScalaUtil._

class RealSimTest extends FunSuite with Matchers with BeforeAndAfterAllConfigMap {
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
  
}
