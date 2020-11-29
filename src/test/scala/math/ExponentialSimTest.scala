
//package rial.tests

//import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.FunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

//import org.scalatest.flatspec.AnyFlatSpec
//import org.scalatest.matchers.should.Matchers
//import scopt.OptionParser

import scala.util.Random
import scala.math._
import rial.math._
import rial.util.ScalaUtil._

class ExponentialSimTest extends FunSuite with BeforeAndAfterAllConfigMap {
//class ExponentialSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

  var err1lsbPos = 0
  var err1lsbNeg = 0
  test("Compare my Pow2F32Sim with pow(2,x) for proper range") {
    for(i <- 1 to n) {
      val x  = ((r.nextDouble()-0.5)*128.0).toFloat
      val z0 = pow(2.0, x).toFloat
      val z0i = java.lang.Float.floatToRawIntBits(z0)
      val xi = java.lang.Float.floatToRawIntBits(x)
      val zi = ExponentialSim.pow2F32Sim(xi).toInt
      val z  = java.lang.Float.intBitsToFloat(zi)
      val errf = z-z0
      val erri = zi-z0i
      if (erri.abs>1) {
        println(f"Error more than 1 LSB : $x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      } else if (erri==1) err1lsbPos+=1
      else if (erri== -1) err1lsbNeg+=1
      assert(erri.abs<=1)
      //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
    }
    println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
  }

  err1lsbPos = 0
  err1lsbNeg = 0
  test("Compare my Pow2F32Sim with pow(2,x) for all range") {
    for(i <- 1 to n) {
      val xi = r.nextInt()
      val x  = java.lang.Float.intBitsToFloat(xi)
      val z0 = normalizeFloat(pow(2.0, x).toFloat)
      val z0i = java.lang.Float.floatToRawIntBits(z0)
      val zi = ExponentialSim.pow2F32Sim(xi).toInt
      val z  = java.lang.Float.intBitsToFloat(zi)
      val errf = z-z0
      val erri = zi-z0i
      if (erri.abs>1) {
        println(f"Error more than 1 LSB : $x%14.7e($xi%08x) : $z0%14.7e($z0i%08x) $z%14.7e($zi%08x) $errf%14.7e $erri%d")
      } else if (erri==1) err1lsbPos+=1
      else if (erri== -1) err1lsbNeg+=1
      assert(erri.abs<=1)
    }
    println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
  }
  
}
