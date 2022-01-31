
//package rial.tests

//import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

//import org.scalatest.flatspec.AnyFlatSpec
//import org.scalatest.matchers.should.Matchers
//import scopt.OptionParser

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

import rial.testUtil.ScalaTestUtil._

class LogSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
//class ExponentialSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

  def log2Test(sim: Log2Sim, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    val spec = sim.spec
    test(s"Log base 2 , format ${spec.toStringShort}, ${generatorStr}") {
      var err1lsbPos = 0
      var err1lsbNeg = 0
      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble
        val z0 = log(x0)/log(2.0)
        val z0r = new RealGeneric(spec, z0)
        //println(f"${x0}%f ${z0}%f ${z0r.value.toLong}%x")
        val zi = sim.eval( x )
        val zd = zi.toDouble
        val errf = zd-z0
        val erri = errorLSB(zi, z0)
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if (z0r.value != zi.value) {
          println(f"x=${x.value.toLong}%x zRef=${z0r.value.toLong}%x z=${zi.value.toLong}%x")
        }
        if (z0.isInfinity) {
          assert(zd.isInfinity)
        } else if (zd.isInfinite) {
          assert(z0r.isInfinite)
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if (erri.abs>=2.0) {
            println(f"Error more than 2 LSB : ${x.toDouble}%14.7e : $z0%14.7e ${zi.toDouble}%14.7e $errf%14.7e $erri%f")
          } else if (erri>=1.0) err1lsbPos+=1
          else if (erri<= -1.0) err1lsbNeg+=1
          //println(f"erri=${erri}")
          assert(erri.abs<=2)
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
    }
  }

  log2Test(new Log2Sim( RealSpec.Float32Spec, 2, 8, 2 ), n, r,
    "2^x Test Within (-128,128)",generateRealWithinPos(128.0,_,_))
  log2Test(new Log2Sim( RealSpec.Float32Spec, 2, 8, 2 ), n, r,
    "2^x Test All range",generateRealFullPos(_,_) )

  log2Test(new Log2Sim( RealSpec.BFloat16Spec, 0, 7, 3 ), n, r,
    "Test Within (-128,128)",generateRealWithinPos(128.0,_,_))
  log2Test(new Log2Sim( RealSpec.BFloat16Spec, 0, 7, 3 ), n, r,
    "Test All range",generateRealFullPos(_,_) )

  /*
  expTest(log2F32TableI, log2F32ExtraBits, RealSpec.Float32Spec, n, r,
    "exp Test Within (-128,128)",generateRealWithin(128.0,_,_))
  expTest(log2F32TableI, log2F32ExtraBits, RealSpec.Float32Spec, n, r,
    "exp Test All range",generateRealFull(_,_) )

   */

}
