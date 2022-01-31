
//package rial.tests


import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}



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


class ExponentialSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
//class ExponentialSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def pow2Test(t: FuncTableInt, extraBits: Int, spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    test(s"Power of 2 , format ${spec.toStringShort}, ${generatorStr}") {
      var err1lsbPos = 0
      var err1lsbNeg = 0
      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble
        val z0 = pow(2.0, x0)
        val z0r = new RealGeneric(spec, z0)
        val zi = ExponentialSim.pow2simGeneric( t, extraBits, x )
        val zd = zi.toDouble
        val errf = zd-z0
        val erri = errorLSB(zi, z0)
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if (z0r.value != zi.value) {
          println(f"${x.value.toLong}%x ${z0r.value.toLong}%x ${zi.value.toLong}%x")
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
          assert(erri.abs<=1)
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
    }
  }

  def expTest(t: FuncTableInt, extraBits: Int, spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    test(s"Exponential, format ${spec.toStringShort}, ${generatorStr}") {
      var err1lsbPos = 0
      var err1lsbNeg = 0
      var lsbErr=4
      for(i <- 1 to n) {
        val x  = generator(spec,r)
        val x0 = x.toDouble
        val z0 = java.lang.Math.exp(x0)
        val z0r = new RealGeneric(spec, z0)
        val xlog2e = x0 /java.lang.Math.log(2.0d)
        val zPow2  = pow(2.0, xlog2e)
        println(f"${x0} ${z0} ${zPow2} ${z0-zPow2}")
        val zi = ExponentialSim.pow2simGeneric( t, extraBits, new RealGeneric(spec, xlog2e) )
        //val zi = ExponentialSim.expSimGeneric( t, extraBits, x )
        val zd = zi.toDouble
        val errf = zd-z0
        val erri = errorLSB(zi, z0)
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if (z0r.value != zi.value) {
          println(f"${x.value.toLong}%x ${z0r.value.toLong}%x ${zi.value.toLong}%x")
        }
        if (z0.isInfinity) {
          assert(zd.isInfinity)
        } else if (zd.isInfinite) {
          assert(z0r.isInfinite)
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if (erri.abs>=lsbErr) {
            println(f"Error more than ${lsbErr} LSB : ${x.toDouble}%14.7e : $z0%14.7e ${zi.toDouble}%14.7e $errf%14.7e $erri%f")
          } else if (erri>=1.0) err1lsbPos+=1
          else if (erri<= -1.0) err1lsbNeg+=1
          assert(erri.abs<lsbErr)
        }
        //println(f"$x%14.7e : $z0%14.7e $z%14.7e $errf%14.7e $erri%d")
      }
      println(f"N=$n%d : 1LSB errors positive $err1lsbPos%d / negative $err1lsbNeg%d")
    }
  }

  // fracW include extra bits added during calc.
  def pow2TableGeneration( order : Int, adrW : Int, fracW : Int ) = {
    val tableD = new FuncTableDouble( x => pow(2.0,x)-1.0, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  val pow2F32ExtraBits = 2
  val pow2F32TableI = pow2TableGeneration( 2, 8, 23+pow2F32ExtraBits )

  pow2Test(pow2F32TableI, pow2F32ExtraBits, RealSpec.Float32Spec, n, r,
    "2^x Test Within (-128,128)",generateRealWithin(128.0,_,_))
  pow2Test(pow2F32TableI, pow2F32ExtraBits, RealSpec.Float32Spec, n, r,
    "2^x Test All range",generateRealFull(_,_) )

  expTest(pow2F32TableI, pow2F32ExtraBits, RealSpec.Float32Spec, n, r,
    "exp Test Within (-128,128)",generateRealWithin(128.0,_,_))
  expTest(pow2F32TableI, pow2F32ExtraBits, RealSpec.Float32Spec, n, r,
    "exp Test All range",generateRealFull(_,_) )

  val pow2BF16ExtraBits = 1
  val pow2BF16TableI = pow2TableGeneration( 0, 8, 7 )

  pow2Test(pow2BF16TableI, 0, RealSpec.BFloat16Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  pow2Test(pow2BF16TableI, 0, RealSpec.BFloat16Spec, n, r,
    "Test All range",generateRealFull(_,_) )

}
