package rial.tests

import org.scalatest.funsuite.AnyFunSuite
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.util.Random
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math._
import rial.util.ScalaUtil._
import rial.arith._
import rial.table._

class SoftPlusSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateReal1to4( spec: RealSpec, r : Random) = {
    val rD : Double = r.nextDouble()*3.0+1.0
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  var counter = 0
  val specialValues = Seq(
      0.0,
     -0.0,
      1.0,
     -1.0,
      2.0,
     -2.0,
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }


  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def softplusTest(t: Seq[FuncTableInt], spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric),
    tolerance: Double) = {
    test(s"softplus(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      val softplus = (x: Double) => {
        if(x > log(Double.MaxValue)) { // exp overflows
          x
        } else {
          log(1.0 + exp(x))
        }
      }

      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = softplus(x0)
        val z0r  = new RealGeneric(spec, z0)
        val zi   = SoftPlusSim.softplusSimGeneric( t, x )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if (errf.abs.abs > tolerance) {
          val xsgn = bit(spec.W-1, x.value).toInt
          val xexp = slice(spec.manW, spec.exW, x.value)
          val xman = x.value & maskSL(spec.manW)

          val ztestsgn = bit(spec.W-1, zi.value).toInt
          val ztestexp = slice(spec.manW, spec.exW, zi.value)
          val ztestman = zi.value & maskSL(spec.manW)

          val zrefsgn = bit(spec.W-1, z0r.value).toInt
          val zrefexp = slice(spec.manW, spec.exW, z0r.value)
          val zrefman = z0r.value & maskSL(spec.manW)

          println(f"x0  = ${x0} (${x.sgn}|${x.ex}|${x.man.toLong.toBinaryString}|)")
          println(f"ref = ${z0}")
          println(f"sim = ${zd}")
          println(f"test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman}(${zrefman.toLong}%x)) = softplus(x = ${xsgn}|${xexp}(${xexp - spec.exBias})|${xman}(${xman.toLong}%x))")
        }

        if (zi.isNaN) {
          assert(x0.isNaN)
        } else if (zi.isInfinite) {
          assert(zi.sgn == 0, "softplus takes positive inf only")
          assert(x.isInfinite, "softplus takes +inf only when x is inf")
        } else {
          assert(errf.abs<=tolerance, f"SoftPlusSimTest: x = ${x0}, zref = ${z0}, zsim = ${zd}, diff = ${errf}")

          if (maxError < errf.abs) {
            maxError    = errf.abs
            xatMaxError = x0
            zatMaxError = zd
          }
        }
      }
      println(f"${generatorStr} Summary")
      if(maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError} where the value is "
              + f"${zatMaxError} != ${softplus(xatMaxError)}, "
              + f"diff = ${zatMaxError - softplus(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- ${k}%4d errors (${log2DownL(k)+1}%2d ULPs) positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val softplusF32TableI = SoftPlusSim.tableGeneration( 2, 8, 23, 23+2 )

  softplusTest(softplusF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1e-5)
  softplusTest(softplusF32TableI, RealSpec.Float32Spec, n, r,
    "Test All range",generateRealFull(_,_), 1e-5)
  softplusTest(softplusF32TableI, RealSpec.Float32Spec, n, r,
    "Test special value",generateSpecialValues(_,_), 1e-5)

//   val softplusBF16TableI = SoftPlusSim.tableGeneration(0, 7, 7, 7 ) // [1,2) + [2,4) + 1.0
// 
//   softplusTest(softplusBF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
//   softplusTest(softplusBF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test All range",generateRealFull(_,_), 1 )
//   softplusTest(softplusBF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test special value",generateSpecialValues(_,_), 1 )
// 
//   val float48Spec = new RealSpec(10, 511, 37)
//   val softplusFP48TableI = SoftPlusSim.tableGeneration(3, 10, 37, 37+4 )
// 
//   softplusTest(softplusFP48TableI, float48Spec, n, r,
//     "Test Within (-128,128)",generateRealWithin(128.0,_,_), 1)
//   softplusTest(softplusFP48TableI, float48Spec, n, r,
//     "Test All range",generateRealFull(_,_), 1 )
//   softplusTest(softplusFP48TableI, float48Spec, n, r,
//     "Test special value",generateSpecialValues(_,_), 1 )
// 
//   val softplusFP64TableI = SoftPlusSim.tableGeneration(3, 10, 52, 52+1 )
// 
//   softplusTest(softplusFP64TableI, RealSpec.Float64Spec, n, r,
//     "Test Within (-128,128)",generateRealWithin(128.0,_,_), 3) // XXX 3!
//   softplusTest(softplusFP64TableI, RealSpec.Float64Spec, n, r,
//     "Test All range",generateRealFull(_,_), 3 )
//   softplusTest(softplusFP64TableI, RealSpec.Float64Spec, n, r,
//     "Test special value",generateSpecialValues(_,_), 3 )
}
