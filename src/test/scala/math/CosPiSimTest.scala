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

class CosPiSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {

  var n = 1000000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRealWithin( from : Double, to : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = from + r.nextDouble() * (to - from)
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  def cosPiTest(t: Seq[FuncTableInt], spec : RealSpec, n : Int, r : Random,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric),
    tolerance : Int) = {
    test(s"cosPi(x), format ${spec.toStringShort}, ${generatorStr}") {
      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = math.cos(math.Pi * x0)
        val z0r  = new RealGeneric(spec, z0)
        val zi   = CosPiSim.cosPiSimGeneric( t, x )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble)
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if (erri.abs > 3) { // error exceeds 2 bit
          val xsgn = bit(spec.W-1, x.value).toInt
          val xexp = slice(spec.manW, spec.exW, x.value)
          val xman = x.value & maskSL(spec.manW)

          val ztestsgn = bit(spec.W-1, zi.value).toInt
          val ztestexp = slice(spec.manW, spec.exW, zi.value)
          val ztestman = zi.value & maskSL(spec.manW)

          val zrefsgn = bit(spec.W-1, z0r.value).toInt
          val zrefexp = slice(spec.manW, spec.exW, z0r.value)
          val zrefman = z0r.value & maskSL(spec.manW)

          println(f"x   = ${x0}, cos(x) = ${math.cos(Pi * x0)}")
          println(f"ref = ${z0}")
          println(f"sim = ${zd}")
          println(f"test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman.toLong.toBinaryString}(${ztestman.toLong}%x)) != ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x))")
        }

        if (x0.isInfinity) {
          assert(zi.isNaN)
        } else if (x0.isNaN) {
          assert(zi.isNaN)
        } else {
          if(erri != 0) {
            val errkey = erri.abs.toLong
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
            maxError = erri.abs
            xatMaxError = x0
            zatMaxError = zd
          }
        }
      }
      println(f"${generatorStr} Summary")
      if(maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
              + f"${zatMaxError} != ${acos(xatMaxError)}, "
              + f"diff = ${zatMaxError - acos(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- ${k}%4d errors (${log2DownL(k)+1}%2d bits in the LSB side) positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val cosPiF32TableI = CosPiSim.cosPiTableGeneration( 2, 8, 23, 23+2 )

  cosPiTest(cosPiF32TableI, RealSpec.Float32Spec, n, r,
     "Test Within [-1, 0]", generateRealWithin(-1.0, 0.0,_,_), 3)
  cosPiTest(cosPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [0, 0.5]", generateRealWithin(0.0, 0.5-pow(2.0, -23),_,_), 3)
  cosPiTest(cosPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [0.5-2^-9, 0.5]", generateRealWithin(0.5 - pow(2.0, -9), 0.5-pow(2.0, -23),_,_), 3)
  cosPiTest(cosPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [0.5-2^-5, 0.5]", generateRealWithin(0.5 - pow(2.0, -5), 0.5-pow(2.0, -23),_,_), 3)

  cosPiTest(cosPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [2^-1, 1]", generateRealWithin(0.5, 1.0,_,_), 3)
  cosPiTest(cosPiF32TableI, RealSpec.Float32Spec, n, r,
    "Test Within [1, 2]", generateRealWithin(1.0, 2.0,_,_), 3)

//   val cosPiBF16TableI = CosPiSim.cosPiTableGeneration(0, 7, 7, 7 ) // [1,2) + [2,4) + 1.0
//
//   cosPiTest(cosPiBF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test Within (-128,128)",generateRealWithin(128.0,_,_))
//   cosPiTest(cosPiBF16TableI, RealSpec.BFloat16Spec, n, r,
//     "Test All range",generateRealFull(_,_) )
}
