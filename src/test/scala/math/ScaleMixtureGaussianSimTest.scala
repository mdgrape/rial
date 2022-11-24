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

class ScaleMixtureGaussianSimTest extends AnyFunSuite with BeforeAndAfterAllConfigMap {
  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*2.0*p
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

  def smgTest(t: FuncTableInt, sgmA: Double, sgmB: Double, spec: RealSpec,
    n: Int, r: Random,
    generatorStr: String, generator: ( (RealSpec, Random) => RealGeneric),
    tolerance: Int) = {
    test(s"scaleMixtureGaussian(x), format ${spec.toStringShort}, ${generatorStr}") {
      counter = 0

      var maxError   = 0.0
      var xatMaxError = 0.0
      var zatMaxError = 0.0

      val errs = collection.mutable.Map.empty[Long, (Int, Int)]

      val sgmA2 = sgmA * sgmA
      val sgmB2 = sgmB * sgmB
      val sgmP2 = (sgmA2 * sgmB2) / (sgmA2 - sgmB2)

      println(f"=============================================================")
      println(f"σ_a^2 = ${sgmA2}")
      println(f"σ_b^2 = ${sgmB2}")
      println(f"σ'^2  = ${sgmP2}")
      println(f"max(1st term) = ${1.0 / (sgmB / sgmA + 1.0)}")
      println(f"2nd term = ${sgmP2 / sgmA2}")
      val refFunc = (w: Double) => {
        val g = exp(-w * w / (2*sgmP2))
        -(w / sgmP2) * (1.0 / (sgmB / (g * sgmA) + 1.0) + sgmP2 / sgmA2)
      }

      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = refFunc(x0)

        val z0r  = new RealGeneric(spec, z0)
        val zi   = ScaleMixtureGaussianSim.scaleMixtureGaussianSimGeneric( t, x, sgmA, sgmB, true )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toLong
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if ((z0r.value - zi.value).abs > tolerance) {
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
          println(f"test(${ztestsgn}|${ztestexp}(${ztestexp - spec.exBias})|${ztestman.toLong.toBinaryString}(${ztestman.toLong}%x)) != " +
            f"ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString}(${zrefman.toLong}%x)) that is the result of " +
            f"smg(x = ${x.toDouble}:${xsgn}|${xexp}(${xexp - spec.exBias})|${xman.toLong.toBinaryString}(${xman.toLong}%x))")
        }

        if (x.isInfinite) {
          assert(zi.isInfinite && x.sgn != zi.sgn)
        } else if (x.isNaN) {
          assert(zi.isNaN)
        } else if (x.isZero) {
          assert(zi.isZero)
        } else {
          assert(erri.abs<=tolerance)

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

          if (maxError < erri.abs) {
            maxError    = erri.abs.toDouble
            xatMaxError = x0
            zatMaxError = zd
          }
        }
      }
      println(f"${generatorStr} Summary")
      if(maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
              + f"${zatMaxError} != ${refFunc(xatMaxError)}, "
              + f"diff = ${zatMaxError - refFunc(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- ${k}%4d errors (${log2DownL(k)+1}%2d ULPs) positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val sgmA16 = exp(-1.0)
  val sgmB16 = exp(-6.0)
  val smg16F32TableI = ScaleMixtureGaussianSim.tableGeneration( 2, 8, 23, 23+3, sgmA16, sgmB16 )

  val transitionPoint = ScaleMixtureGaussianSim.tableDomain(23, sgmA16, sgmB16)

  smgTest(smg16F32TableI, sgmA16, sgmB16, RealSpec.Float32Spec, n, r,
    "Test Within (-1,1)", generateRealWithin(transitionPoint,_,_), pow(2.0, 18).toInt)
//   smgTest(smg16F32TableI, sgmA16, sgmB16, RealSpec.Float32Spec, n, r,
//     "Test All range",generateRealFull(_,_), pow(2.0, 9).toInt )
//   smgTest(smg16F32TableI, sgmA16, sgmB16, RealSpec.Float32Spec, n, r,
//     "Test special value",generateSpecialValues(_,_), pow(2.0, 9).toInt )

  val smg16BF16TableI = ScaleMixtureGaussianSim.tableGeneration( 0, 7, 7, 7+3, sgmA16, sgmB16 )

//   smgTest(smg16BF16TableI, sgmA16, sgmB16, RealSpec.BFloat16Spec, n, r,
//     "Test Within (-128,128)",generateRealWithin(128.0,_,_), pow(2.0, 9))
//   smgTest(smg16BF16TableI, sgmA16, sgmB16, RealSpec.BFloat16Spec, n, r,
//     "Test All range",generateRealFull(_,_), pow(2.0, 9) )
//   smgTest(smg16BF16TableI, sgmA16, sgmB16, RealSpec.BFloat16Spec, n, r,
//     "Test special value",generateSpecialValues(_,_), pow(2.0, 9) )
}
