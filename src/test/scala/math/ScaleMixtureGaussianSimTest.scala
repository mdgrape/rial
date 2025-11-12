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
      println(f"1 / σ_a^2 = ${1.0 / sgmA2}")

      val refFunc = (w: Double) => {
        val g = exp(-(w * w) / (2*sgmP2))
        if(g == 0.0) {
          -(w / sgmA2)
        } else {
          -(w / sgmA2) * ((sgmA2 / sgmP2) / (sgmB / (g * sgmA) + 1.0) + 1.0)
        }
      }

      for(i <- 1 to n) {
        val x    = generator(spec,r)
        val x0   = x.toDouble
        val z0   = refFunc(x0)

        val z0r  = new RealGeneric(spec, z0)
        val zi   = ScaleMixtureGaussianSim.scaleMixtureGaussianSimGeneric( t, x, sgmA, sgmB, false )
        val zd   = zi.toDouble
        val errf = zd - z0r.toDouble
        val erri = errorLSB(zi, z0r.toDouble).toLong
        //println(f"${x.value.toLong}%x $z0 ${zi.toDouble}")
        if ((z0r.value - zi.value).abs > tolerance && !z0.isNaN && !z0.isInfinite && (abs(z0) < pow(2.0, spec.exMax+1))) {
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
            f"smg(x = ${x.toDouble} : ${xsgn}|${xexp}(${xexp - spec.exBias})|${xman.toLong.toBinaryString}(${xman.toLong}%x))")
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
            val errkey = log2DownL(erri.abs) + 1
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
      println(f"${generatorStr} Summary: total ${n} cases")
      if(maxError != 0.0) {
        println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
              + f"${zatMaxError} != ${refFunc(xatMaxError)}, "
              + f"diff = ${zatMaxError - refFunc(xatMaxError)}, x = ${xatMaxError}")
      }
      for(kv <- errs.toSeq.sortBy(_._1)) {
        val (k, (errPos, errNeg)) = kv
        println(f"N=$n%d : +/- ${k}%2d bit errors : positive $errPos%d / negative $errNeg%d")
      }
      println( "---------------------------------------------------------------")
    }
  }

  val sgmA16 = exp(-1.0)
  val sgmB16 = exp(-6.0)
  val smg16TableMaxBitDigit = ScaleMixtureGaussianSim.tableMaxBitDigit(sgmA16, sgmB16)

  val smg16F32TableI = ScaleMixtureGaussianSim.tableGeneration( 2, 8, 23, 23+8, sgmA16, sgmB16 )
  println(f"smg16TableMaxBitDigit = ${smg16TableMaxBitDigit}")
  val transitionPointF32 = pow(2.0, ScaleMixtureGaussianSim.tableDomainDigit(23, sgmA16, sgmB16))

  smgTest(smg16F32TableI, sgmA16, sgmB16, RealSpec.Float32Spec, n, r,
    f"Test Within (-${transitionPointF32},${transitionPointF32}) where table is used",
    generateRealWithin(transitionPointF32,_,_), pow(2.0, smg16TableMaxBitDigit).toInt)
  smgTest(smg16F32TableI, sgmA16, sgmB16, RealSpec.Float32Spec, n, r,
    "Test All range",generateRealFull(_,_), pow(2.0, smg16TableMaxBitDigit).toInt )

  val smg16B16TableI = ScaleMixtureGaussianSim.tableGeneration( 0, 7, 7, 7+10, sgmA16, sgmB16 )
  val transitionPointB16 = pow(2.0, ScaleMixtureGaussianSim.tableDomainDigit(RealSpec.BFloat16Spec.manW, sgmA16, sgmB16))

  smgTest(smg16B16TableI, sgmA16, sgmB16, RealSpec.BFloat16Spec, n, r,
    f"Test Within (-${transitionPointB16},${transitionPointB16}) where table is used",
    generateRealWithin(transitionPointB16,_,_), pow(2.0, smg16TableMaxBitDigit).toInt)
  smgTest(smg16B16TableI, sgmA16, sgmB16, RealSpec.BFloat16Spec, n, r,
    "Test All range",generateRealFull(_,_), pow(2.0, smg16TableMaxBitDigit).toInt )

}
