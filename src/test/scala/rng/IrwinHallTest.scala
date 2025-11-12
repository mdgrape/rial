package rial.tests

import org.scalatest._

import chisel3._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import org.apache.commons.math3.special.Erf

import spire.math.SafeLong

import rial.arith.{RealSpec, RealGeneric}
import rial.rng._
import rial.util.ScalaUtil._

import scala.util.Random
import scala.math._
import scala.language.reflectiveCalls

//
// Test IrwinHall distribution
//
class IrwinHallTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test normal distribution(m=0,s=1), generated based on Irwin-Hall distribution"

  var n = 100000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("100000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRandomUInt(width: Int): BigInt = {
    BigInt(width, r)
  }

  val inv_fact12 = 1.0 / Seq.tabulate(12)(i => i.toDouble+1.0).fold(1.0)(_*_)
  def cdf(xnorm: Double): Double = {
    val x = xnorm + 6.0 // recover Irwin-Hall distribution from normal distr
    Seq.tabulate(floor(x).toInt)(k => {
      val sgn = if(k % 2 == 1) { -1.0 } else { 1.0 }
      val coef = binomial(12, k)
      sgn * coef * pow(x - k, 12)
    }).fold(0.0)(_+_) * inv_fact12
  }
  def binomial(n: Int, k: Int): Double = {
    Seq.tabulate(k)(i => { (n - i).toDouble / (i + 1).toDouble }).fold(1.0)(_*_)
  }

  private def runChiSquared( spec: RealSpec, rndW: Int ) = {
    it should f"IrwinHall(x) spec ${spec.toStringShort} Chi^2 test" in {
      test( new IrwinHall(spec, rndW)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          var zs = scala.collection.mutable.ArrayBuffer.empty[Double]

          // +/- 3 sigma ~ 99.7%.
          // +/- 4 sigma ~ 99.994%.
          // +/- 5 sigma ~ 99.99994%.
          val xmax   =  5.0
          val xmin   = -5.0
          val xrange = xmax - xmin

          // Since we use FP that is represented in the binary numeral system,
          // we need to split the range by a power of 2. Otherwise, the min and
          // max numbers of a bin does not align to the FP precision and the
          // "true" number of FP numbers in a bin differs between each other.
          // This cause significant deviation from the uniform distribution.
          val nbins = 16
          val dx    = xrange / nbins

          val ndeg      = nbins - 1
          val threshold = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          c.io.en.poke(true.B) // always en
          for(t <- 1 to n * 2) {

            c.io.input.valid.poke(true.B)
            for(i <- 0 until c.io.input.rnds.length) {
              c.io.input.rnds(i).poke(generateRandomUInt(rndW).U)
            }

            val valid = c.io.output.valid.peek().litValue == 1
            val zi    = c.io.output.rnd.peek().litValue

            if (valid) {
              val zd = new RealGeneric(spec, zi)

              val manMask = (BigInt(1) << spec.manW) - BigInt(1)
              val zisgn = bit(spec.W-1, zi).toInt
              val ziexp = slice(spec.manW, spec.exW, zi)
              val ziman = zi & manMask
              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                println(f"outlier z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              } else {
                // println(f"z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              zs += zd.toDouble
            }
            // println(f"-----------------------------------------------")
            c.clock.step(1)
          }

          // calculate chi^2

          val nz = zs.length

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx
            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

            val nref = (cdf(maxRange) - cdf(minRange)) * nz

            val term = (nsamples - nref) * (nsamples - nref) / nref
            chi2 += term
          }
          chi2 /= ndeg

//           if(threshold < chi2) {
            println(f"-----------------------------------------------")
            println(f"${nz} real values generated")
            for(i <- 0 until nbins) {
              val minRange = xmin +  i    * dx
              val maxRange = xmin + (i+1) * dx
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))

              val nref = (cdf(maxRange) - cdf(minRange)) * nz
              println("n in [%8.3f, %8.3f) = %8d should be %10.3f".format(minRange, maxRange, nsamples, nref))
            }
            println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")
//           }

          assert(chi2 < threshold)
        }
      }
    }
  }

  runChiSquared(RealSpec.Float32Spec, 32)
}

//
// Test IrwinHall gaussian generator
//
class IrwinHallGeneratorTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Irwin-Hall generator with ready/valid I/F"

  var n = 100000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("100000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRandomUInt(width: Int): BigInt = {
    BigInt(width, r)
  }

  val inv_fact12 = 1.0 / Seq.tabulate(12)(i => i.toDouble+1.0).fold(1.0)(_*_)
  def cdf(xnorm: Double): Double = {
    val x = xnorm + 6.0 // recover Irwin-Hall distribution from normal distr
    Seq.tabulate(floor(x).toInt)(k => {
      val sgn = if(k % 2 == 1) { -1.0 } else { 1.0 }
      val coef = binomial(12, k)
      sgn * coef * pow(x - k, 12)
    }).fold(0.0)(_+_) * inv_fact12
  }
  def binomial(n: Int, k: Int): Double = {
    Seq.tabulate(k)(i => { (n - i).toDouble / (i + 1).toDouble }).fold(1.0)(_*_)
  }

  private def runChiSquared( spec: RealSpec, rndW: Int ) = {
    it should f"IrwinHall(x) spec ${spec.toStringShort} Chi^2 test" in {
      test( new IrwinHallGenerator(spec, rndW)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          var zs = scala.collection.mutable.ArrayBuffer.empty[Double]

          // +/- 3 sigma ~ 99.7%.
          // +/- 4 sigma ~ 99.994%.
          // +/- 5 sigma ~ 99.99994%.
          val xmax   =  5.0
          val xmin   = -5.0
          val xrange = xmax - xmin

          // Since we use FP that is represented in the binary numeral system,
          // we need to split the range by a power of 2. Otherwise, the min and
          // max numbers of a bin does not align to the FP precision and the
          // "true" number of FP numbers in a bin differs between each other.
          // This cause significant deviation from the uniform distribution.
          val nbins = 16
          val dx    = xrange / nbins

          val ndeg      = nbins - 1
          val threshold = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          for(t <- 1 to n * 2) {

            c.io.input.valid.poke(true.B)
            for(i <- 0 until c.io.input.bits.length) {
              c.io.input.bits(i).poke(generateRandomUInt(rndW).U)
            }

            c.io.output.ready.poke(true.B)
            val valid = c.io.output.valid.peek().litValue == 1
            val zi    = c.io.output.bits.peek().litValue

            if (valid) {
              val zd = new RealGeneric(spec, zi)

              val manMask = (BigInt(1) << spec.manW) - BigInt(1)
              val zisgn = bit(spec.W-1, zi).toInt
              val ziexp = slice(spec.manW, spec.exW, zi)
              val ziman = zi & manMask
              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                println(f"outlier z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              } else {
                // println(f"z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              zs += zd.toDouble
            }
            // println(f"-----------------------------------------------")
            c.clock.step(1)
          }

          // calculate chi^2

          val nz = zs.length

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx
            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

            val nref = (cdf(maxRange) - cdf(minRange)) * nz

            val term = (nsamples - nref) * (nsamples - nref) / nref
            chi2 += term
          }
          chi2 /= ndeg

//           if(threshold < chi2) {
            println(f"-----------------------------------------------")
            println(f"${nz} real values generated")
            for(i <- 0 until nbins) {
              val minRange = xmin +  i    * dx
              val maxRange = xmin + (i+1) * dx
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))

              val nref = (cdf(maxRange) - cdf(minRange)) * nz
              println("n in [%8.3f, %8.3f) = %8d should be %10.3f".format(minRange, maxRange, nsamples, nref))
            }
            println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")
//           }

          assert(chi2 < threshold)
        }
      }
    }
  }

  runChiSquared(RealSpec.Float32Spec, 32)
}
