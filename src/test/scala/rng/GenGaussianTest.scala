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
import spire.math.Numeric
import spire.implicits._

import rial.arith._
import rial.math._
import rial.rng._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

class GenGaussianTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test normal distribution, mean = 0, stddev = 1"

  var n = 100000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("100000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  private def runTest(
    rndW: Int,
    spec: RealSpec,
    nOrder: Int,
    polyStage: Int,
    mulStage: Int,
  ) = {
    it should f"GenGaussian(x) chi^2 w/ ${spec.toStringShort}, stage=(${polyStage}, ${mulStage})" in {
      test( new GenRandomGaussian(rndW, spec, nOrder, polyStage, mulStage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          var zs = scala.collection.mutable.ArraySeq.empty[Double]

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

          for(i <- 1 to n) {

            c.io.x.valid.poke(true.B)
            c.io.x.bits .poke(BigInt(rndW, r).U(rndW.W))

            val valid = c.io.valid.peek().litValue == 1
            val zi    = c.io.z.peek().litValue.toBigInt

            if (valid) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"outlier z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              zs = zs :+ zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          val cdf = (x: Double) => {
            0.5 * (1.0 + Erf.erf(x / sqrt(2.0)))
          }

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

  runTest(
    rndW      = 64,
    spec      = RealSpec.Float32Spec,
    nOrder    = 3,
    polyStage = 0,
    mulStage  = 0,
  )

  runTest(
    rndW      = 64,
    spec      = RealSpec.Float32Spec,
    nOrder    = 3,
    polyStage = 3,
    mulStage  = 2,

  )
}
