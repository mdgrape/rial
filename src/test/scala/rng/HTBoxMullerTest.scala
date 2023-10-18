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

//
// Test High-Throughput BoxMuller
//
class HTBoxMullerTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test normal distribution, mean = 0, stddev = 1"

  val r = new Random(123456789)

  def generateRandomUInt(width: Int) = {
    val rnd = SafeLong(r.nextLong)
    rnd & maskSL(width)
  }

  private def runChiSquared( cfg: HTBoxMullerConfig, n: Int, description: String ) = {
    it should f"HTBoxMuller(x) Chi^2 test ${description}" in {
      test( new HTBoxMuller(cfg)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val rndW = cfg.rndW
          val realSpec = cfg.realSpec
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

          for(i <- 1 to n + c.nStage) {

            c.io.out.ready.poke(true.B)
            c.io.in.valid.poke(true.B)
            c.io.in.x.poke(generateRandomUInt(rndW).toBigInt.U(rndW.W))
            c.io.in.y.poke(generateRandomUInt(rndW).toBigInt.U(rndW.W))

            val valid = c.io.out.valid.peek().litValue == 1
            val z1    = c.io.out.z1.peek().litValue
            val z2    = c.io.out.z2.peek().litValue

            if (valid) {
              val z1r = new RealGeneric(realSpec, z1)
              val z2r = new RealGeneric(realSpec, z2)

              if(z1r.toDouble < xmin || xmax <= z1r.toDouble){
                println(f"outlier z = (${z1r.sgn}|${z1r.ex}(${z1r.ex-realSpec.exBias})|${z1r.man.toLong.toBinaryString}) = ${z1r.toDouble}")
              }
              zs = zs :+ z1r.toDouble

              if(z2r.toDouble < xmin || xmax <= z2r.toDouble){
                println(f"outlier z = (${z2r.sgn}|${z2r.ex}(${z2r.ex-realSpec.exBias})|${z2r.man.toLong.toBinaryString}) = ${z2r.toDouble}")
              }
              zs = zs :+ z2r.toDouble
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

  val rndW          = 32
  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3

  val cfgFP32 = HTBoxMullerConfig(
      rndW,
      realSpec = RealSpec.Float32Spec,
      polySpec = new PolynomialSpec(RealSpec.Float32Spec.manW, nOrderFP32, adrWFP32, extraBitsFP32),
      int2floatStge = PipelineStageConfig.atOut(1),
      polyPreStage  = PipelineStageConfig.atOut(1),
      polyCalcStage = PipelineStageConfig.atOut(3),
      polyPostStage = PipelineStageConfig.atOut(1),
      i2fPolyGap    = true,
      preCalcGap    = true,
      tableCalcGap  = true,
      calcPostGap   = true,
      mulStage = PipelineStageConfig.atOut(3)
    )

  runChiSquared(cfgFP32, 5000, "FP32, rndW=32")

  val nOrderBF16    = 0
  val adrWBF16      = 7
  val extraBitsBF16 = 1

  val cfgBF16 = HTBoxMullerConfig(
      rndW,
      realSpec = RealSpec.Float32Spec,
      polySpec = new PolynomialSpec(RealSpec.Float32Spec.manW, nOrderBF16, adrWBF16, extraBitsBF16),
      int2floatStge = PipelineStageConfig.none,
      polyPreStage  = PipelineStageConfig.none,
      polyCalcStage = PipelineStageConfig.none,
      polyPostStage = PipelineStageConfig.none,
      i2fPolyGap    = true,
      preCalcGap    = true,
      tableCalcGap  = true,
      calcPostGap   = true,
      mulStage = PipelineStageConfig.atOut(1)
    )

  runChiSquared(cfgBF16, 5000, "BF16, rndW=32")
}
