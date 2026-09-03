package rial.tests

import org.scalatest._

import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import org.scalatest.funspec.AnyFunSpec

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

class GenGaussianTest extends AnyFunSpec
    with ChiselSim with EmitVcd {


  addOption(CliOption.int("n", "number of test cases (default 100000)"))
  def n: Int = getOption[Int]("n").getOrElse(100000)

  // Use SplittableRandom rather than scala.util.Random (java.util.Random, a
  // 48-bit LCG): the LCG correlates the two 32-bit halves of each 64-bit draw,
  // which here correlates the sign bit (63) with the magnitude-bucket bits and
  // shows up as a spurious per-bin sign asymmetry in the chi^2.
  val sr = new java.util.SplittableRandom(123456789L)
  val mask64 = (BigInt(1) << 64) - 1
  def rnd64(): BigInt = BigInt(sr.nextLong()) & mask64

  private def runTest(
    rndW: Int,
    spec: RealSpec,
    nOrder: Int,
    polyStage: Int,
    mulStage: Int,
  ) = {
    it(f"GenGaussian(x) chi^2 w/ ${spec.toStringShort}, stage=(${polyStage}, ${mulStage})") {
      simulate(new GenRandomGaussian(rndW, spec, nOrder, polyStage, mulStage)) { c =>
        {
          val zs = scala.collection.mutable.ArrayBuffer.empty[Double]

          // +/- 3 sigma ~ 99.7%.
          // +/- 4 sigma ~ 99.994%.
          // +/- 5 sigma ~ 99.99994%.
          val xmax   =  5.0
          val xmin   = -5.0
          val xrange = xmax - xmin

          val nbins = 16
          val dx    = xrange / nbins

          // Pearson's chi^2 is only valid where the expected count is large
          // enough; bins with expected < 5 (the far tails at this n) are
          // excluded so a single Poisson fluctuation there doesn't dominate.
          val minExpected = 5.0
          val threshold   = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          for(i <- 1 to n) {

            c.io.rnd.valid.poke(true.B)
            c.io.rnd.bits .poke(rnd64().U(rndW.W))

            val valid = c.io.z.valid.peek().litValue == 1
            val zi    = c.io.z.bits .peek().litValue.toBigInt

            if (valid) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"outlier z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              zs += zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          val cdf = (x: Double) => {
            0.5 * (1.0 + Erf.erf(x / sqrt(2.0)))
          }

          val nz = zs.length

          var chi2     = 0.0
          var usedBins = 0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx
            val nref     = (cdf(maxRange) - cdf(minRange)) * nz

            if (nref >= minExpected) {
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))
              val term = (nsamples - nref) * (nsamples - nref) / nref
              chi2 += term
              usedBins += 1
            }
          }
          val ndeg = usedBins - 1
          chi2 /= ndeg

//           if(threshold < chi2) {
            println(f"-----------------------------------------------")
            println(f"${nz} real values generated")
            for(i <- 0 until nbins) {
              val minRange = xmin +  i    * dx
              val maxRange = xmin + (i+1) * dx
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))

              val nref = (cdf(maxRange) - cdf(minRange)) * nz
              val used = if (nref >= minExpected) "" else "  (excluded: exp<5)"
              println("n in [%8.3f, %8.3f) = %8d should be %10.3f%s".format(minRange, maxRange, nsamples, nref, used))
            }
            println(f"chi^2 = ${chi2}, dof = ${ndeg}, threshold(<3%%) = ${threshold}")
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

  // runTest(
  //   rndW      = 64,
  //   spec      = RealSpec.Float32Spec,
  //   nOrder    = 3,
  //   polyStage = 3,
  //   mulStage  = 2,
  //
  // )
}
