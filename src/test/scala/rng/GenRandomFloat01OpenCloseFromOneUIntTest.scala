import org.scalatest._

import chisel3._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.arith._
import rial.rng._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.language.reflectiveCalls

//
// Test GenRandomFloat01 using chi-squared test.
//

class Uniform01OpenCloseTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test uniform real distribution (0, 1]"

  val r = new Random(123456789)

  def generateRandomUInt(width: Int): BigInt = {
    val rnd = BigInt(r.nextLong)
    val mask = (BigInt(1) << width) - BigInt(1)
    rnd & mask
  }

  private def runtest ( rndW: Int, spec: RealSpec, stage: PipelineStageConfig ) = {

    val n = 1000
    val total = stage.total
    val pipeconfig = stage.getString

    it should f"uniform01OpenClose(x) pipereg ${pipeconfig} rndW ${rndW} spec ${spec.toStringShort}" in {
      test( new GenRandomFloat01OpenCloseFromOneUInt(rndW, spec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val nstage = stage.total
          val zs = scala.collection.mutable.ArrayBuffer.empty[Double]

          val zmax   = 1.0
          val zmin   = 0.0
          val zrange = zmax - zmin

          // Since we use FP that is represented in the binary numeral system,
          // we need to split the range by a power of 2. Otherwise, the min and
          // max numbers of a bin does not align to the FP precision and the
          // "true" number of FP numbers in a bin differs between each other.
          // This cause significant deviation from the uniform distribution.
          val nbins = 16
          val dz    = zrange / nbins
          val nref  = n.toDouble / nbins

          val ndeg      = nbins - 1
          val threshold = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          for(i <- 1 to n+nstage) {

            c.io.rnd.poke(generateRandomUInt(rndW).U(rndW.W))
            val zi = c.io.z.peek().litValue

            if (i > nstage) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble <= zmin || zmax < zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              assert(zmin < zd.toDouble && zd.toDouble <= zmax)

              zs += zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = zmin +  i    * dz
            val maxRange = zmin + (i+1) * dz

            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

            val term = (nsamples - nref) * (nsamples - nref) / nref
            chi2 += term
          }
          chi2 /= ndeg

          println(f"-----------------------------------------------")
          println(f"rndW = ${rndW}, float = ${spec.W}(ex=${spec.exW}, man=${spec.manW})")
          for(i <- 0 until nbins) {
            val minRange = zmin +  i    * dz
            val maxRange = zmin + (i+1) * dz
            val nsamples = zs.count(a => (minRange <= a && a < maxRange))
            println(f"n in [${minRange}%6f, ${maxRange}%6f) = ${nsamples} should be ${nref}")
          }
          println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")

          assert(chi2 < threshold)
        }
      }
    }
  }

  runtest(32, RealSpec.BFloat16Spec, PipelineStageConfig.none)
  runtest(32, RealSpec.Float32Spec,  PipelineStageConfig.none)
  runtest(64, RealSpec.Float32Spec,  PipelineStageConfig.none)
}

