package rial.tests

import org.scalatest._

import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import org.scalatest.funspec.AnyFunSpec

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

class Uniform01CloseTest extends AnyFunSpec
    with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 4000)"))
  def n: Int = getOption[Int]("n").getOrElse(4000)

  val r = new Random(123456789)

  def generateRandomUInt(width: Int): BigInt = {
    val rnd = BigInt(r.nextLong())
    val mask = (BigInt(1) << width) - BigInt(1)
    rnd & mask
  }

  private def runtest ( rndW: Int, spec: RealSpec, stage: PipelineStageConfig,
    checkHold: Boolean = false ) = {

    val total = stage.total
    val pipeconfig = stage.getString

    it(f"uniform01Close(x) pipereg ${pipeconfig} rndW ${rndW} spec ${spec.toStringShort}") {
      simulate(new GenRandomFloat01Close(rndW, spec, stage)) { c =>
        {
          val nstage = stage.total
          // clock edges since the start of the run: z.valid is rnd.valid
          // delayed by nstage of them
          var cycle = 0
          val zs = scala.collection.mutable.ArrayBuffer.empty[Double]

          // the pipeline registers advance while rnd.valid is high; it stays
          // high for the whole run but for the hold check at the end
          c.io.rnd.valid.poke(true.B)

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

          // check if rnd==0 can be handled
          {
            c.io.rnd.bits.poke(0.U(rndW.W))
            for(i <- 0 until nstage) {
              assert(c.io.z.valid.peek().litToBoolean == (cycle >= nstage),
                f"z.valid is ${c.io.z.valid.peek().litToBoolean}%b after $cycle%d " +
                f"cycle(s) of a $nstage%d stage pipeline")
              c.clock.step(1)
              cycle += 1
            }
            c.io.z.bits.expect(0.U)
            c.io.z.valid.expect(true.B)
          }

          // generate random numbers
          for(i <- 1 to n+nstage) {

            c.io.rnd.bits.poke(generateRandomUInt(rndW).U(rndW.W))
            val zi = c.io.z.bits.peek().litValue
            assert(c.io.z.valid.peek().litToBoolean == (cycle >= nstage),
              f"z.valid is ${c.io.z.valid.peek().litToBoolean}%b after $cycle%d " +
              f"cycle(s) of a $nstage%d stage pipeline")

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
            cycle += 1
          }

          // With valid low the pipeline registers do not advance, so the
          // output holds whatever the datapath carried, whatever the input does.
          if (checkHold) {
            assert(nstage > 0, "the hold check needs at least one pipeline register")
            // fill every stage with the same input, so draining the pipeline
            // does not move the output by itself
            c.io.rnd.bits.poke(generateRandomUInt(rndW).U(rndW.W))
            c.clock.step(nstage)
            val zHold = c.io.z.bits.peek().litValue

            c.io.rnd.valid.poke(false.B)
            c.io.rnd.bits.poke(generateRandomUInt(rndW).U(rndW.W))
            for(i <- 1 to nstage+2) {
              c.clock.step(1)
              assert(c.io.z.bits.peek().litValue == zHold,
                f"z.bits moved while rnd.valid was low after $i%d cycle(s)")
            }
            // the low valid has travelled through the whole chain by now
            assert(!c.io.z.valid.peek().litToBoolean,
              "z.valid is still high after rnd.valid was low for the whole latency")
            c.io.rnd.valid.poke(true.B)
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

  runtest( 32, RealSpec.BFloat16Spec, PipelineStageConfig.none)
  runtest( 32, RealSpec.Float32Spec,  PipelineStageConfig.none)
  runtest( 64, RealSpec.Float32Spec,  PipelineStageConfig.none)
  runtest(128, RealSpec.Float32Spec,  PipelineStageConfig.none)

  // pipelined: one cut chosen by the module, two cuts chosen by the module,
  // and one register at every internal boundary
  runtest( 64, RealSpec.Float32Spec,  PipelineStageConfig.default(1))
  runtest( 64, RealSpec.Float32Spec,  PipelineStageConfig.default(2))
  // this configuration also carries the hold check of `valid`
  runtest( 64, RealSpec.Float32Spec,
    PipelineStageConfig.specify("norm" -> 1, "round" -> 1), checkHold = true)
}

