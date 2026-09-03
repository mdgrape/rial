package rial.tests

import org.scalatest._

import chisel3._
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
import rial.testUtil.RialChiselSim

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

//
// Test GenRandomFloat12 using chi-squared test.
//

class Uniform12Test extends AnyFunSpec
    with RialChiselSim with EmitVcd {


  addOption(CliOption.int("n", "number of test cases (default 1000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000)

  val r = new Random(123456789)

  def generateRandomUInt(width: Int) = {
    val rnd = SafeLong(r.nextLong())
    rnd & maskSL(width)
  }

  private def runtest ( spec : RealSpec, stage: PipelineStageConfig) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it(f"uniform12(x) pipereg $pipeconfig spec ${spec.toStringShort}") {
      simulate(new GenRandomFloat12(spec, stage)) { c =>
        {
          val nstage = stage.total
          val inW = spec.manW
          var zs = scala.collection.mutable.ArraySeq.empty[Double]

          val xmax   = 2.0
          val xmin   = 1.0
          val xrange = xmax - xmin

          // Since we use FP that is represented in the binary numeral system,
          // we need to split the range by a power of 2. Otherwise, the min and
          // max numbers of a bin does not align to the FP precision and the
          // "true" number of FP numbers in a bin differs between each other.
          // This cause significant deviation from the uniform distribution.
          val nbins = 16
          val dx    = xrange / nbins
          val nref  = n.toDouble / nbins

          val ndeg      = nbins - 1
          val threshold = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          for(i <- 1 to n+nstage) {
            c.io.rnd.valid.poke(true.B)
            c.io.rnd.bits.poke(generateRandomUInt(inW).toBigInt.U(inW.W))
            val zi = c.io.z.bits.peek().litValue.toBigInt

            if (i > nstage) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              assert(xmin <= zd.toDouble && zd.toDouble < xmax)

              zs = zs :+ zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx

            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

            val term = (nsamples - nref) * (nsamples - nref) / nref
            chi2 += term
          }
          chi2 /= ndeg

          if(threshold < chi2) {
            println(f"-----------------------------------------------")
            for(i <- 0 until nbins) {
              val minRange = xmin +  i    * dx
              val maxRange = xmin + (i+1) * dx
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))
              println(f"n in [${minRange}%6f, ${maxRange}%6f) = ${nsamples} should be ${nref}")
            }
            println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")
          }

          assert(chi2 < threshold)
        }
      }
    }
  }

  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none)
  runtest(RealSpec.Float32Spec,  PipelineStageConfig.none)
  runtest(RealSpec.Float64Spec,  PipelineStageConfig.none)
}


//
// Test GenRandomFloat01 using chi-squared test.
//

class Uniform01Test extends AnyFunSpec
    with RialChiselSim with EmitVcd {


  addOption(CliOption.int("n", "number of test cases (default 1000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000)

  val r = new Random(123456789)

  def generateRandomUInt(width: Int) = {
    val rnd = SafeLong(r.nextLong())
    rnd & maskSL(width)
  }

  private def runtest ( spec : RealSpec, stage: PipelineStageConfig) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it(f"uniform01(x) pipereg $pipeconfig spec ${spec.toStringShort}") {
      simulate(new GenRandomFloat01(64, spec, stage)) { c =>
        {
          val nstage = stage.total
          val inW = 64
          var zs = scala.collection.mutable.ArraySeq.empty[Double]

          val xmax   = 1.0
          val xmin   = 0.0
          val xrange = xmax - xmin

          // Since we use FP that is represented in the binary numeral system,
          // we need to split the range by a power of 2. Otherwise, the min and
          // max numbers of a bin does not align to the FP precision and the
          // "true" number of FP numbers in a bin differs between each other.
          // This cause significant deviation from the uniform distribution.
          val nbins = 16
          val dx    = xrange / nbins
          val nref  = n.toDouble / nbins

          val ndeg      = nbins - 1
          val threshold = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          for(i <- 1 to n+nstage) {
            c.io.rnd.valid.poke(true.B)
            c.io.rnd.bits.poke(generateRandomUInt(inW).toBigInt.U(inW.W))
            val zi = c.io.z.bits.peek().litValue.toBigInt

            if (i > nstage) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              assert(xmin <= zd.toDouble && zd.toDouble < xmax)

              zs = zs :+ zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx

            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

            val term = (nsamples - nref) * (nsamples - nref) / nref
            chi2 += term
          }
          chi2 /= ndeg

          if(threshold < chi2) {
            println(f"-----------------------------------------------")
            for(i <- 0 until nbins) {
              val minRange = xmin +  i    * dx
              val maxRange = xmin + (i+1) * dx
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))
              println(f"n in [${minRange}%6f, ${maxRange}%6f) = ${nsamples} should be ${nref}")
            }
            println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")
          }

          assert(chi2 < threshold)
        }
      }
    }
  }

  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none)
  runtest(RealSpec.Float32Spec,  PipelineStageConfig.none)
  runtest(RealSpec.Float64Spec,  PipelineStageConfig.none)
}
