package rial.tests

import org.scalatest._
import org.scalatest.funspec.AnyFunSpec

import rial.dyn._
import rial.util.PipelineStageConfig

import chisel3._
import chisel3.util._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import scala.collection.mutable.Queue
import scala.util.Random
import scala.math._

class DynMulTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))

  // split x into elements of elemW bits. [0] is the LSB side.
  def toElems(x: BigInt, width: Int, elemW: Int, signed: Boolean): Seq[BigInt] = {
    Seq.tabulate(width / elemW)(i => {
      val u = (x >> (elemW * i)) & ((BigInt(1) << elemW) - 1)
      if (signed && u.testBit(elemW-1)) { u - (BigInt(1) << elemW) } else { u }
    })
  }

  def fromElems(xs: Seq[BigInt], elemW: Int): BigInt = {
    val mask = (BigInt(1) << elemW) - 1
    xs.zipWithIndex.foldLeft(BigInt(0))( (acc, x) => acc | ((x._1 & mask) << (elemW * x._2)) )
  }

  // the product of each element, twice as wide as the element and at the
  // element's own window of the output
  def mulSoft(
    a: BigInt, b: BigInt, width: Int, elemW: Int, signed: Boolean
  ): BigInt = {
    val ae = toElems(a, width, elemW, signed)
    val be = toElems(b, width, elemW, signed)
    fromElems(ae.zip(be).map{ case(x, y) => x * y }, 2 * elemW)
  }

  // an operand whose every element is one of the values the sign handling and
  // the recoding of the multiplier are most sensitive to: zero, all ones, the
  // most negative value, one, the most positive value and the two alternating
  // bit patterns, which give the longest runs of nonzero digits
  def special(width: Int, elemW: Int, rng: Random): BigInt = {
    val alternating = Seq.tabulate(elemW)(i => BigInt(i % 2) << i).sum
    val values = Seq(
      BigInt(0), (BigInt(1) << elemW) - 1, BigInt(1) << (elemW-1), BigInt(1),
      (BigInt(1) << (elemW-1)) - 1, alternating, alternating >> 1)
    fromElems(Seq.fill(width / elemW)(values(rng.nextInt(values.length))), elemW)
  }

  def operand(width: Int, elemW: Int, rng: Random): BigInt = {
    if (rng.nextInt(4) == 0) {
      special(width, elemW, rng)
    } else {
      BigInt(width, rng)
    }
  }

  def runtest(laneW: Int, nLanes: Int, signed: Boolean,
    stage: PipelineStageConfig, seed: Int) = {
    it(f"DynMul(laneW=${laneW}, nLanes=${nLanes}, signed=${signed}, " +
       f"stage=[${stage.getString}])") {
      val n = getOption[Int]("n").getOrElse(1000)
      val width = laneW * nLanes
      val elemWs = Seq.tabulate(log2Ceil(nLanes)+1)(m => laneW << m)

      simulate(new DynMul(laneW, nLanes, signed, stage)) {
        c => {
          val rng = new Random(seed)
          // the pipeline registers advance while in.valid is high, which it is
          // for the whole run
          c.io.in.valid.poke(true.B)
          val latency = c.getStage
          assert(latency == stage.resolve(DynMul.stagePoints).total)
          // clock edges since the start of the run: out.valid is in.valid
          // delayed by latency of them
          var cycle = 0

          for(elemW <- elemWs) {
            val mode = DynLane.modeOf(elemW, laneW)
            val q = new Queue[(BigInt, BigInt, BigInt)]

            for(i <- 0 until n + latency) {
              val ain = operand(width, elemW, rng)
              val bin = operand(width, elemW, rng)

              q += ((ain, bin, mulSoft(ain, bin, width, elemW, signed)))

              c.io.in.bits.a   .poke(ain.U(width.W))
              c.io.in.bits.b   .poke(bin.U(width.W))
              c.io.in.bits.mode.poke(mode.U)

              val out = c.io.out.bits.prod.peek().litValue
              assert(c.io.out.valid.peek().litToBoolean == (cycle >= latency),
                f"out.valid is ${c.io.out.valid.peek().litToBoolean}%b after " +
                f"$cycle%d cycle(s) of a $latency%d stage pipeline")

              c.clock.step(1)
              cycle += 1

              if(latency <= i) {
                val (x, y, ref) = q.dequeue()
                assert(out == ref,
                  f"elemW=${elemW}: ${x}%x * ${y}%x = ${out}%x, but ref = ${ref}%x")
              }
            }
            // the elements of the queue that are still in flight belong to the
            // next mode, which pokes its own inputs from now on
            q.clear()
          }
        }
      }
    }
  }

  describe("DynMul") {
    // int8x4 / int16x2 / int32
    runtest(8, 4, true,  PipelineStageConfig.none, 123456789)
    runtest(8, 4, false, PipelineStageConfig.none, 123456789)

    // other configurations
    runtest( 4, 8, true,  PipelineStageConfig.none, 123456789)
    runtest(16, 2, true,  PipelineStageConfig.none, 123456789)
    runtest( 8, 1, true,  PipelineStageConfig.none, 123456789)
    runtest( 8, 1, false, PipelineStageConfig.none, 123456789)
    runtest( 2, 4, true,  PipelineStageConfig.none, 123456789)
    runtest( 2, 4, false, PipelineStageConfig.none, 123456789)
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by DynMul itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("DynMul with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.specify("in" -> 1, "reduce" -> 1, "out" -> 1))) {
      runtest(8, 4, true, stage, 123456789)
    }
  }
}
