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

class DynAddTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))

  // split x into elements of elemW bits. [0] is the LSB side.
  def toElems(x: BigInt, width: Int, elemW: Int): Seq[BigInt] = {
    Seq.tabulate(width / elemW)(i =>
      (x >> (elemW * i)) & ((BigInt(1) << elemW) - 1))
  }

  def fromElems(xs: Seq[BigInt], elemW: Int): BigInt = {
    val mask = (BigInt(1) << elemW) - 1
    xs.zipWithIndex.foldLeft(BigInt(0))( (acc, x) => acc | ((x._1 & mask) << (elemW * x._2)) )
  }

  def addSoft(a: BigInt, b: BigInt, width: Int, elemW: Int): BigInt = {
    val ae = toElems(a, width, elemW)
    val be = toElems(b, width, elemW)
    fromElems(ae.zip(be).map{ case(x, y) => x + y }, elemW)
  }

  // an operand whose every element is one of the values the carries are most
  // sensitive to: zero, all ones, the most negative value and one
  def special(width: Int, elemW: Int, rng: Random): BigInt = {
    val values = Seq(
      BigInt(0), (BigInt(1) << elemW) - 1, BigInt(1) << (elemW-1), BigInt(1))
    fromElems(Seq.fill(width / elemW)(values(rng.nextInt(values.length))), elemW)
  }

  def operand(width: Int, elemW: Int, rng: Random): BigInt = {
    if (rng.nextInt(4) == 0) {
      special(width, elemW, rng)
    } else {
      BigInt(width, rng)
    }
  }

  def runtest(laneW: Int, nLanes: Int, stage: PipelineStageConfig, seed: Int) = {
    it(f"DynAdd(laneW=${laneW}, nLanes=${nLanes}, stage=[${stage.getString}])") {
      val n = getOption[Int]("n").getOrElse(1000)
      val width = laneW * nLanes
      val elemWs = Seq.tabulate(log2Ceil(nLanes)+1)(m => laneW << m)

      simulate(new DynAdd(laneW, nLanes, stage)) {
        c => {
          val rng = new Random(seed)
          // the pipeline registers advance while in.valid is high, which it is
          // for the whole run
          c.io.in.valid.poke(true.B)
          val latency = c.getStage
          assert(latency == stage.resolve(DynAdd.stagePoints).total)
          // clock edges since the start of the run: out.valid is in.valid
          // delayed by latency of them
          var cycle = 0

          for(elemW <- elemWs) {
            val mode = DynLane.modeOf(elemW, laneW)
            val q = new Queue[(BigInt, BigInt, BigInt)]

            for(i <- 0 until n + latency) {
              val ain = operand(width, elemW, rng)
              val bin = operand(width, elemW, rng)

              q += ((ain, bin, addSoft(ain, bin, width, elemW)))

              c.io.in.bits.a   .poke(ain.U(width.W))
              c.io.in.bits.b   .poke(bin.U(width.W))
              c.io.in.bits.mode.poke(mode.U)

              val out = c.io.out.bits.sum.peek().litValue
              assert(c.io.out.valid.peek().litToBoolean == (cycle >= latency),
                f"out.valid is ${c.io.out.valid.peek().litToBoolean}%b after " +
                f"$cycle%d cycle(s) of a $latency%d stage pipeline")

              c.clock.step(1)
              cycle += 1

              if(latency <= i) {
                val (x, y, ref) = q.dequeue()
                assert(out == ref,
                  f"elemW=${elemW}: ${x}%x + ${y}%x = ${out}%x, but ref = ${ref}%x")
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

  describe("DynAdd compared to software impl") {
    // int8x4 / int16x2 / int32
    runtest(8, 4, PipelineStageConfig.none, 123456789)

    // other configurations
    runtest( 4, 8, PipelineStageConfig.none, 123456789)
    runtest(16, 2, PipelineStageConfig.none, 123456789)
    runtest( 8, 1, PipelineStageConfig.none, 123456789)
    runtest( 1, 4, PipelineStageConfig.none, 123456789)
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by DynAdd itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("DynAdd with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.specify("in" -> 1, "out" -> 1))) {
      runtest(8, 4, stage, 123456789)
    }
  }
}
