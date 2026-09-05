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

class DynCompareTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))

  // split x into elements of elemW bits. [0] is the LSB side.
  def toElems(x: BigInt, width: Int, elemW: Int, signed: Boolean): Seq[BigInt] = {
    Seq.tabulate(width / elemW)(i => {
      val u = (x >> (elemW * i)) & ((BigInt(1) << elemW) - 1)
      if (signed && u.testBit(elemW-1)) { u - (BigInt(1) << elemW) } else { u }
    })
  }

  // expected per-lane flags. the result of an element is broadcast to all the
  // lanes it occupies.
  def compareSoft(
    a: BigInt, b: BigInt, width: Int, laneW: Int, elemW: Int, signed: Boolean
  ): Seq[(Boolean, Boolean, Boolean)] = {
    val ae = toElems(a, width, elemW, signed)
    val be = toElems(b, width, elemW, signed)
    Seq.tabulate(width / laneW)(i => {
      val e = i / (elemW / laneW)
      (ae(e) == be(e), ae(e) < be(e), ae(e) > be(e))
    })
  }

  def runtest(laneW: Int, nLanes: Int, signed: Boolean,
    stage: PipelineStageConfig, seed: Int) = {
    it(f"DynCompare(laneW=${laneW}, nLanes=${nLanes}, signed=${signed}, " +
       f"stage=[${stage.getString}])") {
      val n = getOption[Int]("n").getOrElse(1000)
      val width = laneW * nLanes
      val elemWs = Seq.tabulate(log2Ceil(nLanes)+1)(m => laneW << m)

      simulate(new DynCompare(laneW, nLanes, signed, stage)) {
        c => {
          val rng = new Random(seed)
          // the pipeline registers advance while in.valid is high, which it is
          // for the whole run
          c.io.in.valid.poke(true.B)
          val latency = c.getStage
          assert(latency == stage.resolve(DynCompare.stagePoints).total)
          // clock edges since the start of the run: out.valid is in.valid
          // delayed by latency of them
          var cycle = 0

          for(elemW <- elemWs) {
            val mode = DynLane.modeOf(elemW, laneW)
            val q = new Queue[(BigInt, BigInt, Seq[(Boolean, Boolean, Boolean)])]

            for(i <- 0 until n + latency) {
              // make equal elements appear often enough
              val ain = BigInt(width, rng)
              val bin = if (rng.nextInt(4) == 0) { ain } else { BigInt(width, rng) }

              q += ((ain, bin, compareSoft(ain, bin, width, laneW, elemW, signed)))

              c.io.in.bits.a   .poke(ain.U(width.W))
              c.io.in.bits.b   .poke(bin.U(width.W))
              c.io.in.bits.mode.poke(mode.U)

              val out = Seq.tabulate(nLanes)(l => (
                c.io.out.bits.eq(l).peek().litValue == 1,
                c.io.out.bits.lt(l).peek().litValue == 1,
                c.io.out.bits.gt(l).peek().litValue == 1
              ))
              assert(c.io.out.valid.peek().litToBoolean == (cycle >= latency),
                f"out.valid is ${c.io.out.valid.peek().litToBoolean}%b after " +
                f"$cycle%d cycle(s) of a $latency%d stage pipeline")

              c.clock.step(1)
              cycle += 1

              if(latency <= i) {
                val (x, y, ref) = q.dequeue()
                for(l <- 0 until nLanes) {
                  assert(out(l) == ref(l),
                    f"elemW=${elemW}, lane=${l}: cmp(${x}%x, ${y}%x) = ${out(l)}, but ref = ${ref(l)}")
                }
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

  describe("DynCompare compared to software impl") {
    // int8x4 / int16x2 / int32
    runtest(8, 4, true,  PipelineStageConfig.none, 123456789)
    runtest(8, 4, false, PipelineStageConfig.none, 123456789)

    // other configurations
    runtest( 4, 8, true,  PipelineStageConfig.none, 123456789)
    runtest( 4, 8, false, PipelineStageConfig.none, 123456789)
    runtest(16, 4, true,  PipelineStageConfig.none, 123456789)
    runtest( 8, 1, true,  PipelineStageConfig.none, 123456789)
    runtest( 1, 4, true,  PipelineStageConfig.none, 123456789)
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by DynCompare itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("DynCompare with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.specify("in" -> 1, "sub" -> 1, "out" -> 1))) {
      runtest(8, 4, true, stage, 123456789)
    }
  }
}
