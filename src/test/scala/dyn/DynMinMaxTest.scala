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

class DynMinMaxTest extends AnyFunSpec with ChiselSim with EmitVcd {

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

  def minmaxSoft(
    a: BigInt, b: BigInt, width: Int, elemW: Int, signed: Boolean
  ): (BigInt, BigInt) = {
    val ae = toElems(a, width, elemW, signed)
    val be = toElems(b, width, elemW, signed)
    val mn = ae.zip(be).map{ case(x, y) => x.min(y) }
    val mx = ae.zip(be).map{ case(x, y) => x.max(y) }
    (fromElems(mn, elemW), fromElems(mx, elemW))
  }

  def runtest(laneW: Int, nLanes: Int, signed: Boolean,
    stage: PipelineStageConfig, seed: Int) = {
    it(f"DynMinMax(laneW=${laneW}, nLanes=${nLanes}, signed=${signed}, " +
       f"stage=[${stage.getString}])") {
      val n = getOption[Int]("n").getOrElse(1000)
      val width = laneW * nLanes
      val elemWs = Seq.tabulate(log2Ceil(nLanes)+1)(m => laneW << m)

      simulate(new DynMinMax(laneW, nLanes, signed, stage)) {
        c => {
          val rng = new Random(seed)
          // the pipeline registers advance while in.valid is high, which it is
          // for the whole run
          c.io.in.valid.poke(true.B)
          val latency = c.getStage
          assert(latency == stage.resolve(DynMinMax.stagePoints).total)
          // clock edges since the start of the run: out.valid is in.valid
          // delayed by latency of them
          var cycle = 0

          for(elemW <- elemWs) {
            val mode = DynLane.modeOf(elemW, laneW)
            val q = new Queue[(BigInt, BigInt, BigInt, BigInt)]

            for(i <- 0 until n + latency) {
              // make equal elements appear often enough
              val ain = BigInt(width, rng)
              val bin = if (rng.nextInt(4) == 0) { ain } else { BigInt(width, rng) }

              val (mnref, mxref) = minmaxSoft(ain, bin, width, elemW, signed)
              q += ((ain, bin, mnref, mxref))

              c.io.in.bits.a   .poke(ain.U(width.W))
              c.io.in.bits.b   .poke(bin.U(width.W))
              c.io.in.bits.mode.poke(mode.U)

              val mnout = c.io.out.bits.min.peek().litValue
              val mxout = c.io.out.bits.max.peek().litValue
              assert(c.io.out.valid.peek().litToBoolean == (cycle >= latency),
                f"out.valid is ${c.io.out.valid.peek().litToBoolean}%b after " +
                f"$cycle%d cycle(s) of a $latency%d stage pipeline")

              c.clock.step(1)
              cycle += 1

              if(latency <= i) {
                val (x, y, mn, mx) = q.dequeue()
                assert(mnout == mn,
                  f"elemW=${elemW}: min(${x}%x, ${y}%x) = ${mnout}%x, but ref = ${mn}%x")
                assert(mxout == mx,
                  f"elemW=${elemW}: max(${x}%x, ${y}%x) = ${mxout}%x, but ref = ${mx}%x")
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

  describe("DynMinMax compared to software impl") {
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
  //  - default(n) : the registers are placed by DynMinMax itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("DynMinMax with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.specify("in" -> 1, "cmp" -> 1, "out" -> 1))) {
      runtest(8, 4, true, stage, 123456789)
    }
  }
}
