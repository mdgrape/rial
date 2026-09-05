package rial.tests

import org.scalatest._
import org.scalatest.funspec.AnyFunSpec

import rial.ecc._
import rial.util.PipelineStageConfig

import chisel3._
import chisel3.util._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import scala.collection.mutable.{Queue, ArrayBuffer}
import scala.util.Random
import scala.math._

class FletchersChecksumTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))

  // Ai = Ai-1 + di
  // Bi = Bi-1 + Ai
  def fletcher16Soft(w: Int, data: BigInt): Int = {

    val ds = Seq.tabulate((w+7) / 8)(i => {
      ((data >> (8 * i)) & 0x00FF).toInt
    })

    var a = 0
    var b = 0
    for(d <- ds) {
      a = a + d
      b = b + a
    }
    val amod = a % 255
    val bmod = b % 255

    (bmod << 8) + amod
  }

  def runtest(dataW: Int, stage: PipelineStageConfig, seed: Int) = {
    it(f"Fletcher16(dataW=${dataW}, stage=[${stage.getString}])") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new Fletcher(16, dataW, stage)) {
        c => {
          val rng = new Random(seed)
          // the pipeline registers advance while in.valid is high, which it is
          // for the whole run
          c.io.in.valid.poke(true.B)
          val latency = c.getStage
          assert(latency == stage.resolve(Fletcher.stagePoints).total)

          val q = new Queue[(BigInt, Int)]
          for(i <- 0 until n + latency) {

            val xin = BigInt(dataW, rng).abs
            val zref = fletcher16Soft(dataW, xin)
            q += ((xin, zref))

            c.io.in.bits.poke(xin.U(dataW.W))
            val zvalid = c.io.out.valid.peek().litToBoolean
            val zout   = c.io.out.bits.peek().litValue

            // println(f"xin=${xin}%x, zout=${zout}%x, zref=${zref}%x")

            // out.valid is in.valid delayed by the latency
            assert(zvalid == (latency <= i),
              f"out.valid is ${zvalid}%b after $i%d cycle(s) of a " +
              f"$latency%d stage pipeline")

            c.clock.step(1)

            if(latency <= i) {
              val (x, z) = q.dequeue()
              assert(z == zout, f"checksum(${x}%x) = ${z}%x, but checksumSoft(${x}%x) = ${zout}%x")
            }
          }
        }
      }
    }
  }

  describe("Fletcher's checksum compared to software impl") {
    for (dataW <- Seq(32, 64, 128)) {
      runtest(dataW, PipelineStageConfig.none, 123456789)
    }
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by Fletcher itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("Fletcher's checksum with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.specify("in" -> 1, "sum" -> 1, "mod" -> 1, "out" -> 1))) {
      runtest(64, stage, 123456789)
    }
  }
}
