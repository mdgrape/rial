package rial.tests

import org.scalatest.funspec.AnyFunSpec
import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import rial.rng._
import rial.arith.RealSpec

import scala.util.Random
import scala.language.reflectiveCalls

// Bit-exact match between the IrwinHall hardware and IrwinHallSim, driving both
// with the same 12 random words. The pipeline has nStage latency, so the output
// at cycle t corresponds to the inputs poked at cycle t - nStage.
class IrwinHallSimMatchTest extends AnyFunSpec with ChiselSim with EmitVcd {
  addOption(CliOption.int("n", "number of random inputs (default 10000)"))

  it("IrwinHall HW == IrwinHallSim.irwinHall (bit-exact)") {
    val n    = getOption[Int]("n").getOrElse(10000)
    val spec = RealSpec.Float32Spec
    val rndW = 32
    simulate(new IrwinHall(spec, rndW)) { c =>
      val rng    = new Random(20260530)
      val nStage = c.nStage
      val inputs = scala.collection.mutable.ArrayBuffer.empty[Seq[BigInt]]

      c.io.en.poke(true.B)
      c.io.input.valid.poke(true.B)
      for (t <- 0 until n + nStage) {
        val rnds = Seq.fill(12)(BigInt(rndW, rng))
        inputs += rnds
        for (i <- 0 until 12) c.io.input.rnds(i).poke(rnds(i).U(rndW.W))

        if (t >= nStage) {
          assert(c.io.output.valid.peek().litValue == 1)
          val hw  = c.io.output.rnd.peek().litValue
          val sim = IrwinHallSim.irwinHall(inputs(t - nStage), rndW, spec).value.toBigInt
          assert(hw == sim, f"t=$t hw=0x$hw%08x sim=0x$sim%08x in=${inputs(t - nStage).map(_.toString(16))}")
        }
        c.clock.step(1)
      }
    }
  }
}
