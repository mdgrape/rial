package rial.tests

import org.scalatest.funspec.AnyFunSpec
import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import rial.rng._
import rial.arith.RealSpec
import rial.util.PipelineStageConfig

import scala.util.Random

// Bit-exact match between the GenRandomFloat* hardware and GenRandomFloatSim,
// driving both with the same random input bits.
class GenRandomFloatSimMatchTest extends AnyFunSpec with ChiselSim with EmitVcd {
  addOption(CliOption.int("n", "number of random inputs (default 10000)"))

  it("GenRandomFloat12 HW == GenRandomFloatSim.genFloat12 (bit-exact)") {
    val n    = getOption[Int]("n").getOrElse(10000)
    val spec = RealSpec.Float32Spec
    simulate(new GenRandomFloat12(spec, PipelineStageConfig.none)) { c =>
      val rng = new Random(20260526)
      c.io.rnd.valid.poke(true.B)
      for (_ <- 0 until n) {
        val rnd = BigInt(spec.manW, rng)            // manW-bit mantissa input
        c.io.rnd.bits.poke(rnd.U(spec.manW.W))
        val hw  = c.io.z.bits.peek().litValue
        val sim = GenRandomFloatSim.genFloat12(rnd, spec).value.toBigInt
        assert(hw == sim, f"rnd=0x$rnd%x hw=0x$hw%08x sim=0x$sim%08x")
      }
    }
  }

  it("GenRandomFloat01 HW == GenRandomFloatSim.genFloat01 (bit-exact)") {
    val n    = getOption[Int]("n").getOrElse(10000)
    val spec = RealSpec.Float32Spec
    val rndW = 32
    simulate(new GenRandomFloat01(rndW, spec, PipelineStageConfig.none)) { c =>
      val rng = new Random(20260527)
      c.io.rnd.valid.poke(true.B)
      for (_ <- 0 until n) {
        val rnd = BigInt(rndW, rng)
        c.io.rnd.bits.poke(rnd.U(rndW.W))
        val hw  = c.io.z.bits.peek().litValue
        val sim = GenRandomFloatSim.genFloat01(rnd, rndW, spec).value.toBigInt
        assert(hw == sim, f"rnd=0x$rnd%08x hw=0x$hw%08x sim=0x$sim%08x")
      }
    }
  }
}
