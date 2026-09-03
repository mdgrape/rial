package rial.tests

import org.scalatest._

import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

import chisel3._
import chisel3.experimental.BundleLiterals._
import chisel3.util._

import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import org.scalatest.funspec.AnyFunSpec

import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig
import rial.rng._
import rial.testUtil.RialChiselSim

class SquaresTest extends AnyFunSpec
  with RialChiselSim with EmitVcd {


  addOption(CliOption.int("n", "number of test cases (default 1000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000)

  // verify the HW against the software model (SquaresSim) used by RngSim.
  def squares(cnt: BigInt, key: BigInt): BigInt = SquaresSim.squares32(cnt, key)

  def runTest(key: BigInt) = {
    simulate(new Squares32) {
      c => {
        c.io.en.poke(true.B)

        c.io.input.valid.poke(true.B)
        c.io.input.key.poke(key.U)

        for (i <- 0 to n-1) {
          c.io.input.count.poke(i.U)

          if(i >= 4) {
            assert(c.io.output.valid.peek().litValue == 1)
            assert(c.io.output.rand.peek().litValue == squares(BigInt(i-4), key))
          } else {
            assert(c.io.output.valid.peek().litValue == 0)
          }

          c.clock.step(1)
        }
      }
    }
  }

  it("Squares should be equal to software impl") {
    runTest(BigInt("2467cb532b5ce8d1", 16))
  }
}

class Squares32GeneratorTest extends AnyFunSpec
  with RialChiselSim with EmitVcd {


  addOption(CliOption.int("n", "number of test cases (default 1000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000)

  def squares(cnt: BigInt, key: BigInt): BigInt = SquaresSim.squares32(cnt, key)

  def runTest(key: BigInt) = {
    simulate(new Squares32Generator) {
      c => {
        c.io.initialized.expect(false.B)

        c.io.init.en   .poke(true.B)
        c.io.init.key  .poke(key.U)
        c.io.init.count.poke(  0.U)

        c.io.rand.ready.poke(false.B)

        c.clock.step(1) // --------------------------------------------

        c.io.initialized.expect(true.B)

        c.io.init.en   .poke(false.B)
        c.io.init.key  .poke(0.U)
        c.io.init.count.poke(0.U)

        c.io.rand.ready.poke(true.B)
        for (i <- 0 to n-1) {
          c.io.rand.valid.expect( (i >= c.nStage).B )

          if(c.io.rand.valid.peek().litValue == 1) {
            c.io.rand.bits.expect(squares(BigInt(i-c.nStage), key).U)
          }
          c.clock.step(1)
        }
      }
    }
  }

  it("Squares should be equal to software impl") {
    runTest(BigInt("2467cb532b5ce8d1", 16))
  }
}
