package rial.tests

import org.scalatest.funspec.AnyFunSpec

import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd

import rial.arith._
import rial.util.PipelineStageConfig

import scala.collection.mutable.Queue
import scala.util.Random

//
// SplitMultiplier: a product cut into slices is the product the Scala one is.
//
class SplitMultiplierTest extends AnyFunSpec with ChiselSim with EmitVcd {

  val r = new Random(20260101)

  /** Corner values plus `n` random ones of `w` bits. */
  def operands(w: Int, n: Int): Seq[BigInt] = {
    val mask = (BigInt(1) << w) - 1
    val corner = Seq(BigInt(0), BigInt(1), mask, mask - 1,
      BigInt(1) << (w - 1), (BigInt(1) << (w - 1)) - 1)
    corner ++ (1 to n).map(_ => BigInt(w, r))
  }

  /** Two's complement reading of `v`, `w` bits wide. */
  def signedOf(v: BigInt, w: Int): BigInt = {
    if (v.testBit(w - 1)) {
      v - (BigInt(1) << w)
    } else {
      v
    }
  }

  def mulTest(xW: Int, yW: Int, sliceW: Int, stage: PipelineStageConfig,
    signed: Boolean, n: Int) = {

    simulate(new SplitMultiplier(xW, yW, sliceW, stage, signed)) { c =>
      c.io.in.valid.poke(true.B)
      val nstage = c.getStage
      val outW   = xW + yW
      val mask   = (BigInt(1) << outW) - 1
      val xs = operands(xW, n)
      val ys = operands(yW, n)
      val q  = new Queue[(BigInt, BigInt, BigInt)]
      val cases = for (x <- xs; y <- ys) yield (x, y)
      for ((x, y) <- cases ++ Seq.fill(nstage)((BigInt(0), BigInt(0)))) {
        val prod = if (signed) {
          (signedOf(x, xW) * signedOf(y, yW)) & mask
        } else {
          x * y
        }
        c.io.in.bits.x.poke(x.U(xW.W))
        c.io.in.bits.y.poke(y.U(yW.W))
        q.enqueue((x, y, prod))
        if (q.length > nstage) {
          val (xe, ye, ze) = q.dequeue()
          c.io.out.valid.expect(true.B)
          assert(c.io.out.bits.peek().litValue == ze,
            f"SplitMultiplier(${xW}, ${yW}, ${sliceW}, signed=${signed}): " +
            f"x=${xe}%x y=${ye}%x -> ${c.io.out.bits.peek().litValue}%x, " +
            f"expected ${ze}%x")
        }
        c.clock.step(1)
      }
    }
  }

  describe("Test an unsigned product") {
    it("split at 17, a register at every point") {
      mulTest(24, 24, 17,
        PipelineStageConfig.specify("in" -> 1, "pp" -> 1, "out" -> 1), false, 20)
    }
    it("split at 17 of an operand that needs three slices") {
      mulTest(39, 24, 17, PipelineStageConfig.specify("pp" -> 1), false, 12)
    }
    it("three slices with a register between the hops of the chain") {
      mulTest(39, 24, 17,
        PipelineStageConfig.specify("pp" -> 1, "chain1" -> 1, "out" -> 1),
        false, 12)
    }
    it("four slices with both chain points registered") {
      mulTest(56, 12, 17,
        PipelineStageConfig.specify("in" -> 1, "pp" -> 1, "chain1" -> 1,
          "chain2" -> 1, "out" -> 1), false, 8)
    }
    it("split at 17 of the right hand operand, which is the wider one") {
      mulTest(12, 39, 17, PipelineStageConfig.specify("pp" -> 1), false, 12)
    }
    it("split at a width the operands already fit") {
      mulTest(12, 12, 17, PipelineStageConfig.default(2), false, 20)
    }
  }

  describe("Test a signed product") {
    it("split at 17, a register at every point") {
      mulTest(18, 34, 17,
        PipelineStageConfig.specify("in" -> 1, "pp" -> 1, "out" -> 1), true, 12)
    }
    it("split at 17 of an operand that needs three slices") {
      mulTest(40, 18, 17, PipelineStageConfig.specify("pp" -> 1), true, 10)
    }
    it("split at 9, an odd slicing of both signs") {
      mulTest(20, 13, 9, PipelineStageConfig.default(2), true, 12)
    }
  }

  describe("Test the points a product offers") {
    it("chains only what an unsigned product of more than two slices sums") {
      assert(SplitMultiplier.chainPoints(39, 24, 17) == 1)
      assert(SplitMultiplier.chainPoints(56, 12, 17) == 2)
      assert(SplitMultiplier.chainPoints(24, 24, 17) == 0)
      assert(SplitMultiplier.chainPoints(39, 24, 17, signed = true) == 0)
    }
    it("names the chain points between the partial products and the output") {
      assert(SplitMultiplier.stagePoints().map(_.name) == Seq("pp", "out"))
      assert(SplitMultiplier.stagePoints(2).map(_.name) ==
        Seq("pp", "chain1", "chain2", "out"))
    }
  }
}
