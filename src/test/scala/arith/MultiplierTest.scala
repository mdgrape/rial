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
// Multiplier: every structure gives the same product as the Scala one.
//
class MultiplierTest extends AnyFunSpec with ChiselSim with EmitVcd {

  val r = new Random(20260101)

  /** Corner values plus `n` random ones of `w` bits. */
  def operands(w: Int, n: Int, signed: Boolean): Seq[BigInt] = {
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

  def mulTest(xW: Int, yW: Int, impl: MulImpl, stage: PipelineStageConfig,
    signed: Boolean, n: Int) = {

    simulate(new Multiplier(xW, yW, impl, stage, signed)) { c =>
      c.io.in.valid.poke(true.B)
      val nstage = c.getStage
      val outW   = xW + yW
      val mask   = (BigInt(1) << outW) - 1
      val xs = operands(xW, n, signed)
      val ys = operands(yW, n, signed)
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
            f"Multiplier(${xW}, ${yW}, ${impl.tag}, signed=${signed}): " +
            f"x=${xe}%x y=${ye}%x -> ${c.io.out.bits.peek().litValue}%x, expected ${ze}%x")
        }
        c.clock.step(1)
      }
    }
  }

  describe("Test an unsigned product") {
    it("plain, combinational") {
      mulTest(24, 24, MulImpl.plain, PipelineStageConfig.none, false, 20)
    }
    it("split at 17, a register at every point") {
      mulTest(24, 24, MulImpl.split(17),
        PipelineStageConfig.specify("in" -> 1, "pp" -> 1, "out" -> 1), false, 20)
    }
    it("split at 17 of an operand that needs three slices") {
      mulTest(39, 24, MulImpl.split(17),
        PipelineStageConfig.specify("pp" -> 1), false, 12)
    }
    it("split at 17 of the right hand operand, which is the wider one") {
      mulTest(12, 39, MulImpl.split(17),
        PipelineStageConfig.specify("pp" -> 1), false, 12)
    }
    it("split at a width the operands already fit") {
      mulTest(12, 12, MulImpl.split(17), PipelineStageConfig.default(2), false, 20)
    }
  }

  describe("Test a signed product") {
    it("plain, combinational") {
      mulTest(18, 18, MulImpl.plain, PipelineStageConfig.none, true, 20)
    }
    it("split at 17, a register at every point") {
      mulTest(18, 34, MulImpl.split(17),
        PipelineStageConfig.specify("in" -> 1, "pp" -> 1, "out" -> 1), true, 12)
    }
    it("split at 17 of an operand that needs three slices") {
      mulTest(40, 18, MulImpl.split(17),
        PipelineStageConfig.specify("pp" -> 1), true, 10)
    }
    it("split at 9, an odd slicing of both signs") {
      mulTest(20, 13, MulImpl.split(9), PipelineStageConfig.default(2), true, 12)
    }
  }
}
