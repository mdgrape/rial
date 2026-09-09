package rial.tests

import org.scalatest.funspec.AnyFunSpec

import chisel3._
import circt.stage.ChiselStage
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
            f"x=${xe}%x y=${ye}%x -> ${c.io.out.bits.peek().litValue}%x, " +
            f"expected ${ze}%x")
        }
        c.clock.step(1)
      }
    }
  }

  /** `(x * y + c) >> shift` as the module reads it, `outW` bits wide. */
  def mulAddOf(x: BigInt, y: BigInt, cv: BigInt, xW: Int, yW: Int, cW: Int,
    shift: Int, outW: Int, signed: Boolean): BigInt = {

    val full = if (signed) {
      signedOf(x, xW) * signedOf(y, yW) + signedOf(cv, cW)
    } else {
      x * y + cv
    }
    (full >> shift) & ((BigInt(1) << outW) - 1)
  }

  def mulAddTest(xW: Int, yW: Int, cW: Int, shift: Int, impl: MulImpl,
    stage: PipelineStageConfig, signed: Boolean, n: Int) = {

    simulate(new Multiplier(xW, yW, impl, stage, signed, Some(cW),
      Some(shift))) { dut =>

      dut.io.in.valid.poke(true.B)
      val nstage = dut.getStage
      val outW   = Multiplier.outWidth(xW, yW, Some(cW), shift)
      val xs = operands(xW, n, signed)
      val ys = operands(yW, n, signed)
      val cs = operands(cW, n, signed)
      val q  = new Queue[(BigInt, BigInt, BigInt, BigInt)]
      val cases = for (x <- xs; y <- ys; cv <- cs) yield (x, y, cv)
      val flush = Seq.fill(nstage)((BigInt(0), BigInt(0), BigInt(0)))
      for ((x, y, cv) <- cases ++ flush) {
        val want = mulAddOf(x, y, cv, xW, yW, cW, shift, outW, signed)
        dut.io.in.bits.x.poke(x.U(xW.W))
        dut.io.in.bits.y.poke(y.U(yW.W))
        dut.io.in.bits.c.get.poke(cv.U(cW.W))
        q.enqueue((x, y, cv, want))
        if (q.length > nstage) {
          val (xe, ye, ce, ze) = q.dequeue()
          dut.io.out.valid.expect(true.B)
          assert(dut.io.out.bits.peek().litValue == ze,
            f"Multiplier(${xW}, ${yW}, ${impl.tag}, signed=${signed}, " +
            f"cW=${cW}, shift=${shift}): x=${xe}%x y=${ye}%x c=${ce}%x -> " +
            f"${dut.io.out.bits.peek().litValue}%x, expected ${ze}%x")
        }
        dut.clock.step(1)
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
    it("three slices with a register between the hops of the chain") {
      mulTest(39, 24, MulImpl.split(17),
        PipelineStageConfig.specify("pp" -> 1, "chain1" -> 1, "out" -> 1),
        false, 12)
    }
    it("four slices with both chain points registered") {
      mulTest(56, 12, MulImpl.split(17),
        PipelineStageConfig.specify("in" -> 1, "pp" -> 1, "chain1" -> 1,
          "chain2" -> 1, "out" -> 1), false, 8)
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

  describe("Test a multiply-add, which is the macc structure") {
    it("one multiplication and unsigned, with a shift") {
      mulAddTest(12, 12, 20, 7, MulImpl.macc(17),
        PipelineStageConfig.specify("out" -> 1), false, 3)
    }
    it("one multiplication and signed, with a shift") {
      mulAddTest(12, 12, 20, 7, MulImpl.macc(17),
        PipelineStageConfig.specify("in" -> 1, "out" -> 1), true, 3)
    }
    it("one multiplication and signed, with no shift at all") {
      mulAddTest(12, 12, 26, 0, MulImpl.macc(17), PipelineStageConfig.none,
        true, 3)
    }
    it("unsigned, the addend on the sum") {
      mulAddTest(24, 24, 40, 9, MulImpl.macc(17),
        PipelineStageConfig.specify("pp" -> 1, "out" -> 1), false, 3)
    }
    it("unsigned over three slices") {
      mulAddTest(39, 24, 50, 11, MulImpl.macc(17),
        PipelineStageConfig.specify("pp" -> 1, "chain1" -> 1), false, 2)
    }
    it("signed, the addend on the tiles") {
      mulAddTest(21, 20, 43, 20, MulImpl.macc(17),
        PipelineStageConfig.specify("pp" -> 1, "out" -> 1), true, 3)
    }
    it("signed, the shift below the cut") {
      mulAddTest(21, 20, 43, 5, MulImpl.macc(17),
        PipelineStageConfig.specify("pp" -> 1), true, 3)
    }
    it("signed, no shift, so the tiles meet at full width") {
      mulAddTest(21, 20, 43, 0, MulImpl.macc(17),
        PipelineStageConfig.specify("pp" -> 1, "out" -> 1), true, 3)
    }
    it("signed, an operand of two full slices") {
      mulAddTest(19, 33, 52, 18, MulImpl.macc(17),
        PipelineStageConfig.specify("pp" -> 1), true, 2)
    }
    it("signed, an operand that fits one slice") {
      mulAddTest(19, 12, 30, 18, MulImpl.macc(17),
        PipelineStageConfig.specify("pp" -> 1, "out" -> 1), true, 3)
    }
    it("a slice width of 9 and signed, an odd cut") {
      mulAddTest(14, 13, 26, 13, MulImpl.macc(9),
        PipelineStageConfig.default(2), true, 3)
    }
    it("refuses an operand wider than two slices") {
      assertThrows[IllegalArgumentException] {
        Multiplier.cutBit(19, 40, MulImpl.macc(17), true, 18)
      }
    }
    it("refuses a multiply-add with no addend") {
      assertThrows[IllegalArgumentException] {
        ChiselStage.emitCHIRRTL(new Multiplier(19, 24, MulImpl.macc(17),
          PipelineStageConfig.none, true))
      }
    }
    it("refuses an addend on a structure that has none") {
      for (impl <- Seq(MulImpl.plain, MulImpl.split(17))) {
        assertThrows[IllegalArgumentException] {
          ChiselStage.emitCHIRRTL(new Multiplier(19, 24, impl,
            PipelineStageConfig.none, true, Some(41), Some(18)))
        }
      }
    }
    it("reads its own name") {
      assert(MulImpl.parse("macc:17") == MulImpl.macc(17))
      assert(MulImpl.macc(17).tag == "macc17")
      assert(MulImpl.macc(17).isMacc && MulImpl.macc(17).isSliced)
      assert(Multiplier.zeroAddendW(MulImpl.macc(17)).isDefined)
      assert(Multiplier.zeroAddendW(MulImpl.split(17)).isEmpty)
      assertThrows[IllegalArgumentException] {
        MulImpl.parse("macc")
      }
    }
  }
}
