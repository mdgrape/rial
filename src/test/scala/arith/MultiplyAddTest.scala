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
// MultiplyAdd: (x * y >> shift) + c is the number the Scala one is.
//
class MultiplyAddTest extends AnyFunSpec with ChiselSim with EmitVcd {

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

  /** `(x * y >> shift) + c` as the module reads it, `outW` bits wide. */
  def mulAddOf(x: BigInt, y: BigInt, cv: BigInt, xW: Int, yW: Int, cW: Int,
    shift: Int, outW: Int, signed: Boolean): BigInt = {

    val prod = if (signed) {
      signedOf(x, xW) * signedOf(y, yW)
    } else {
      x * y
    }
    val addend = if (signed) {
      signedOf(cv, cW)
    } else {
      cv
    }
    ((prod >> shift) + addend) & ((BigInt(1) << outW) - 1)
  }

  def mulAddTest(xW: Int, yW: Int, cW: Int, shift: Int, sliceW: Int,
    stage: PipelineStageConfig, signed: Boolean, n: Int,
    small: Boolean = false) = {

    simulate(new MultiplyAdd(xW, yW, cW, shift, sliceW, stage, signed,
      small)) { dut =>

      dut.io.in.valid.poke(true.B)
      val nstage = dut.getStage
      val outW   = MultiplyAdd.outWidth(xW, yW, cW, shift)
      val xs = operands(xW, n)
      val ys = operands(yW, n)
      val cs = operands(cW, n)
      val q  = new Queue[(BigInt, BigInt, BigInt, BigInt)]
      val cases = for (x <- xs; y <- ys; cv <- cs) yield (x, y, cv)
      val flush = Seq.fill(nstage)((BigInt(0), BigInt(0), BigInt(0)))
      for ((x, y, cv) <- cases ++ flush) {
        val want = mulAddOf(x, y, cv, xW, yW, cW, shift, outW, signed)
        dut.io.in.bits.x.poke(x.U(xW.W))
        dut.io.in.bits.y.poke(y.U(yW.W))
        dut.io.in.bits.c.poke(cv.U(cW.W))
        q.enqueue((x, y, cv, want))
        if (q.length > nstage) {
          val (xe, ye, ce, ze) = q.dequeue()
          dut.io.out.valid.expect(true.B)
          assert(dut.io.out.bits.peek().litValue == ze,
            f"MultiplyAdd(${xW}, ${yW}, ${cW}, ${shift}, ${sliceW}, " +
            f"signed=${signed}, small=${small}): " +
            f"x=${xe}%x y=${ye}%x c=${ce}%x -> " +
            f"${dut.io.out.bits.peek().litValue}%x, expected ${ze}%x")
        }
        dut.clock.step(1)
      }
    }
  }

  describe("Test an unsigned multiply-add") {
    it("one multiplication, with a shift") {
      mulAddTest(12, 12, 20, 7, 17,
        PipelineStageConfig.specify("out" -> 1), false, 3)
    }
    it("the addend on the sum") {
      mulAddTest(24, 24, 40, 9, 17,
        PipelineStageConfig.specify("pp" -> 1, "out" -> 1), false, 3)
    }
    it("over three slices") {
      mulAddTest(39, 24, 50, 11, 17,
        PipelineStageConfig.specify("in" -> 1, "pp" -> 1), false, 2)
    }
    it("the addend below the product") {
      mulAddTest(24, 24, 12, 9, 17,
        PipelineStageConfig.specify("pp" -> 1), false, 3)
    }
    it("the small form is the same multiply-add") {
      mulAddTest(39, 24, 50, 11, 17,
        PipelineStageConfig.specify("pp" -> 1), false, 2, true)
    }
  }

  for (small <- Seq(false, true)) {
    describe(s"Test a signed multiply-add, small = $small") {
      it("one multiplication, with a shift") {
        mulAddTest(12, 12, 20, 7, 17,
          PipelineStageConfig.specify("in" -> 1, "out" -> 1), true, 3, small)
      }
      it("one multiplication, with no shift at all") {
        mulAddTest(12, 12, 26, 0, 17, PipelineStageConfig.none, true, 3, small)
      }
      it("one multiplication, the addend below the product") {
        mulAddTest(12, 12, 8, 5, 17,
          PipelineStageConfig.specify("pp" -> 1), true, 3, small)
      }
      it("an operand one bit wider than a slice") {
        // the widths of a step of a fixed point function of the math library,
        // whose addend is wider than the shifted product
        mulAddTest(18, 22, 27, 17, 17,
          PipelineStageConfig.specify("pp" -> 1), true, 2, small)
      }
      it("the addend on the sum") {
        mulAddTest(21, 20, 43, 20, 17,
          PipelineStageConfig.specify("pp" -> 1, "out" -> 1), true, 3, small)
      }
      it("the shift below the slice the operand is cut at") {
        mulAddTest(21, 20, 43, 5, 17,
          PipelineStageConfig.specify("pp" -> 1), true, 3, small)
      }
      it("no shift, so the sum is the result") {
        mulAddTest(21, 20, 43, 0, 17,
          PipelineStageConfig.specify("pp" -> 1, "out" -> 1), true, 3, small)
      }
      it("both operands wider than a slice") {
        // the widths of the last step of the gaussian table
        mulAddTest(23, 23, 23, 22, 17,
          PipelineStageConfig.specify("pp" -> 1, "out" -> 1), true, 2, small)
      }
      it("both operands wider than a slice, a register at every point") {
        mulAddTest(23, 23, 23, 22, 17,
          PipelineStageConfig.specify("in" -> 1, "pp" -> 1, "out" -> 1),
          true, 2, small)
      }
      it("an operand of two full slices") {
        mulAddTest(19, 33, 52, 18, 17,
          PipelineStageConfig.specify("pp" -> 1), true, 2, small)
      }
      it("an operand that fits one slice") {
        mulAddTest(19, 12, 30, 18, 17,
          PipelineStageConfig.specify("pp" -> 1, "out" -> 1), true, 3, small)
      }
      it("a slice width of 9, an odd cut") {
        mulAddTest(14, 13, 26, 13, 9, PipelineStageConfig.default(2), true, 3,
          small)
      }
    }
  }

  describe("Test the width and the names of a multiply-add") {
    it("is one bit above both the shifted product and the addend") {
      assert(MultiplyAdd.outWidth(24, 24, 12, 9) == 24 + 24 - 9 + 1)
      assert(MultiplyAdd.outWidth(24, 24, 40, 9) == 40 + 1)
    }
    it("reads its own name") {
      assert(MulImpl.parse("macc:17") == MulImpl.macc(17))
      assert(MulImpl.macc(17).tag == "macc17")
      assert(MulImpl.macc(17).isMacc && MulImpl.macc(17).isSliced)
      assertThrows[IllegalArgumentException] {
        MulImpl.parse("macc")
      }
    }
    it("reads the names of the small forms") {
      assert(MulImpl.parse("splitSmall:17") == MulImpl.splitSmall(17))
      assert(MulImpl.parse("maccSmall:17") == MulImpl.maccSmall(17))
      assert(MulImpl.splitSmall(17).tag == "splitSmall17")
      assert(MulImpl.splitSmall(17).isSplit && MulImpl.splitSmall(17).isSmall)
      assert(MulImpl.maccSmall(17).isMacc && MulImpl.maccSmall(17).isSliced)
      assert(MulImpl.maccSmall(17).isSmall && !MulImpl.macc(17).isSmall)
      assert(!MulImpl.split(17).isSmall)
      assertThrows[IllegalArgumentException] {
        MulImpl.parse("splitSmall")
      }
    }
  }
}
