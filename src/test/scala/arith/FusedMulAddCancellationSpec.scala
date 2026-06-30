// Minimal regression for the FusedMulAddFPGeneric exact-cancellation defect.
//
// When x*y + z cancels exactly to zero, the fused multiply-add must return +0.
// Before the fix the near path renormalizes a zero significand into a spurious
// tiny normal:
//   BF16  1.0 * 1.0 + (-1.0)  ->  0x3880   (~6.1e-5)   instead of 0x0000
//   FP32  1.0 * 1.0 + (-1.0)  ->  0x48800000 (262144)  instead of 0x00000000
//
// Place at src/test/scala/arith/FusedMulAddCancellationSpec.scala and run:
//   sbt "testOnly rial.tests.FusedMulAddCancellationSpec"
package rial.tests

import chisel3._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import rial.arith.{FusedMulAddFPGeneric, RealSpec, RoundSpec}
import rial.util.PipelineStageConfig

class FusedMulAddCancellationSpec extends AnyFlatSpec with ChiselScalatestTester with Matchers {

  /** Combinational BF16/FP32 FMA wrapper driving raw bit encodings. */
  private class Fma(s: RealSpec) extends Module {
    val io = IO(new Bundle {
      val x = Input(UInt(s.W.W)); val y = Input(UInt(s.W.W))
      val z = Input(UInt(s.W.W)); val w = Output(UInt(s.W.W))
    })
    val u = Module(
      new FusedMulAddFPGeneric(s, s, s, s, RoundSpec.roundToEven, PipelineStageConfig.none)
    )
    u.io.x := io.x; u.io.y := io.y; u.io.z := io.z; io.w := u.io.w
  }

  private def w(c: Fma, x: Long, y: Long, z: Long): BigInt = {
    c.io.x.poke(BigInt(x).U(c.io.x.getWidth.W))
    c.io.y.poke(BigInt(y).U(c.io.y.getWidth.W))
    c.io.z.poke(BigInt(z).U(c.io.z.getWidth.W))
    c.io.w.peek().litValue // combinational (PipelineStageConfig.none)
  }

  behavior of "FusedMulAddFPGeneric exact cancellation (x*y + z == 0 -> +0)"

  it should "return +0 for BF16 1.0*1.0 + (-1.0) (and keep non-cancelling cases correct)" in {
    test(new Fma(RealSpec.BFloat16Spec)).withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
      // The essential bug: exact cancellation must be +0. Pre-fix: 0x3880.
      w(c, 0x3f80, 0x3f80, 0xbf80) shouldBe BigInt(0x0000)
      // A couple more exact-cancellation shapes (residue used to scale with operands).
      w(c, 0x4000, 0x4000, 0xc080) shouldBe BigInt(0x0000) //  2.0 * 2.0 + (-4.0)
      w(c, 0xbf80, 0x3f80, 0x3f80) shouldBe BigInt(0x0000) // -1.0 * 1.0 + ( 1.0)
      // Control: a non-cancelling near case must be untouched by the fix.
      w(c, 0x3e80, 0x4000, 0x3f80) shouldBe BigInt(0x3fc0) // 0.25 * 2.0 + 1.0 = 1.5
    }
  }

  it should "return +0 for FP32 1.0*1.0 + (-1.0)" in {
    test(new Fma(RealSpec.Float32Spec)).withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
      // Pre-fix: 0x48800000 (= 262144).
      w(c, 0x3f800000L, 0x3f800000L, 0xbf800000L) shouldBe BigInt(0)
    }
  }
}
