// Regression for the FusedMulAddFPGeneric 1-ULP fused mis-round on EFFECTIVE
// SUBTRACTION (product and addend opposite signs). The shifted-out part of the
// smaller operand was handled by a sign-agnostic sticky / lost borrow, turning a
// round-to-even tie the wrong way. Pre-fix examples:
//   far-prod   0x88f0 * 0xefd8 + 0x81a1 -> 0x394b  (should be 0x394a)
//   far-addend 0xa1b7 * 0xb728 + 0x9be4 -> 0x9bdd  (should be 0x9bdc)
package rial.tests

import chisel3._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import rial.arith.{FusedMulAddFPGeneric, RealGeneric, RealSpec, RoundSpec}
import rial.util.PipelineStageConfig
import spire.math.SafeLong
import scala.util.Random

class FusedMulAddRoundingSpec extends AnyFlatSpec with ChiselScalatestTester with Matchers {
  private val s = RealSpec.BFloat16Spec

  private class Fma extends Module {
    val io = IO(new Bundle {
      val x = Input(UInt(16.W)); val y = Input(UInt(16.W))
      val z = Input(UInt(16.W)); val w = Output(UInt(16.W))
    })
    val u = Module(new FusedMulAddFPGeneric(s, s, s, s, RoundSpec.roundToEven, PipelineStageConfig.none))
    u.io.x := io.x; u.io.y := io.y; u.io.z := io.z; io.w := u.io.w
  }

  private def ideal(x: Int, y: Int, z: Int): Int =
    new RealGeneric(s, SafeLong(x))
      .fmadd(s, RoundSpec.roundToEven, new RealGeneric(s, SafeLong(y)), new RealGeneric(s, SafeLong(z)))
      .value.toBigInt.toInt & 0xffff

  private def hw(c: Fma, x: Int, y: Int, z: Int): Int = {
    c.io.x.poke(x.U(16.W)); c.io.y.poke(y.U(16.W)); c.io.z.poke(z.U(16.W))
    c.io.w.peek().litValue.toInt & 0xffff
  }

  private def randNormal(r: Random): Int = {
    val sgn = r.nextInt(2); val exp = 1 + r.nextInt(0xFE); val man = r.nextInt(0x80)
    (sgn << 15) | (exp << 7) | man
  }

  behavior of "FusedMulAddFPGeneric round-once correctness on effective subtraction"

  it should "round the documented far-prod / far-addend tie cases down (match RealGeneric.fmadd)" in {
    test(new Fma).withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
      hw(c, 0x88f0, 0xefd8, 0x81a1) shouldBe 0x394a // far-prod   (pre-fix 0x394b)
      hw(c, 0xa1b7, 0xb728, 0x9be4) shouldBe 0x9bdc // far-addend (pre-fix 0x9bdd)
    }
  }

  it should "be 0-ULP vs RealGeneric.fmadd over a normal-finite sweep (non-cancelling)" in {
    test(new Fma).withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
      val r = new Random(2024)
      var bad = 0
      val ex = scala.collection.mutable.ArrayBuffer[String]()
      var i = 0
      while (i < 8000) {
        val x = randNormal(r); val y = randNormal(r); val z = randNormal(r)
        val ref = ideal(x, y, z)
        val isZero = (ref & 0x7fff) == 0 // exact cancellation is a separate defect
        if (((ref >> 7) & 0xFF) != 0xFF && !isZero) {
          val g = hw(c, x, y, z)
          if (g != ref) { bad += 1; if (ex.size < 10) ex += f"x=0x$x%04x y=0x$y%04x z=0x$z%04x hw=0x$g%04x ideal=0x$ref%04x" }
        }
        i += 1
      }
      ex.foreach(info(_))
      withClue(s"$bad non-cancelling mismatches vs RealGeneric.fmadd: ") { bad shouldBe 0 }
    }
  }
}
