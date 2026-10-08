package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._

/** The facade the JNI binding builds its tables through.
 *
 * Every preset of [[rial.math.MathSim]] names a [[rial.math.MathFuncConfig]]
 * preset, and the C++ side mirrors its numbers, so the two have to end up at
 * the same table and polynomial widths.
 */
class MathSimTest extends AnyFunSpec {

  describe("MathSim.standardFP32") {

    val sim = MathSim.standardFP32
    val cfg = MathFuncConfig.standardFP32

    it("supports the same functions as MathFuncConfig.standardFP32") {
      assert(sim.config.funcs.toSet == cfg.funcs.toSet)
    }

    it("reads every table at the width that config reads it at") {
      for (fn <- cfg.funcs) {
        assert(sim.config.bitWOf(fn) == cfg.bitWOf(fn),
          f"${FuncKind.getString(fn)}: ${sim.config.bitWOf(fn)} != " +
          f"${cfg.bitWOf(fn)}")
      }
    }
  }

  describe("MathSim: a table address width per function") {

    it("reads every table at addrW when no width is given") {
      // the position is as wide as a flat address of 8 bits leaves, so that
      // the table reaches every position of it
      val sim = new MathSim(8, 0x7F, 23, 2, 8, Array.empty[Int], 15, 26,
        Array("sqrt", "reciprocal"), Double.NaN, Double.NaN)
      assert(sim.config.bitWOf(Sqrt).addrW == 8)
      assert(sim.config.bitWOf(Reciprocal).addrW == 8)
    }

    it("gives each name its own width") {
      val sim = new MathSim(8, 0x7F, 23, 2, 8, Array(5, 7), 18, 26,
        Array("sqrt", "reciprocal"), Double.NaN, Double.NaN)
      assert(sim.config.bitWOf(Sqrt).addrW == 5)
      assert(sim.config.bitWOf(Reciprocal).addrW == 7)
    }

    it("hands a width to the function that owns the table") {
      val sim = new MathSim(8, 0x7F, 23, 2, 8, Array(6, 6), 18, 26,
        Array("sin", "cos"), Double.NaN, Double.NaN)
      assert(sim.config.bitWOf(Sin).addrW == 6)
      assert(sim.config.bitWOf(Cos).addrW == 6)
    }

    it("rejects two widths for one table") {
      assertThrows[IllegalArgumentException] {
        new MathSim(8, 0x7F, 23, 2, 8, Array(6, 7), 18, 26,
          Array("sin", "cos"), Double.NaN, Double.NaN)
      }
    }

    it("gives the sqrt table its width when acos names it") {
      val sim = new MathSim(8, 0x7F, 23, 2, 8, Array(5, 7), 18, 26,
        Array("sqrt", "acos"), Double.NaN, Double.NaN)
      assert(sim.config.bitWOf(Sqrt).addrW       == 5)
      assert(sim.config.bitWOf(ACosPhase1).addrW == 5)
      assert(sim.config.bitWOf(ACosPhase2).addrW == 7)
    }

    it("rejects a width per function that is not one per function") {
      assertThrows[IllegalArgumentException] {
        new MathSim(8, 0x7F, 23, 2, 8, Array(5), 18, 26,
          Array("sqrt", "reciprocal"), Double.NaN, Double.NaN)
      }
    }
  }
}
