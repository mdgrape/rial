package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.arith.RealSpec
import rial.math.MathFuncConfig
import rial.rng.GenRandomGaussianSim
import rial.table.{FuncTable, PolynomialGenerator}

/** Tests of the range check of the table generator.
 *
 * The evaluator of a math function reads a negative result as all ones, i.e.
 * as ~1.0, so generating a table whose polynomial could go below zero must
 * fail. Every default configuration generates its tables while it is built, so
 * a configuration that builds is one no row of which can go below zero.
 */
class PolynomialGeneratorTest extends AnyFunSpec {

  describe("the default configurations") {

    val configs = Seq(
      "simpleFP32"         -> (() => MathFuncConfig.simpleFP32),
      "simpleBF16"         -> (() => MathFuncConfig.simpleBF16),
      "standardFP32"       -> (() => MathFuncConfig.standardFP32),
      "standardBF16"       -> (() => MathFuncConfig.standardBF16),
      "standardScaledFP32" -> (() => MathFuncConfig.standardScaledFP32),
      "standardScaledBF16" -> (() => MathFuncConfig.standardScaledBF16),
      "standardFixedFP32"  -> (() => MathFuncConfig.standardFixedFP32),
      "standardFixedBF16"  -> (() => MathFuncConfig.standardFixedBF16),
      "allFP32"            -> (() => MathFuncConfig.allFP32),
      "allBF16"            -> (() => MathFuncConfig.allBF16),
    )

    configs.foreach { case (name, cfg) =>
      it(f"no row of a table of ${name} goes below zero") {
        info(f"${cfg().romLayout.rows} rows")
      }
    }

    it("no row of the probit table of GenRandomGaussian is out of range") {
      // it is read by an evaluator that clamps a negative result to zero, so
      // the rows that go slightly below zero near the shallow end are allowed
      GenRandomGaussianSim.genGaussian(BigInt(1) << 40, 64, RealSpec.Float32Spec)
    }
  }

  describe("a table that goes below zero") {

    // x - 0.3 is negative on most of [0, 1), and of order one, so the
    // polynomial of the interval is that function itself
    def generator(negativeIsClamped: Boolean): PolynomialGenerator = {
      val ft = new FuncTable(1)
      ft.addInterval((x: Double) => { x - 0.3 }, 0.0, 1.0)
      new PolynomialGenerator(ft, 18, 20, negativeIsClamped)
    }

    it("is refused, with the interval and its x range") {
      val gen = generator(false)
      val e = intercept[IllegalArgumentException] {
        gen.generateTable(gen.getBitWidth())
      }
      info(e.getMessage)
      assert(e.getMessage.contains("interval 0"))
      assert(e.getMessage.contains("x in [0.0, 1.0)"))
    }

    it("is accepted where the evaluator clamps a negative result to zero") {
      val gen = generator(true)
      gen.generateTable(gen.getBitWidth())
    }
  }
}
