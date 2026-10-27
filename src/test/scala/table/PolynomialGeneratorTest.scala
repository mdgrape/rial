package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.arith.RealSpec
import rial.math.MathFuncConfig
import rial.rng.GenRandomGaussianSim
import rial.table.{FindMinMax, FuncTable, PolynomialGenerator}

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

  describe("a polynomial with two turning points inside its interval") {

    // 0.5 + 0.3 (t^3 - t) + 0.05 t^2 with t = 2x - 1 in [-1, 1], which stays
    // inside [0.35, 0.65] and turns at t = -0.64 and t = 0.52; every
    // coefficient is nonzero, so each has a width
    lazy val gen = {
      val ft = new FuncTable(3)
      ft.addInterval((x: Double) => {
        val t = 2.0 * x - 1.0
        0.5 + 0.3 * (t * t * t - t) + 0.05 * t * t
      }, 0.0, 1.0)
      new PolynomialGenerator(ft, 18, 20)
    }
    lazy val bitW  = gen.getBitWidth()
    lazy val table = gen.generateTable(bitW)

    def fits(v: BigInt, w: Int): Boolean = {
      -(BigInt(1) << (w - 1)) <= v && v < (BigInt(1) << (w - 1))
    }

    it("keeps every Horner intermediate inside calcBits at every position") {
      val cs = table.coeffs(0)
      assert(fits(BigInt(cs(bitW.nOrder)), bitW.calcBits(bitW.nOrder)))

      for (xu <- 0 until (1 << bitW.inBitW)) {
        // the steps of PolynomialEval.eval
        val x   = BigInt(xu) - (BigInt(1) << bitW.inPrec)
        var sum = BigInt(cs(bitW.nOrder))
        for (k <- (bitW.nOrder - 1) to 0 by -1) {
          val xminW = math.min(bitW.calcBits(k + 1), bitW.inPrec)
          sum = (((x >> (bitW.inPrec - xminW)) * sum) >> xminW) + cs(k)
          assert(fits(sum, bitW.calcBits(k)),
            f"x=${xu}, step ${k}: ${sum} does not fit ${bitW.calcBits(k)} bits")
        }
        assert(0 <= sum && sum < (BigInt(1) << bitW.outBitW),
          f"x=${xu}: ${sum} is out of [0, 2^${bitW.outBitW})")
      }
    }

    it("has its extrema at the turning points") {
      val ext = 2.0 / (3.0 * math.sqrt(3.0))
      val (lo, hi) = FindMinMax.polynomial(Seq(0.0, -1.0, 0.0, 1.0), -1.0, 1.0)
      assert(math.abs(lo + ext) < 1e-12)
      assert(math.abs(hi - ext) < 1e-12)

      // (t^2 - 1/4)^2 = t^4 - t^2/2 + 1/16, with a zero leading coefficient
      val (lo4, hi4) = FindMinMax.polynomial(
        Seq(1.0 / 16, 0.0, -0.5, 0.0, 1.0, 0.0), -1.0, 1.0)
      assert(math.abs(lo4) < 1e-12)
      assert(math.abs(hi4 - 9.0 / 16) < 1e-12)
    }
  }
}
