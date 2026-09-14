package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._
import rial.arith.RealSpec

/** The exp polynomial must never round up into the hidden bit.
 *
 *  ExpPostProcess rounds the polynomial result to `manW` bits and builds the
 *  correction term from `2^manW + zman`; a rounding that carries would take
 *  the hidden bit from the wrong place.  ExpTableCoeff forbids it at
 *  elaboration, and these checks run the same computation without a
 *  simulation.
 */
class ExpTableCarryTest extends AnyFunSpec {

  def threshold(cfg: MathFuncConfig, spec: RealSpec): BigInt = {
    val padding = cfg.bitW.outBitW - spec.manW
    require(padding > 0)
    (BigInt(1) << cfg.bitW.outBitW) - (BigInt(1) << (padding - 1))
  }

  // Pow2 reads the same table as Exp, at the argument the caller has scaled.
  val configs = Seq(
    ("standardFP32", MathFuncConfig.standardFP32, RealSpec.Float32Spec,  Exp),
    ("standardBF16", MathFuncConfig.standardBF16, RealSpec.BFloat16Spec, Exp),
    ("allFP32",      MathFuncConfig.allFP32,      RealSpec.Float32Spec,  Exp),
    ("allBF16",      MathFuncConfig.allBF16,      RealSpec.BFloat16Spec, Exp),
    ("standardFixedFP32", MathFuncConfig.standardFixedFP32,
                          RealSpec.Float32Spec,  Pow2),
    ("standardFixedBF16", MathFuncConfig.standardFixedBF16,
                          RealSpec.BFloat16Spec, Pow2),
  )

  describe("ExpTableCoeff.maxPolynomial stays below the carry threshold") {

    for ((name, cfg, spec, fn) <- configs) {
      it(f"$name") {
        val maxPoly = ExpTableCoeff.maxPolynomial(cfg.coeffs(fn), spec)
        val limit   = threshold(cfg, spec)
        val outBitW = cfg.bitW.outBitW

        info(f"$name: max = 2^$outBitW - ${(BigInt(1) << outBitW) - maxPoly}" +
             f", threshold = 2^$outBitW - ${(BigInt(1) << outBitW) - limit}" +
             f", margin = ${limit - maxPoly}")

        assert(maxPoly < limit)
      }
    }
  }

  describe("ExpTableCoeff.maxPolynomial reproduces the measured maxima") {

    it("standardFP32 reaches 2^26 - 12") {
      assert(ExpTableCoeff.maxPolynomial(MathFuncConfig.standardFP32.coeffs(Exp),
             RealSpec.Float32Spec) == (BigInt(1) << 26) - 12)
    }

    it("standardBF16 reaches 2^8 - 3") {
      assert(ExpTableCoeff.maxPolynomial(MathFuncConfig.standardBF16.coeffs(Exp),
             RealSpec.BFloat16Spec) == (BigInt(1) << 8) - 3)
    }
  }
}
