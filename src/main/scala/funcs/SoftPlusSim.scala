//% @file SoftPlusSim.scala
//
// Simulator for softplus(x) = log(1 + exp(x)).
//
// Two sub-tables, selected by x.sgn:
//   pos table (x >= 0): approximates log(1+exp(x01*rangeMax)) / (2*rangeMax),
//                       so the polynomial result has hidden exponent
//                       (rangeMaxLog2 + 1) baked into the final result.
//   neg table (x <  0): approximates log(1 + exp(-x01*rangeMax))  in [0, ln(2)].
//
// `rangeMaxLog2` is chosen so that beyond rangeMax the function value rounds
// to either `x` itself (x>0) or to 0 (x<0) within FP precision.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import scala.math._

import spire.math.SafeLong
import spire.implicits._

import rial.approx._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric

class SoftPlusTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  val rangeMaxLog2 = SoftPlusTable.rangeMaxLog2(outBitW)
  val rangeMax     = 1 << rangeMaxLog2

  // pos branch: log(1 + exp(x01*rangeMax)) / (2*rangeMax)  in [log(2)/(2*rangeMax), 0.5)
  private val funcPos = (x01: Double) => {
    val x = x01 * rangeMax
    val z = log(1.0 + exp(x)) / (rangeMax * 2.0)
    assert(z <= 1.0, f"funcPos: x = ${x}, z = ${z}")
    z
  }
  // neg branch: log(1 + exp(-x01*rangeMax))  in (0, log(2)]
  private val funcNeg = (x01: Double) => {
    val x = -x01 * rangeMax
    val z = log(1.0 + exp(x))
    assert(z <= 1.0, f"funcNeg: x = ${x}, z = ${z}")
    z
  }

  private val nIntervals = 1 << addrW
  private val ftPos = FuncTable(funcPos, nOrder, nIntervals, 0.0, 1.0)
  private val ftNeg = FuncTable(funcNeg, nOrder, nIntervals, 0.0, 1.0)
  private val genPos = new PolynomialGenerator(ftPos, inBitW, outBitW)
  private val genNeg = new PolynomialGenerator(ftNeg, inBitW, outBitW)

  def getBitWidth(): PolynomialBitWidth = {
    genPos.getBitWidth().merge(genNeg.getBitWidth())
  }

  def genPosCoeffs(bitW: PolynomialBitWidth): PolynomialCoefficients = { genPos.generateTable(bitW) }
  def genNegCoeffs(bitW: PolynomialBitWidth): PolynomialCoefficients = { genNeg.generateTable(bitW) }

  /** Concatenated form: [pos] ++ [neg].  Consumed by the HW (SoftPlusTableCoeff)
   *  which slices it back into two independent ROMs and selects with x.sgn. */
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    val cPos = genPosCoeffs(bitW)
    val cNeg = genNegCoeffs(bitW)
    PolynomialCoefficients(bitW, cPos.coeffs ++ cNeg.coeffs)
  }
}

object SoftPlusTable {

  /** Smallest power-of-two threshold above which softplus(|x|) is
   *  indistinguishable from `x` (x>0) or 0 (x<0) within `precW` bits of
   *  precision.  Mirrors the original math/ port.
   */
  def rangeMaxLog2(precW: Int): Int = {
    val log2  = (x: Double) => log(x) / log(2.0)
    val delta = pow(2.0, -precW)
    // |softplus(x) - x|       < delta  (for x > 0)
    //  softplus(-|x|)         < delta  (for x < 0)
    //   both reduce to:   x > log(1 / (exp(delta) - 1))
    val threshold = log(1.0 / (exp(delta) - 1.0))
    ceil(log2(threshold)).toInt
  }

  def fp32(): SoftPlusTable = {
    val addrW     = 8
    val extraBits = 3
    new SoftPlusTable(
      nOrder  = 2,
      addrW   = addrW,
      inBitW  = (RealSpec.Float32Spec.manW - addrW),
      outBitW = (RealSpec.Float32Spec.manW + extraBits),
    )
  }
  def bf16(): SoftPlusTable = {
    val addrW     = RealSpec.BFloat16Spec.manW
    val extraBits = 0
    new SoftPlusTable(
      nOrder  = 0,
      addrW   = addrW,
      inBitW  = (RealSpec.BFloat16Spec.manW - addrW),
      outBitW = (RealSpec.BFloat16Spec.manW + extraBits),
    )
  }

  private lazy val fp32T   = fp32()
  private lazy val fp32Pos = fp32T.genPosCoeffs(fp32T.getBitWidth())
  private lazy val fp32Neg = fp32T.genNegCoeffs(fp32T.getBitWidth())
  def softplusDefaultFP32(x: Float): Float = {
    SoftPlusSimGeneric(fp32Pos, fp32Neg,
      RealGeneric.fromFloat(RealSpec.Float32Spec, x)).toFloat
  }
}

/** An interface for software simulator.
 *
 *  Takes the two sub-tables separately (mirrors LogSimGeneric's signature).
 */
object SoftPlusSimGeneric {

  def apply(
    tPos: PolynomialCoefficients,
    tNeg: PolynomialCoefficients,
    x:    RealGeneric,
  ): RealGeneric = {

    if (x.isNaN)                    { return RealGeneric.nan (x.spec)    }
    if (x.isInfinite && x.sgn == 0) { return RealGeneric.inf (x.spec, 0) }
    if (x.isInfinite && x.sgn == 1) { return RealGeneric.zero(x.spec)    }

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val addrW     = tPos.bitW.addrW
    val outPrec   = tPos.bitW.outPrec
    val calcW     = outPrec
    val extraBits = outPrec - manW
    assert(tNeg.bitW.addrW   == addrW)
    assert(tNeg.bitW.outPrec == outPrec)

    val rangeMaxLog2 = SoftPlusTable.rangeMaxLog2(outPrec)
    val rangeMaxEx   = exBias + rangeMaxLog2 - 1

    // ---- saturation: large |x| ----
    if (x.ex > rangeMaxEx) {
      if (x.sgn == 0) { return x                       }
      else            { return RealGeneric.zero(spec)  }
    }

    // ---- scale |x| into [0, rangeMax) fixed-point ----
    val xScaled = x.manW1 >> (rangeMaxEx - x.ex + 1)

    // ---- polynomial table query ----
    val t = if (x.sgn == 0) tPos else tNeg
    val xW   = manW - addrW
    val addr = slice(xW, addrW, xScaled).toInt
    val xiv  = slice(0,  xW,    xScaled).toBigInt
    val zScaled0 = PolynomialEval.eval(addr, xiv, t)

    val zScaled = if (zScaled0 < 0) {
      println(f"WARNING (${this.getClass.getName}) : polynomial value negative at x=${x.toDouble}")
      SafeLong(0)
    } else if (zScaled0 >= (SafeLong(1) << calcW)) {
      println(f"WARNING (${this.getClass.getName}) : polynomial range overflow at x=${x.toDouble}")
      maskSL(calcW)
    } else {
      SafeLong(zScaled0)
    }

    // ---- normalize and round ----
    val zShift = (calcW + 1) - binaryWidthSL(zScaled)
    val zmanW1 = zScaled << zShift
    val zEx    = if (x.sgn == 0) {
      exBias - zShift + (rangeMaxLog2 + 1) // pos table embeds factor of 2*rangeMax
    } else {
      exBias - zShift
    }
    val zmanW1Rounded = if (extraBits > 0) {
      (zmanW1 >> extraBits) + bit(extraBits - 1, zmanW1)
    } else {
      zmanW1
    }

    val zMan = if (zmanW1Rounded >= (SafeLong(1) << (manW + 1))) {
      println(f"WARNING (${this.getClass.getName}) : rounding overflow at x=${x.toDouble}, zmanW1 = ${zmanW1Rounded.toLong.toBinaryString}")
      maskSL(manW)
    } else {
      slice(0, manW, zmanW1Rounded)
    }

    new RealGeneric(spec, 0, zEx, zMan)
  }
}
