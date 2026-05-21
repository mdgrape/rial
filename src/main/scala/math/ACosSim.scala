//% @file ACosSim.scala
//
// Simulators for acos(x).
//
// acos is computed in two phases that pipeline-share the rest of the math
// block:
//
//   phase 1:  x -> y = sqrt(1 - |x|)           (uses the Sqrt sub-table)
//   phase 2:  y -> z = acos(1 - y^2) (= acos(|x|)), then negate when x < 0
//                                              (uses ACosPhase2's own table)
//
// Phase-1's output carries x's sign so phase-2 knows whether to do the
// `pi - z` step. Special values are passed through phase-1 with NaN-boxing
// (ex == all-ones, mantissa encoding the case) and decoded by phase-2.
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math._

import spire.math.SafeLong
import spire.implicits._

import rial.table._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding._

// Phase-2's table approximates  acos(1 - t^2)/t - 1 for t in [0, 1).
//  - t -> 0+: acos(1)/0 ... limit = sqrt(2), so f(0) = sqrt(2) - 1.
//  - t -> 1 : acos(0)/1 = pi/2, so f(1) = pi/2 - 1.
class ACosPhase2Table(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  private val func = (t: Double) => {
    if (t == 0.0) { sqrt(2.0) - 1.0 }
    else          { acos(1.0 - t * t) / t - 1.0 }
  }

  private val nIntervals = 1 << addrW
  private val ft  = FuncTable(func, nOrder, nIntervals, 0.0, 1.0)
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW)

  def getBitWidth(): PolynomialBitWidth = {gen.getBitWidth()}
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    gen.generateTable(bitW)
  }
}

// -------------------------------------------------------------------------
// phase 1
// -------------------------------------------------------------------------

object ACosPhase1SimGeneric {

  /** Compute y = sqrt(1 - |x|), preserving x's sign in the result. Special
   *  cases (x = 0, |x| >= 1, nan, inf) are encoded with NaN-boxing in the
   *  exponent so phase 2 can decode and emit the right constant.
   *
   *  @param tSqrt  the SqrtTable's coefficients (shared with the Sqrt block)
   *  @param x      input
   */
  def apply(tSqrt: PolynomialCoefficients, x: RealGeneric): RealGeneric = {

    val spec   = x.spec
    val manW   = spec.manW
    val exBias = spec.exBias

    val zsgn = x.sgn

    // -- special values (NaN-boxing) --
    //   mantissa MSB = 1 : ordinary NaN
    //   mantissa MSB = 0 : special:
    //                       1 -> x == 0   -> y = pi/2
    //                       2 -> x >= +1  -> y = 0
    //                       3 -> x <= -1  -> y = pi
    if (x.isNaN) {
      return new RealGeneric(spec, zsgn, maskI(spec.exW), SafeLong(1) << (spec.manW - 1))
    }
    if (x.isInfinite) {
      return new RealGeneric(spec, zsgn, maskI(spec.exW), SafeLong(1) << (spec.manW - 1))
    }
    if (x.ex == 0) {
      return new RealGeneric(spec, zsgn, maskI(spec.exW), SafeLong(1))
    }
    if (exBias <= x.ex) {
      // |x| >= 1; round to ±1.
      if (x.sgn == 0) {
        return new RealGeneric(spec, zsgn, maskI(spec.exW), SafeLong(2))
      } else {
        return new RealGeneric(spec, zsgn, maskI(spec.exW), SafeLong(3))
      }
    }

    assert(x.ex < exBias)
    assert(-1.0 < x.toDouble && x.toDouble < 1.0)

    // -- compute 1 - |x| (left in normalized form for sqrt) --
    val xShift      = exBias - x.ex
    val xShifted    = if (xShift > 1 + manW + 2) SafeLong(0) else (x.manW1 << 2) >> xShift
    val xSubtracted = (SafeLong(1) << (manW + 2)) - xShifted

    val xNormalizeShift = (1 + manW + 2) - binaryWidthSL(xSubtracted)
    val xNormalized     = xSubtracted << xNormalizeShift
    assert(bit(manW + 2, xNormalized) == 1 && xNormalized < (SafeLong(1) << (manW + 3)))

    val xRounded          = (xNormalized >> 2) + bit(1, xNormalized)
    val xRoundedMoreThan2 = bit(manW + 1, xRounded)
    val xConvertedEx      = exBias - xNormalizeShift + xRoundedMoreThan2
    val xConvertedMan     = slice(0, manW, xRounded)

    // 1 - |x| == 1 (i.e. |x| ~ 0)
    if (xConvertedEx == exBias) {
      assert(xConvertedMan == SafeLong(0))
      return new RealGeneric(spec, zsgn, maskI(spec.exW), SafeLong(1))
    }

    // -- run sqrt on (1 - |x|) reusing the Sqrt sub-table --
    val y    = new RealGeneric(spec, 0, xConvertedEx, xConvertedMan)
    val sqrt = SqrtSimGeneric(tSqrt, y)

    new RealGeneric(spec, zsgn, sqrt.ex, sqrt.man)
  }
}

// -------------------------------------------------------------------------
// phase 2
// -------------------------------------------------------------------------

object ACosPhase2SimGeneric {

  /** Compute acos(1 - x^2) from phase-1's output (with the original sign).
   *  Decodes the NaN-boxed special cases produced by phase-1.
   *
   *  @param tACos  ACosPhase2Table's coefficients
   *  @param x      phase-1's output (a RealGeneric with sgn = original x.sgn)
   */
  def apply(tACos: PolynomialCoefficients, x: RealGeneric): RealGeneric = {

    val spec   = x.spec
    val manW   = spec.manW
    val exBias = spec.exBias

    val addrW     = tACos.bitW.addrW
    val outPrec   = tACos.bitW.outPrec
    val extraBits = outPrec - manW
    val calcW     = manW + extraBits

    // -- decode phase-1's NaN-boxed special cases --
    if (x.ex == maskI(spec.exW)) {
      if      (x.man == SafeLong(1)) { return new RealGeneric(spec, Pi * 0.5) }
      else if (x.man == SafeLong(2)) { return new RealGeneric(spec, 0.0)      }
      else if (x.man == SafeLong(3)) { return new RealGeneric(spec, Pi)       }
      else                           { return RealGeneric.nan(spec)           }
    }

    assert(0.0 <= abs(x.toDouble) && abs(x.toDouble) <= 1.0,
      f"0.0 <= abs(x(${x.toDouble})) <= 1.0")

    // -- table lookup of acos(1 - t^2)/t - 1 --
    val xAligned = if (x.ex == 0) {
      SafeLong(0)
    } else {
      x.manW1 >> (exBias - x.ex)
    }

    val xW   = manW - addrW
    val addr = slice(xW, addrW, xAligned).toInt
    val xiv  = slice(0,  xW,    xAligned).toBigInt

    val polynomial = SafeLong(PolynomialEval.eval(addr, xiv, tACos))

    // -- multiply (acos(1-t^2)/t) * t to get acos(1-t^2) --
    val lhsManW1 = polynomial + (SafeLong(1) << calcW) // value in [1, 2), hidden 1 at calcW
    val rhsManW1 = x.manW1

    val zProd                = lhsManW1 * rhsManW1
    val zMoreThan2           = bit((calcW + 1) + (manW + 1) - 1, zProd)
    val zRoundBit            = calcW + zMoreThan2
    val zRounded             = roundBySpec(RoundSpec.roundToEven, zRoundBit, zProd)
    val zMoreThan2AfterRound = bit(manW + 1, zRounded)
    val zManW1               = zRounded >> zMoreThan2AfterRound

    val lhsEx = exBias            // (lhs has hidden 1 -> ex = 0 unbiased)
    val rhsEx = x.ex
    val zEx   = lhsEx + rhsEx - exBias + zMoreThan2 + zMoreThan2AfterRound

    if (x.sgn == 0) {
      return new RealGeneric(spec, 0, zEx, slice(0, manW, zManW1))
    }

    // -- for negative x, return pi - z --
    val pi        = new RealGeneric(spec, Pi)
    val zShift    = Seq(1 + manW + 2, exBias - zEx).min
    assert(zShift >= 0)
    val piAligned = pi.manW1 << 3
    val zAligned  = (zManW1 << 2) >> zShift
    val negZ      = piAligned - zAligned

    val negZNormalizeShift = (1 + manW + 3) - binaryWidthSL(negZ)
    val negZNormalized     = negZ << negZNormalizeShift
    assert(bit(manW + 3, negZNormalized) == 1)

    val negZRounded          = slice(3, manW, negZNormalized) + bit(2, negZNormalized)
    val negZRoundedMoreThan2 = bit(manW, negZRounded)
    val negZEx               = (exBias + 1 - negZNormalizeShift) + negZRoundedMoreThan2
    val negZMan              = slice(0, manW, negZRounded)

    new RealGeneric(spec, 0, negZEx, negZMan)
  }
}
