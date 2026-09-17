//% @file PreProcMultiplier.scala
//
// Sizing of the shared pre-process multiplier of the math functions that need
// an `x * constant` step before the polynomial stage (sin/cos, exp, ...).
//
// The product itself is written by the orchestrating module, which computes
// the minimum lhsW required by the configured functions and only builds it
// when needed.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import rial.arith.RealSpec
import rial.math.FuncKind._

/** Which functions need the shared pre-process multiplier, and how wide its
 *  left hand operand has to be for them.
 *
 * @see [[rial.arith.MulImpl]]
 */
object PreProcMultiplier {

  /** True iff any configured function requires the pre-process multiplier. */
  def isNeeded(funcs: Seq[FuncKind]): Boolean = {
    funcs.exists(FuncKind.needPreMult)
  }

  /** Minimum `lhsW` that covers every pre-mult-needing function in `funcs`.
   *  Returns 0 if none of them need the multiplier (in which case the caller
   *  should not instantiate it at all).
   */
  def lhsW(funcs: Seq[FuncKind], spec: RealSpec): Int = {
    val ws = funcs.collect {
      case Sin | Cos => SinCosPreMulArgs.lhsW(spec)
      case Exp       => ExpPreMulArgs.lhsW(spec)
      case ScaleMixtureGaussian => ScaleMixtureGaussianPreMulArgs.lhsW(spec)
    }
    if (ws.isEmpty) { 0 } else { ws.max }
  }
}
