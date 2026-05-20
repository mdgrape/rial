//% @file PreProcMultiplier.scala
//
// Shared pre-process multiplier for math functions that need an
// `x * constant` step before the polynomial stage (sin/cos, exp, ...).
//
// Sized externally: the orchestrating module computes the minimum lhsW
// required by the configured functions (via `PreProcMultiplier.lhsW`)
// and only instantiates this module when needed.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import chisel3._
import chisel3.util._

import rial.arith.RealSpec
import rial.funcs.FuncKind._

private[rial] class PreProcMultiplier(
  val lhsW:  Int,
  val manW:  Int,
  val stage: Int,
) extends Module {

  val outW = lhsW + (1 + manW)

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val lhs = UInt(lhsW.W)
      val rhs = UInt((1 + manW).W)
    }))
    val out = Valid(UInt(outW.W))
  })

  // upstream zeroes the operands when the multiplier isn't selected
  val prod = io.in.bits.lhs * io.in.bits.rhs
  io.out := Pipe(io.in.valid, prod, stage)
}

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
