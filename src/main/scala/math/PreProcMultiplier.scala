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
package rial.math

import chisel3._
import chisel3.util._

import rial.arith.RealSpec
import rial.math.FuncKind._
import rial.util._

private[rial] class PreProcMultiplier(
  val lhsW:  Int,
  val manW:  Int,
  val stage: PipelineStageConfig,
) extends Module {

  val outW = lhsW + (1 + manW)

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = PreProcMultiplier.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val lhs = UInt(lhsW.W)
      val rhs = UInt((1 + manW).W)
    }))
    val out = Valid(UInt(outW.W))
  })

  // --------------------------- stage "in" ---------------------------
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val lhsIn = cutIn(io.in.bits.lhs)
  val rhsIn = cutIn(io.in.bits.rhs)

  // upstream zeroes the operands when the multiplier isn't selected
  val prodComb = lhsIn * rhsIn

  // -------------------------- stage "prod" --------------------------
  private val cutProd = new SharedPipe("prod", cutIn.valid, stages.at("prod"))

  val prod = cutProd(prodComb)

  // --------------------------- stage "out" --------------------------
  io.out := Pipe(cutProd.valid, prod, stages.at("out"))
}

object PreProcMultiplier {

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.math.PreProcMultiplier]].
   *
   *  - `"prod"` the partial products
   *  - `"out"`  the summed product
   *
   * A wide product is several multiplier blocks whose outputs are summed, so
   * the two boundaries let the sum have a pipeline stage of its own.
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("prod", 3.0),
    StagePoint("out",  1.0))

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
