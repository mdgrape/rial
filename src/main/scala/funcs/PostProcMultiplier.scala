//% @file PostProcMultiplier.scala
//
// Shared post-process multiplier for math functions whose polynomial table
// produces some quantity that must be multiplied by a runtime value before
// the result becomes the final mantissa (sin/cos: poly*y, log: poly*(1-x),
// acos2: poly*x, ...).
//
// Performs the (1+fracW)*(1+manW) multiply and the round-to-manW directly,
// returning the rounded mantissa together with a 1-bit exponent-increment
// signal that the per-function post-process consumes.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import chisel3._
import chisel3.util._

import rial.arith.RoundSpec
import rial.arith.FloatChiselUtil
import rial.funcs.FuncKind._

private[rial] class PostProcMultiplierOutput(val manW: Int) extends Bundle {
  val zman  = UInt(manW.W)
  val exInc = UInt(1.W)
}

private[rial] class PostProcMultiplier(
  val fracW:     Int,
  val manW:      Int,
  val roundSpec: RoundSpec,
  val stage:     Int,
) extends Module {

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val lhs = UInt((1 + fracW).W) // polynomial result with hidden 1
      val rhs = UInt((1 + manW).W)  // runtime value with hidden 1
    }))
    val out = Valid(new PostProcMultiplierOutput(manW))
  })

  val prod      = io.in.bits.lhs * io.in.bits.rhs
  val moreThan2 = prod((manW + 1) + (fracW + 1) - 1)

  val roundBits = fracW
  val sticky    = prod(roundBits - 2, 0).orR | (moreThan2 & prod(roundBits - 1))
  val round     = Mux(moreThan2, prod(roundBits),     prod(roundBits - 1))
  val shifted   = Mux(moreThan2, prod((manW + 1) + (fracW + 1) - 2, (fracW + 1)    ),
                                 prod((manW + 1) + (fracW + 1) - 3, (fracW + 1) - 1))
  val lsb       = shifted(0)
  val roundInc  = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky)
  val rounded   = shifted +& roundInc
  val moreThan2AfterRound = rounded(manW)

  val outBits = Wire(new PostProcMultiplierOutput(manW))
  outBits.zman  := rounded(manW - 1, 0)
  outBits.exInc := moreThan2 | moreThan2AfterRound

  io.out := Pipe(io.in.valid, outBits, stage)
}

object PostProcMultiplier {

  /** True iff any configured function requires the post-process multiplier. */
  def isNeeded(funcs: Seq[FuncKind]): Boolean = {
    funcs.exists(FuncKind.needPostMult)
  }
}
