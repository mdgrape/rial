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
package rial.math

import chisel3._
import chisel3.util._

import rial.arith.FloatChiselUtil
import rial.arith.MulImpl
import rial.arith.RoundSpec
import rial.arith.SplitMultiplier
import rial.math.FuncKind._
import rial.util._

private[rial] class PostProcMultiplierOutput(val manW: Int) extends Bundle {
  val zman  = UInt(manW.W)
  val exInc = UInt(1.W)
}

private[rial] class PostProcMultiplier(
  val fracW:     Int,
  val manW:      Int,
  val roundSpec: RoundSpec,
  val stage:     PipelineStageConfig,
  val mulImpl:   MulImpl = MulImpl.plain,
) extends Module {

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = PostProcMultiplier.stagePoints(mulImpl)

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val lhs = UInt((1 + fracW).W) // polynomial result with hidden 1
      val rhs = UInt((1 + manW).W)  // runtime value with hidden 1
    }))
    val out = Valid(new PostProcMultiplierOutput(manW))
  })

  // --------------------------- stage "in" ---------------------------
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val lhsIn = cutIn(io.in.bits.lhs)
  val rhsIn = cutIn(io.in.bits.rhs)

  // ------------------------- stage "prod.*" -------------------------
  require(!mulImpl.isMacc,
    "the post-process multiplier has no addend to fold into a multiply-add; " +
    "ask for MulImpl.split to slice its product")

  val prod = if (mulImpl.isSplit) {
    val mult = Module(new SplitMultiplier(
      1 + fracW, 1 + manW, mulImpl.sliceW.get, stages.sub("prod")))

    mult.io.in.valid  := cutIn.valid
    mult.io.in.bits.x := lhsIn
    mult.io.in.bits.y := rhsIn
    mult.io.out.bits
  } else {
    // the product as it is written by hand: `MulImpl.plain` asks nothing of
    // the operands, so there is nothing to shape
    val cutProdIn  = new SharedPipe("prodIn",  cutIn.valid,     stages.at("prod.in"))
    val cutProdOut = new SharedPipe("prodOut", cutProdIn.valid, stages.at("prod.out"))
    cutProdOut(cutProdIn(lhsIn) * cutProdIn(rhsIn))
  }

  // latency of the product; the valid that leaves it is the one the rounding
  // below and the output register run on
  private val cutProd = new SharedPipe("prod", cutIn.valid, stages.latency("prod"))

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

  // --------------------------- stage "out" ---------------------------
  io.out := Pipe(cutProd.valid, outBits, stages.at("out"))
}

object PostProcMultiplier {

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.math.PostProcMultiplier]].
   *
   *  - `"prod.*"` the boundaries of the product: with `MulImpl.split` those
   *               of [[rial.arith.SplitMultiplier]], and `"prod.out"` alone
   *               for a plain product; the group name `"prod"` stands for
   *               `"prod.out"`
   *  - `"out"`    the rounded mantissa
   */
  def stagePoints(mulImpl: MulImpl): Seq[StageElem] = Seq(
    StageGroup("prod", if (mulImpl.isSplit) {
      SplitMultiplier.stagePoints()
    } else {
      Seq(StagePoint("out", 4.0))
    }),
    StagePoint("out", 2.0))

  /** True iff any configured function requires the post-process multiplier. */
  def isNeeded(funcs: Seq[FuncKind]): Boolean = {
    funcs.exists(FuncKind.needPostMult)
  }
}
