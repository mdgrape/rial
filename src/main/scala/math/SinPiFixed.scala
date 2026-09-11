//% @file SinPiFixed.scala
//
// sin(π x) and cos(π x) without either shared multiplier.
//
// The argument reduction is the one of sin/cos for an argument the caller has
// already divided by π, so the pre-process multiplier stays out of it (see
// SinCos.scala). What differs is what the table holds: sin(π y) itself, not
// sin(π y)/y, so the polynomial output is already the result as a fixed point
// value in [0, 1). Turning it into a float takes a leading-zero count and a
// shift instead of the post-process multiplier, which is what the module
// below does.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{FloatChiselUtil, RealSpec, RoundSpec}
import rial.util._

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

/** Turns [[SinCosPreProcess]]'s info and the polynomial output, a fixed point
 *  value in [0, 1), into the final sin(π x) / cos(π x) bit pattern. */
class SinPiFixedPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val fracW = bitW.outBitW

  /** bits of the normalized polynomial output below the mantissa. */
  val extraBits = fracW - 1 - manW

  require(extraBits >= 0,
    f"outBitW(${fracW}) must hold a hidden one above a manW(${manW}) mantissa")

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = SinPiFixedPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val pre  = Flipped(Valid(new SinCosPreProcInfo(spec)))
    val zres = Flipped(Valid(UInt(fracW.W)))
    val z    = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.pre.valid && io.zres.valid,
    stages.at("in"))

  val info = cutIn(io.pre.bits)
  val zres = cutIn(io.zres.bits)

  val countComb = PriorityEncoder(Reverse(zres))
  val zzeroComb = !zres.orR

  // --------------------------------------------------------------------------
  // stage "lead": the position of the leading one

  private val cutLead = new SharedPipe("lead", cutIn.valid, stages.at("lead"))

  val count     = cutLead(countComb)
  val zresLead  = cutLead(zres)
  val zzeroLead = cutLead(zzeroComb)
  val infoLead  = cutLead(info)

  // a value in [0, 1) with `count` leading zeros lies in
  // [2^-(count+1), 2^-count), so its exponent is -(count + 1).
  val zshiftedComb = (zresLead << count)(fracW - 1, 0)
  val zex0Comb     = (exBias - 1).U(exW.W) - count

  // --------------------------------------------------------------------------
  // stage "shift": the normalized polynomial output

  private val cutShift = new SharedPipe("shift", cutLead.valid, stages.at("shift"))

  val zshifted   = cutShift(zshiftedComb)
  val zex0       = cutShift(zex0Comb)
  val zzeroShift = cutShift(zzeroLead)
  val infoShift  = cutShift(infoLead)

  val zmanRoundedComb = WireDefault(0.U(manW.W))
  val zexRoundedComb  = WireDefault(0.U(exW.W))
  if (extraBits == 0) {
    zmanRoundedComb := zshifted(fracW - 2, 0)
    zexRoundedComb  := zex0
  } else {
    val lsb    = zshifted(extraBits)
    val round  = zshifted(extraBits - 1)
    val sticky = if (extraBits == 1) {
      false.B
    } else {
      zshifted(extraBits - 2, 0).orR
    }
    val roundInc = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
      lsb, round, sticky)
    // the carry of the rounding leaves a zero mantissa one exponent higher
    val zman0 = zshifted(fracW - 2, extraBits) +& roundInc
    zmanRoundedComb := zman0(manW - 1, 0)
    zexRoundedComb  := zex0 + zman0(manW)
  }

  // --------------------------------------------------------------------------
  // stage "round": the rounded mantissa and the exponent it leaves

  private val cutRound = new SharedPipe("round", cutShift.valid, stages.at("round"))

  val zmanRounded = cutRound(zmanRoundedComb)
  val zexRounded  = cutRound(zexRoundedComb)
  val zzero       = cutRound(zzeroShift)
  val infoRound   = cutRound(infoShift)

  // --------------------------------------------------------------------------
  // stage "out"

  val znan = infoRound.znan
  val zone = infoRound.zone
  val zsgn = infoRound.zsgn

  val zex = Mux(znan,  Fill(exW, 1.U(1.W)),
            Mux(zone,  exBias.U(exW.W),
            Mux(zzero, 0.U(exW.W),
                       zexRounded)))
  val zman = Mux(znan,          Cat(1.U(1.W), 0.U((manW - 1).W)),
             Mux(zone || zzero, 0.U(manW.W),
                                zmanRounded))

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := Pipe(cutRound.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.SinPiFixedPostProcess]]. */
object SinPiFixedPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"lead"`  the position of the leading one of the polynomial output
   *  - `"shift"` the normalized polynomial output and its exponent
   *  - `"round"` the rounded mantissa
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("lead",  1.0),
    StagePoint("shift", 1.0),
    StagePoint("round", 1.0),
    StagePoint("out",   1.0))
}
