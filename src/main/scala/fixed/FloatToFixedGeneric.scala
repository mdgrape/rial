//% @file FloatToFixedGeneric.scala
//
// Fixed point -> floating point
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.fixed

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import rial.arith._
import rial.fixed._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._

/** Floating point number -> Fixed point number
 *
 * {{{
 * class FloatToFixedGeneric(...) extends Module {
 *   val io = IO(new Bundle{
 *     val x = Input (UInt(xSpec.W.W))
 *     val y = Input (UInt(w.W)) // defined only if `yintW` is defined
 *     val z = Output(UInt(zSpec.W.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * First input is the input floating point number.
 * Second input (optional) is the current fraction width.
 * If yintW is not given (be `None`), it uses zSpec.fracW that will be hard-coded.
 *
 * 0 means the fixed point value is an integer (w/o fractional part).
 * Since, in most cases, fixed point is signed, it results a weird result if
 * fracW > totalW in a fixed point value.
 * Take care of the width of the second input, y.
 *
 * The latency of this module can be configured via [[rial.util.PipelineStageConfig]].
 * The converter declares the boundaries of its datapath, so the registers can
 * be placed inside it - either automatically or by name:
 *
 *  - `"in"`     the module input
 *  - `"unpack"` after decompose / classify
 *  - `"align"`  after the alignment shift
 *  - `"round"`  after the rounding increment and the negation, and after the
 *               exponent range check that runs in parallel with them
 *  - `"out"`    the module output
 *
 * {{{
 * // three registers, placed by the module
 * new FloatToFixedGeneric(xSpec, None, zSpec, RoundSpec.roundToEven, 3)
 * // hand-picked positions
 * new FloatToFixedGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.specify("unpack" -> 1, "align" -> 1, "round" -> 1))
 * // everything at the output, as before
 * new FloatToFixedGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.atOut(3))
 * }}}
 *
 * @param xSpec spec of the input floating point number.
 * @param yintW width of the second input, fracWidth value. If None, zSpec.fracW is used instead.
 * @param zSpec spec of the output fixed point number.
 * @param roundSpec spec of the rounding
 * @param stage Determines how many registers needed.
 *
 * For parameters, see also:
 * @see [[rial.arith.RealSpec]]
 * @see [[rial.fixed.FixedSpec]]
 * @see [[rial.arith.RoundSpec]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class FloatToFixedGeneric(
  xSpec : RealSpec, yintW: Option[Int], zSpec : FixedSpec,
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  /** Boundaries of the datapath at which pipeline registers can be placed.
   * The weights are rough, target-independent estimates of the delay of the
   * logic that ends at each boundary. */
  val stagePoints = FloatToFixedGeneric.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getParam = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x = Input (UInt(xSpec.W.W))
    val y = yintW.map(w => Input(UInt(w.W)))
    val z = Output(UInt(zSpec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  val xIn = ShiftRegister(io.x, stages.at("in"))
  val yIn = io.y.map(y => ShiftRegister(y, stages.at("in")))

  //----------------------------------------------------------------------
  // S1 : decompose / classify
  val (xsgn,  xex,  xmanW1) = FloatChiselUtil.decomposeWithLeading1(xSpec, xIn)
  val (xzero, xinf, xnan)   = FloatChiselUtil.checkValue(xSpec, xIn)

  //--------------------------- stage "unpack" ---------------------------
  private def cutUnpack[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("unpack"))
  }

  val xmanW1Unpack = cutUnpack(xmanW1)
  val xexUnpack    = cutUnpack(xex)
  val xsgnUnpack   = cutUnpack(xsgn)
  val xzeroUnpack  = cutUnpack(xzero)
  val xinfUnpack   = cutUnpack(xinf)
  val xnanUnpack   = cutUnpack(xnan)
  val yUnpack      = yIn.map(y => cutUnpack(y))

  //----------------------------------------------------------------------
  // S2 : alignment shift
  val zTmp      = Cat(xmanW1Unpack, Fill(zSpec.W, 0.U(1.W)))

  // one extra bit so a full-width shift (= zSpec.W, needed when the
  // scaled magnitude lies in [0.5, 1) and only the round/sticky bits
  // survive) does not wrap to 0. Larger shifts are covered by the
  // isZZero guard below.
  val zTmpShiftBits = log2Up(zSpec.W) + 1
  val zTmpShift     = if(yintW.isEmpty) {
    ((zSpec.intW-1 + xSpec.exBias).U - xexUnpack)(zTmpShiftBits-1, 0)
  } else {
    ((zSpec.W-1 + xSpec.exBias).U - yUnpack.get - xexUnpack)(zTmpShiftBits-1, 0)
  }
  val zTmpShifted = zTmp >> zTmpShift

  //--------------------------- stage "align" ----------------------------
  private def cutAlign[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("align"))
  }

  val zTmpShiftedAlign = cutAlign(zTmpShifted)
  val xsgnAlign        = cutAlign(xsgnUnpack)

  //----------------------------------------------------------------------
  // S3 : rounding increment and negation; in parallel, the exponent range
  //      check (see the special case section below)
  val zTmpLSB    = zTmpShiftedAlign(xSpec.manW+1).asBool
  val zTmpRound  = zTmpShiftedAlign(xSpec.manW  ).asBool
  val zTmpSticky = zTmpShiftedAlign(xSpec.manW-1, 0).orR
  val zTmpInc    = FloatChiselUtil.roundIncBySpec(roundSpec, zTmpLSB, zTmpRound, zTmpSticky)

  val zAbs = zTmpShiftedAlign(zSpec.W + xSpec.manW, xSpec.manW+1) + zTmpInc.asUInt
  val z    = Mux(xsgnAlign === 0.U, zAbs, ~zAbs + 1.U)

  // -----------------------------------------------------------------------
  // special cases; this path runs in parallel with the alignment shift and
  // the rounding, so it crosses both the "align" and the "round" boundary

  val zPosInf = if(zSpec.signed) {
    ((1L<<(zSpec.W-1))-1L).U(zSpec.W.W)
  } else {
    ((1L<<zSpec.W)-1L).U(zSpec.W.W)
  }
  val zNegInf = if(zSpec.signed) {
    (1L<<(zSpec.W-1L)).U(zSpec.W.W) // ?
  } else {
    0.U(zSpec.W.W)
  }

  val zExMaxBiased = if(yintW.isEmpty) {
    (xSpec.exBias + zSpec.intW - 1).U
  } else {
    (xSpec.exBias + zSpec.W - 1).U - yUnpack.get // intW = totalW - fracW
  }
  val zExMinBiased = if(yUnpack.isEmpty) {
    val zExMin = -(zSpec.fracW + 1)
    max(0, xSpec.exBias + zExMin).U
  } else {
    xSpec.exBias.U - yUnpack.get - 1.U
  }

  val isZInf    = xinfUnpack  || zExMaxBiased < xexUnpack
  val isZZero   = xzeroUnpack || xexUnpack < zExMinBiased
  val isZPosInf = isZInf && (xsgnUnpack === 0.U)
  val isZNegInf = isZInf && (xsgnUnpack === 1.U)

  val xnanAlign      = cutAlign(xnanUnpack)
  val isZPosInfAlign = cutAlign(isZPosInf)
  val isZNegInfAlign = cutAlign(isZNegInf)
  val isZZeroAlign   = cutAlign(isZZero)

  //--------------------------- stage "round" ----------------------------
  private def cutRound[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("round"))
  }

  val zRound         = cutRound(z)
  val xnanRound      = cutRound(xnanAlign)
  val isZPosInfRound = cutRound(isZPosInfAlign)
  val isZNegInfRound = cutRound(isZNegInfAlign)
  val isZZeroRound   = cutRound(isZZeroAlign)

  //----------------------------------------------------------------------
  // S4 : Final (special case muxes)
  val z0 = MuxCase(zRound, Seq(
    xnanRound      -> zPosInf,
    isZPosInfRound -> zPosInf,
    isZNegInfRound -> zNegInf,
    isZZeroRound   -> 0.U(zSpec.W.W)
    ))

  io.z := ShiftRegister(z0, stages.at("out"))
}

/** Pipeline stage points of [[rial.fixed.FloatToFixedGeneric]].
 *
 * @see [[rial.util.StagePoint]]
 */
object FloatToFixedGeneric {
  /** Ordered boundaries of the datapath, with a rough, target-independent
   *  estimate of the delay of the logic ending at each. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("unpack", 2), StagePoint("align", 3),
    StagePoint("round",  3), StagePoint("out",   1))
}
