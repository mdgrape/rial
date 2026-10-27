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
 *     val in  = Flipped(Valid(new FloatToFixedGeneric.Input(xSpec, yintW)))
 *     val out = Valid(UInt(zSpec.W.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * First input, `in.bits.x`, is the input floating point number.
 * Second input (optional), `in.bits.y`, is the current fraction width. It is
 * defined only if `yintW` is defined.
 * If yintW is not given (be `None`), it uses zSpec.fracW that will be hard-coded.
 *
 * 0 means the fixed point value is an integer (w/o fractional part).
 * Since, in most cases, fixed point is signed, it results a weird result if
 * fracW > totalW in a fixed point value.
 * Take care of the width of the second input, y.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"unpack"` after the operand unpacking
 *  - `"align"`  after the alignment shift
 *  - `"round"`  after the rounding
 *  - `"out"`    the module output
 *
 * [[getStage]] gives the resulting latency, and `out.valid` is `in.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // three registers, placed by the module
 * new FloatToFixedGeneric(xSpec, None, zSpec, RoundSpec.roundToEven, 3)
 * // hand-picked positions
 * new FloatToFixedGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.specify("unpack" -> 1, "align" -> 1, "round" -> 1))
 * // all the registers at the module output
 * new FloatToFixedGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.atOut(3))
 * // combinational
 * new FloatToFixedGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.none)
 * }}}
 *
 * @param xSpec spec of the input floating point number.
 * @param yintW width of the second input, fracWidth value. If None, zSpec.fracW is used instead.
 * @param zSpec spec of the output fixed point number.
 * @param roundSpec spec of the rounding
 * @param stage pipeline registers of this module.
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

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = FloatToFixedGeneric.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getParam = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val in  = Flipped(Valid(new FloatToFixedGeneric.Input(xSpec, yintW)))
    val out = Valid(UInt(zSpec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val xIn = cutIn(io.in.bits.x)
  val yIn = io.in.bits.y.map(y => cutIn(y))

  //----------------------------------------------------------------------
  // S1 : decompose / classify
  val (xsgn,  xex,  xmanW1) = FloatChiselUtil.decomposeWithLeading1(xSpec, xIn)
  val (xzero, xinf, xnan)   = FloatChiselUtil.checkValue(xSpec, xIn)

  //--------------------------- stage "unpack" ---------------------------
  private val cutUnpack = new SharedPipe("unpack", cutIn.valid, stages.at("unpack"))

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
  // the result only reads bit xSpec.manW of the shifted value and the bits
  // above it, so the shifter runs on the matching window of zTmp
  val zTmpShifted = zTmp(zSpec.W + xSpec.manW, xSpec.manW) >> zTmpShift

  // everything the shift pushes below that window is the sticky, and it is a
  // mantissa bit: bit i of the mantissa lands below bit xSpec.manW exactly
  // when the shift is larger than i + zSpec.W - xSpec.manW
  val zTmpSticky = (0 to xSpec.manW).map(i => {
    val dropAt = i + zSpec.W - xSpec.manW
    if (dropAt < 0) {
      xmanW1Unpack(i).asBool
    } else {
      xmanW1Unpack(i).asBool && zTmpShift > dropAt.U
    }
  }).reduce(_ || _)

  //--------------------------- stage "align" ----------------------------
  private val cutAlign = new SharedPipe("align", cutUnpack.valid, stages.at("align"))

  val zTmpShiftedAlign = cutAlign(zTmpShifted)
  val zTmpStickyAlign  = cutAlign(zTmpSticky)
  val xsgnAlign        = cutAlign(xsgnUnpack)

  //----------------------------------------------------------------------
  // S3 : rounding increment and negation; in parallel, the exponent range
  //      check (see the special case section below)
  val zTmpLSB   = zTmpShiftedAlign(1).asBool
  val zTmpRound = zTmpShiftedAlign(0).asBool
  val zTmpInc   = FloatChiselUtil.roundIncBySpec(roundSpec, zTmpLSB, zTmpRound,
                                                 zTmpStickyAlign)

  // one adder for the rounding increment and the two's complement:
  // sgn=0 gives w + inc, sgn=1 gives ~w + (1 - inc) == -(w + inc) mod 2^W
  val zTmpWindow = zTmpShiftedAlign(zSpec.W, 1)
  val zSum = (zTmpWindow ^ Fill(zSpec.W, xsgnAlign)) +& (xsgnAlign ^ zTmpInc.asUInt)
  val z    = zSum(zSpec.W-1, 0)
  // a positive value that rounds up past the largest representable one
  val zRoundOver = xsgnAlign === 0.U && (if(zSpec.signed) {
    z(zSpec.W-1).asBool
  } else {
    zSum(zSpec.W).asBool
  })

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

  // the sign bit takes one bit of intW
  val zSgnW = if(zSpec.signed) {
    1
  } else {
    0
  }
  val zExMaxBiased = if(yintW.isEmpty) {
    (xSpec.exBias + zSpec.intW - 1 - zSgnW).U
  } else {
    (xSpec.exBias + zSpec.W - 1 - zSgnW).U - yUnpack.get // intW = totalW - fracW
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
  private val cutRound = new SharedPipe("round", cutAlign.valid, stages.at("round"))

  val zRound          = cutRound(z)
  val zRoundOverRound = cutRound(zRoundOver)
  val xnanRound       = cutRound(xnanAlign)
  val isZPosInfRound  = cutRound(isZPosInfAlign)
  val isZNegInfRound  = cutRound(isZNegInfAlign)
  val isZZeroRound    = cutRound(isZZeroAlign)

  //----------------------------------------------------------------------
  // S4 : Final (special case muxes)
  val z0 = MuxCase(zRound, Seq(
    xnanRound       -> zPosInf,
    isZPosInfRound  -> zPosInf,
    isZNegInfRound  -> zNegInf,
    isZZeroRound    -> 0.U(zSpec.W.W),
    zRoundOverRound -> zPosInf
    ))

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutRound.valid, stages.at("out"))

  io.out.bits  := cutOut(z0)
  io.out.valid := cutOut.valid
}

/** Input bundle and pipeline stage points of [[rial.fixed.FloatToFixedGeneric]].
 *
 * @see [[rial.util.StagePoint]]
 */
object FloatToFixedGeneric {

  /** Input of [[rial.fixed.FloatToFixedGeneric]]: the floating point number
   * and, only when `yW` is given, the fraction width of the result.
   *
   * @param xSpec spec of the input floating point number
   * @param yW    width of the fraction width input, if there is one
   */
  class Input(xSpec: RealSpec, yW: Option[Int]) extends Bundle {
    val x = UInt(xSpec.W.W)
    val y = yW.map(w => UInt(w.W))
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.fixed.FloatToFixedGeneric]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("unpack", 2), StagePoint("align", 3),
    StagePoint("round",  3), StagePoint("out",   1))
}
