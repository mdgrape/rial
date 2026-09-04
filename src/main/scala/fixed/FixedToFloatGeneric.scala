//% @file FixedToFloatGeneric.scala
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

/** Fixed point number -> floating point number.
 *
 * {{{
 * class FixedToFloatGeneric(...) extends Module {
 *   val io = IO(new Bundle{
 *     val in  = Flipped(Valid(new FixedToFloatGeneric.Input(xSpec, yintW)))
 *     val out = Valid(UInt(zSpec.W.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * First input, `in.bits.x`, is the input fixed point number.
 * Second input (optional), `in.bits.y`, is the current fraction width. It is
 * defined only if `yintW` is defined.
 * If yintW is not given (be `None`), it uses xSpec.fracW that will be hard-coded.
 *
 * 0 means the fixed point value is an integer (w/o fractional part).
 * Since, in most cases, fixed point is signed, it results a weird result if
 * fracW > totalW in a fixed point value.
 * Take care of the width of the second input, y.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`    the module input
 *  - `"abs"`   after the sign extraction and the absolute value
 *  - `"norm"`  after the leading zero count and the normalisation shift
 *  - `"round"` after the rounding increment
 *  - `"out"`   the module output
 *
 * [[getStage]] gives the resulting latency, and `out.valid` is `in.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // three registers, placed by the module
 * new FixedToFloatGeneric(xSpec, None, zSpec, RoundSpec.roundToEven, 3)
 * // hand-picked positions
 * new FixedToFloatGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.specify("abs" -> 1, "norm" -> 1, "round" -> 1))
 * // all the registers at the module output
 * new FixedToFloatGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.atOut(3))
 * // combinational
 * new FixedToFloatGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.none)
 * }}}
 *
 * @param xSpec spec of the input fixed point number.
 * @param yintW width of the second input, fracWidth value. If None, xSpec.fracW is used instead.
 * @param zSpec spec of the output floating point number.
 * @param roundSpec spec of the rounding
 * @param stage pipeline registers of this module.
 *
 * For parameters, see also:
 * @see [[rial.fixed.FixedSpec]]
 * @see [[rial.arith.RealSpec]]
 * @see [[rial.arith.RoundSpec]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class FixedToFloatGeneric(
  xSpec : FixedSpec, yintW : Option[Int], zSpec : RealSpec,
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = FixedToFloatGeneric.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getParam = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val in  = Flipped(Valid(new FixedToFloatGeneric.Input(xSpec, yintW)))
    val out = Valid(UInt(zSpec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val xIn = cutIn(io.in.bits.x)
  val yIn = io.in.bits.y.map(y => cutIn(y))

  //----------------------------------------------------------------------
  // S1 : extract sign and make x absolute value

  val xzero = ~xIn.orR
  val zSgn  = Wire(UInt(1.W))
  val xabs  = Wire(UInt(xSpec.W.W))
  if (xSpec.signed) {
    zSgn := xIn(xSpec.W-1)
    xabs := Mux(xIn(xSpec.W-1), ~xIn + 1.U, xIn)
  } else {
    zSgn := 0.U(1.W)
    xabs := xIn
  }

  //---------------------------- stage "abs" -----------------------------
  private val cutAbs = new SharedPipe("abs", cutIn.valid, stages.at("abs"))

  val xabsAbs  = cutAbs(xabs)
  val zSgnAbs  = cutAbs(zSgn)
  val xzeroAbs = cutAbs(xzero)
  val yAbs     = yIn.map(y => cutAbs(y))

  //----------------------------------------------------------------------
  // S2 : leading zero count and normalisation shift

  val xclz = PriorityEncoder(Reverse(xabsAbs))
  val xShifted = (xabsAbs << xclz)(xSpec.W-1, 0)

  //---------------------------- stage "norm" ----------------------------
  private val cutNorm = new SharedPipe("norm", cutAbs.valid, stages.at("norm"))

  val xShiftedNorm = cutNorm(xShifted)
  val xclzNorm     = cutNorm(xclz)
  val zSgnNorm     = cutNorm(zSgnAbs)
  val xzeroNorm    = cutNorm(xzeroAbs)
  val yNorm        = yAbs.map(y => cutNorm(y))

  //----------------------------------------------------------------------
  // S3 : mantissa selection and rounding

  val zexInc = Wire(UInt(1.W))
  val zman0 = Wire(UInt(zSpec.manW.W))
  if(xSpec.W < zSpec.manW+1) {
    // need to add least bits (0)
    val padding = zSpec.manW+1 - xSpec.W
    assert(0 < padding)
    zman0 := Cat(xShiftedNorm(xSpec.W-2, 0), Fill(padding, 0.U(1.W)))
    zexInc := 0.U // no rounding, so no carry into the exponent

  } else if(xSpec.W == zSpec.manW+1) {
    // no padding needed
    zman0 := xShiftedNorm(zSpec.manW-1, 0) // remove the leading 1
    zexInc := 0.U // no rounding, so no carry into the exponent

  } else if(xSpec.W == zSpec.manW+2) {
    // rounding, but sticky bit is always 0

    val sticky = false.B
    val round  = xShiftedNorm(0).asBool
    val lsb    = xShiftedNorm(1).asBool
    val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky)

    val manRounded = xShiftedNorm(xSpec.W-2, 1) +& inc.asUInt
    zexInc := manRounded(zSpec.manW)
    zman0  := manRounded(zSpec.manW-1, 0)

  } else {
    assert(xSpec.W > zSpec.manW+2)
    // need to round.
    val lsbBit   = xSpec.W - zSpec.manW - 1
    val roundBit = xSpec.W - zSpec.manW - 2
    assert(0 < lsbBit  )
    assert(0 < roundBit)
    val lsb      = xShiftedNorm(lsbBit       ).asBool
    val round    = xShiftedNorm(roundBit     ).asBool
    val sticky   = xShiftedNorm(roundBit-1, 0).orR
    val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky)

    val manRounded = xShiftedNorm(xSpec.W-2, xSpec.W-zSpec.manW-1) +& inc.asUInt
    zexInc := manRounded(zSpec.manW)
    zman0  := manRounded(zSpec.manW-1, 0)
  }

  //--------------------------- stage "round" ----------------------------
  private val cutRound = new SharedPipe("round", cutNorm.valid, stages.at("round"))

  val zman0Round  = cutRound(zman0)
  val zexIncRound = cutRound(zexInc)
  val xclzRound   = cutRound(xclzNorm)
  val zSgnRound   = cutRound(zSgnNorm)
  val xzeroRound  = cutRound(xzeroNorm)
  val yRound      = yNorm.map(y => cutRound(y))

  //----------------------------------------------------------------------
  // S4 : exponent and the special case muxes

  val zex0 = if (yRound.isEmpty) {
    (zSpec.exBias + xSpec.W - 1 - xSpec.fracW).U - xclzRound + zexIncRound
  } else {
    (zSpec.exBias + xSpec.W - 1).U - yRound.get - xclzRound + zexIncRound
  }

//   printf("zman0 = %b\n", zman0Round)
//   printf("zex0  = %b\n", zex0 )

  val zEx  = Wire(UInt(zSpec.exW.W))
  val zMan = Wire(UInt(zSpec.manW.W))

  if(yRound.isEmpty && xSpec.intW < zSpec.exMax) { // never be inf.
    zEx  := Mux(xzeroRound, 0.U(zSpec.exW.W),  zex0)
    zMan := Mux(xzeroRound, 0.U(zSpec.manW.W), zman0Round)
  } else {
    val zinf = maskI(zSpec.exW).U < zex0 // zex0 is biased, don't use exMax
    zEx  := Mux(zinf, zSpec.exMax.U(zSpec.exW.W), Mux(xzeroRound, 0.U(zSpec.exW.W),  zex0))
    zMan := Mux(xzeroRound || zinf, 0.U(zSpec.manW.W), zman0Round)
  }
//   printf("zEx   = %b\n", zEx  )
//   printf("zMan  = %b\n", zMan )

  val z0 = if (zSpec.disableSign) {Cat(zEx, zMan)} else {Cat(zSgnRound, zEx, zMan)}
//   printf("z0    = %b(w=%d)\n", z0 , z0.getWidth.U)

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutRound.valid, stages.at("out"))

  io.out.bits  := cutOut(z0)
  io.out.valid := cutOut.valid
}

/** Input bundle and pipeline stage points of [[rial.fixed.FixedToFloatGeneric]].
 *
 * @see [[rial.util.StagePoint]]
 */
object FixedToFloatGeneric {

  /** Input of [[rial.fixed.FixedToFloatGeneric]]: the fixed point number and,
   * only when `yW` is given, its fraction width.
   *
   * @param xSpec spec of the input fixed point number
   * @param yW    width of the fraction width input, if there is one
   */
  class Input(xSpec: FixedSpec, yW: Option[Int]) extends Bundle {
    val x = UInt(xSpec.W.W)
    val y = yW.map(w => UInt(w.W))
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.fixed.FixedToFloatGeneric]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("abs",   2), StagePoint("norm", 3),
    StagePoint("round", 2), StagePoint("out",  2))
}
