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
 *     val x = Input (UInt(xSpec.W.W))
 *     val y = Input (UInt(w.W)) // defined only if `yintW` is defined
 *     val z = Output(UInt(zSpec.W.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * First input is the input fixed point number.
 * Second input (optional) is the current fraction width.
 * If yintW is not given (be `None`), it uses xSpec.fracW that will be hard-coded.
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
 *  - `"in"`    the module input
 *  - `"abs"`   after the sign extraction and the absolute value
 *  - `"norm"`  after the leading zero count and the normalisation shift
 *  - `"round"` after the rounding increment
 *  - `"out"`   the module output
 *
 * {{{
 * // three registers, placed by the module
 * new FixedToFloatGeneric(xSpec, None, zSpec, RoundSpec.roundToEven, 3)
 * // hand-picked positions
 * new FixedToFloatGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.specify("abs" -> 1, "norm" -> 1, "round" -> 1))
 * // everything at the output, as before
 * new FixedToFloatGeneric(xSpec, None, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.atOut(3))
 * }}}
 *
 * @param xSpec spec of the input fixed point number.
 * @param yintW width of the second input, fracWidth value. If None, xSpec.fracW is used instead.
 * @param zSpec spec of the output floating point number.
 * @param roundSpec spec of the rounding
 * @param stage Determines how many registers needed.
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

  /** Boundaries of the datapath at which pipeline registers can be placed.
   * The weights are rough, target-independent estimates of the delay of the
   * logic that ends at each boundary. */
  val stagePoints = FixedToFloatGeneric.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getParam = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x = Input (UInt(xSpec.W.W))
    val y = yintW.map(w => Input (UInt(w.W)))
    val z = Output(UInt(zSpec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  val xIn = ShiftRegister(io.x, stages.at("in"))
  val yIn = io.y.map(y => ShiftRegister(y, stages.at("in")))

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
  private def cutAbs[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("abs"))
  }

  val xabsAbs  = cutAbs(xabs)
  val zSgnAbs  = cutAbs(zSgn)
  val xzeroAbs = cutAbs(xzero)
  val yAbs     = yIn.map(y => cutAbs(y))

  //----------------------------------------------------------------------
  // S2 : leading zero count and normalisation shift

  val xclz = PriorityEncoder(Reverse(xabsAbs))
  val xShifted = (xabsAbs << xclz)(xSpec.W-1, 0)

  //---------------------------- stage "norm" ----------------------------
  private def cutNorm[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("norm"))
  }

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

  } else if(xSpec.W == zSpec.manW+1) {
    // no padding needed
    zman0 := xShiftedNorm(zSpec.manW-1, 0) // remove the leading 1

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
  private def cutRound[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("round"))
  }

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

  io.z := ShiftRegister(z0, stages.at("out"))
}

/** Pipeline stage points of [[rial.fixed.FixedToFloatGeneric]].
 *
 * @see [[rial.util.StagePoint]]
 */
object FixedToFloatGeneric {
  /** Ordered boundaries of the datapath, with a rough, target-independent
   *  estimate of the delay of the logic ending at each. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("abs",   2), StagePoint("norm", 3),
    StagePoint("round", 2), StagePoint("out",  2))
}
