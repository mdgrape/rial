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
 * @param xSpec spec of the input floating point number.
 * @param yintW width of the second input, fracWidth value. If None, zSpec.fracW is used instead.
 * @param zSpec spec of the output fixed point number.
 * @param roundSpec spec of the rounding
 * @param stage Determines how many registers needed.
 */
class FloatToFixedGeneric(
  xSpec : RealSpec, yintW: Option[Int], zSpec : FixedSpec,
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  def getParam = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x = Input (UInt(xSpec.W.W))
    val y = yintW.map(w => Input(UInt(w.W)))
    val z = Output(UInt(zSpec.W.W))
  })

  val (xsgn,  xex,  xmanW1) = FloatChiselUtil.decomposeWithLeading1(xSpec, io.x)
  val (xzero, xinf, xnan)   = FloatChiselUtil.checkValue(xSpec, io.x)

  val zTmpWidth = xSpec.manW + 1 + zSpec.W
  val zTmp      = Cat(xmanW1, Fill(zSpec.W, 0.U(1.W)))
  val zTmpEx    = if(yintW.isEmpty) {
    zSpec.intW - 1
  } else {
    (zSpec.W - 1).U - io.y.get // intW = totalW - fracW
  }

  val zTmpShiftBits = log2Up(zSpec.W)
  val zTmpShift     = if(yintW.isEmpty) {
    ((zSpec.intW-1 + xSpec.exBias).U - xex)(zTmpShiftBits-1, 0)
  } else {
    ((zSpec.W-1 + xSpec.exBias).U - io.y.get - xex)(zTmpShiftBits-1, 0)
  }
  val zTmpShifted = zTmp >> zTmpShift

  val zTmpLSB    = zTmpShifted(xSpec.manW+1).asBool
  val zTmpRound  = zTmpShifted(xSpec.manW  ).asBool
  val zTmpSticky = zTmpShifted(xSpec.manW-1, 0).orR
  val zTmpInc    = FloatChiselUtil.roundIncBySpec(roundSpec, zTmpLSB, zTmpRound, zTmpSticky)

  val zAbs = zTmpShifted(zSpec.W + xSpec.manW, xSpec.manW+1)
  val z    = Mux(xsgn === 0.U, zAbs, ~zAbs + 1.U)

  // -----------------------------------------------------------------------

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
    (xSpec.exBias + zSpec.W - 1).U - io.y.get // intW = totalW - fracW
  }
  val zExMinBiased = if(io.y.isEmpty) {
    val zExMin = -(zSpec.fracW + 1)
    max(0, xSpec.exBias + zExMin).U
  } else {
    xSpec.exBias.U - io.y.get - 1.U
  }

  val isZInf    = xinf  || zExMaxBiased < xex
  val isZZero   = xzero || xex < zExMinBiased
  val isZPosInf = isZInf && (xsgn === 0.U)
  val isZNegInf = isZInf && (xsgn === 1.U)

  val z0 = MuxCase(z, Seq(
    xnan      -> zPosInf,
    isZPosInf -> zPosInf,
    isZNegInf -> zNegInf,
    isZZero   -> 0.U(zSpec.W.W)
    ))

  io.z := ShiftRegister(z0, nStage)
}
