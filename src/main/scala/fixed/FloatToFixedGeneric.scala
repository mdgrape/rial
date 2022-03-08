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

class FloatToFixedGeneric(
  xSpec : RealSpec, zSpec : FixedSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  def getParam = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x = Input (UInt(xSpec.W.W))
    val z = Output(UInt(zSpec.W.W))
  })

  val (xsgn,  xex,  xmanW1) = FloatChiselUtil.decomposeWithLeading1(xSpec, io.x)
  val (xzero, xinf, xnan)   = FloatChiselUtil.checkValue(xSpec, io.x)

  val zTmpWidth = xSpec.manW + 1 + zSpec.W
  val zTmp      = Cat(xmanW1, Fill(zSpec.W, 0.U(1.W)))
  val zTmpEx    = zSpec.intW - 1

  val zTmpShiftBits = log2Up(zSpec.W)
  val zTmpShift     = ((zSpec.intW-1 + xSpec.exBias).U - xex)(zTmpShiftBits-1, 0)
  val zTmpShifted   = zTmp >> zTmpShift

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

  val zExMax =   zSpec.intW  - 1
  val zExMin = -(zSpec.fracW + 1)
  val zExMaxBiased = xSpec.exBias + zExMax
  val zExMinBiased = max(0, xSpec.exBias + zExMin)

  val isZInf    = xinf  || zExMaxBiased.U < xex
  val isZZero   = xzero || xex < zExMinBiased.U
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
