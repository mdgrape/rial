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
 * First input is the input fixed point number.
 * Second input (optional) is the current fraction width.
 * If yintW is not given (be `None`), it uses xSpec.fracW that will be hard-coded.
 *
 * 0 means the fixed point value is an integer (w/o fractional part).
 * Since, in most cases, fixed point is signed, it results a weird result if
 * fracW > totalW in a fixed point value.
 * Take care of the width of the second input, y.
 *
 * @param xSpec spec of the input fixed point number.
 * @param yintW width of the second input, fracWidth value. If None, xSpec.fracW is used instead.
 * @param zSpec spec of the output floating point number.
 * @param roundSpec spec of the rounding
 * @param stage Determines how many registers needed.
 */
class FixedToFloatGeneric(
  xSpec : FixedSpec, yintW : Option[Int], zSpec : RealSpec,
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  def getParam = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x = Input (UInt(xSpec.W.W))
    val y = yintW.map(w => Input (UInt(w.W)))
    val z = Output(UInt(zSpec.W.W))
  })

  // extract sign and make x absolute value

  val xzero = ~io.x.orR
  val zSgn  = Wire(UInt(1.W))
  val xabs  = Wire(UInt(xSpec.W.W))
  if (xSpec.signed) {
    zSgn := io.x(xSpec.W-1)
    xabs := Mux(io.x(xSpec.W-1), ~io.x + 1.U, io.x)
  } else {
    zSgn := 0.U(1.W)
    xabs := io.x
  }

  val xclz = PriorityEncoder(Reverse(xabs))
  val xShifted = (xabs << xclz)(xSpec.W-1, 0)

  val zexInc = Wire(UInt(1.W))
  val zman0 = Wire(UInt(zSpec.manW.W))
  if(xSpec.W < zSpec.manW+1) {
    // need to add least bits (0)
    val padding = zSpec.manW+1 - xSpec.W
    assert(0 < padding)
    zman0 := Cat(xShifted(xSpec.W-2, 0), Fill(padding, 0.U(1.W)))

  } else if(xSpec.W == zSpec.manW+1) {
    // no padding needed
    zman0 := xShifted(zSpec.manW-1, 0) // remove the leading 1

  } else if(xSpec.W == zSpec.manW+2) {
    // rounding, but sticky bit is always 0

    val sticky = false.B
    val round  = xShifted(0).asBool
    val lsb    = xShifted(1).asBool
    val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky)

    val manRounded = xShifted(xSpec.W-2, 1) +& inc.asUInt
    zexInc := manRounded(zSpec.manW)
    zman0  := manRounded(zSpec.manW-1, 0)

  } else {
    assert(xSpec.W > zSpec.manW+2)
    // need to round.
    val lsbBit   = xSpec.W - zSpec.manW - 1
    val roundBit = xSpec.W - zSpec.manW - 2
    assert(0 < lsbBit  )
    assert(0 < roundBit)
    val lsb      = xShifted(lsbBit       ).asBool
    val round    = xShifted(roundBit     ).asBool
    val sticky   = xShifted(roundBit-1, 0).orR
    val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky)

    val manRounded = xShifted(xSpec.W-2, xSpec.W-zSpec.manW-1) +& inc.asUInt
    zexInc := manRounded(zSpec.manW)
    zman0  := manRounded(zSpec.manW-1, 0)
  }
  val zex0 = if (io.y.isEmpty) {
    (zSpec.exBias + xSpec.W - 1 - xSpec.fracW).U - xclz + zexInc
  } else {
    (zSpec.exBias + xSpec.W - 1).U - io.y.get - xclz + zexInc
  }

//   printf("zman0 = %b\n", zman0)
//   printf("zex0  = %b\n", zex0 )

  val zEx  = Wire(UInt(zSpec.exW.W))
  val zMan = Wire(UInt(zSpec.manW.W))

  if(io.y.isEmpty && xSpec.intW < zSpec.exMax) { // never be inf.
    zEx  := Mux(xzero, 0.U(zSpec.exW.W),  zex0)
    zMan := Mux(xzero, 0.U(zSpec.manW.W), zman0)
  } else {
    val zinf = maskI(zSpec.exW).U < zex0 // zex0 is biased, don't use exMax
    zEx  := Mux(zinf, zSpec.exMax.U(zSpec.exW.W), Mux(xzero, 0.U(zSpec.exW.W),  zex0))
    zMan := Mux(xzero || zinf, 0.U(zSpec.manW.W), zman0)
  }
//   printf("zEx   = %b\n", zEx  )
//   printf("zMan  = %b\n", zMan )

  val z0 = if (zSpec.disableSign) {Cat(zEx, zMan)} else {Cat(zSgn, zEx, zMan)}
//   printf("z0    = %b(w=%d)\n", z0 , z0.getWidth.U)

  io.z := ShiftRegister(z0, nStage)
}
