//% @file ScaledFixedToFloatGeneric.scala
//
// Fixed point * floating point -> floating point
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

class ScaledFixedToFloatGeneric(
  xSpec : FixedSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  def getParam = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x = Input (UInt(xSpec.W.W))
    val y = Input (UInt(ySpec.W.W))
    val z = Output(UInt(zSpec.W.W))
  })

  val (ysgn,  yex,  ymanW1) = FloatChiselUtil.decomposeWithLeading1(ySpec, io.y)
  val (yzero, yinf, ynan)   = FloatChiselUtil.checkValue(ySpec, io.y)

  val xzero = ~io.x.orR
  val xsgn  = Wire(UInt(1.W))
  val xabs  = Wire(UInt(xSpec.W.W))
  if (xSpec.signed) {
    xsgn := io.x(xSpec.W-1)
    xabs := Mux(io.x(xSpec.W-1), ~io.x + 1.U, io.x)
  } else {
    xsgn := 0.U(1.W)
    xabs := io.x
  }

  val xclz = PriorityEncoder(Reverse(xabs))
  val xShifted = (xabs << xclz)(xSpec.W-1, 0)

  // -------------------------------------------------------------------------
  // mantissa

  val prod = xShifted * ymanW1
  val prodWidth = xSpec.W + ySpec.manW+1
  val prodMoreThan2 = prod(prodWidth - 1)
  assert(prodWidth == prod.getWidth)

  //     decimal point
  //     |
  //     v  xSpec.W + y.manW + 1 = prodWidth
  // .---+------------------------------.
  // | | |          L|R|S.....          |
  //  ^'-------------' '----------------'
  //  1 :   z.manW    :     sticky
  //  : :             :
  //  : :             roundBit = prodWidth - 1 - manW
  //  : prodWidth - 2
  //  prodWidth - 1

  val roundBit = prodWidth - zSpec.manW - 1

  val lsb    = Mux(prodMoreThan2, prod(roundBit),   prod(roundBit-1))
  val round  = Mux(prodMoreThan2, prod(roundBit-1), prod(roundBit-2))
  val sticky = (prodMoreThan2 & prod(roundBit-2)) | prod(roundBit-3, 0).orR

  val man0 = Mux(prodMoreThan2, prod(prodWidth-2, roundBit),
                                prod(prodWidth-3, roundBit-1))
  val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb.asBool, round.asBool, sticky)
  val manRounded = man0 +& inc.asUInt

  val moreThan2AfterRound = manRounded(zSpec.manW)

  // -------------------------------------------------------------------------
  // exponent

  val maxEx = xSpec.intW-1 + ySpec.exMax
  val minEx = -xSpec.fracW + ySpec.exMin

  val zExInc = (prodMoreThan2 | moreThan2AfterRound).asUInt
  val yExNobias = yex.zext - ySpec.exBias.S
  val xExNobias = (xSpec.intW - 1).S - xclz.zext
  val zExNobias = yExNobias +& xExNobias + zExInc.zext
  val zEx0      = zExNobias + zSpec.exBias.S

  val znan  = ynan || (yinf && xzero)
  val zinf  = (zSpec.exMax.S < zExNobias) || yinf
  val zzero = yzero || xzero || zEx0 <= 0.S

  val zEx = Mux(zzero,       0.U(zSpec.exW),
            Mux(zinf | znan, Fill(zSpec.exW, 1.U(1.W)),
                             zEx0(zSpec.exW-1, 0)))

  // -------------------------------------------------------------------------
  // merge

  val zMan = Mux(zzero,       0.U(zSpec.manW),
             Mux(zinf | znan, Cat(znan.asUInt, Fill(zSpec.manW-1, 0.U(1.W))),
                              manRounded(zSpec.manW-1, 0)))

  val z0 = if (zSpec.disableSign) {
    Cat(zEx, zMan)
  } else {
    val zSgn = xsgn ^ ysgn
    Cat(zSgn, zEx, zMan)
  }
//   printf("z0    = %b(w=%d)\n", z0 , z0.getWidth.U)

  io.z := ShiftRegister(z0, nStage)
}
