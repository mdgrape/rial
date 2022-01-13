//% @file exponential.scala
//
// exp function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.FloatChiselUtil

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class ExpPreProcess(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW  = spec.exW
  val manW = spec.manW

  val exBias = spec.exBias

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order
  val extraBits = polySpec.extraBits

  val padding   = extraBits

  val io = IO(new Bundle {
    val x         = Input (UInt(spec.W.W))
    val adr       = Output(UInt(adrW.W))
    val dx        = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }

    // integer part of x, to calculate zex
    val xint      = Output(UInt(exW.W))

    // LSB of fractional part of x, for zman correction
    val xfracLSBs = if(padding != 0) {Some(Output(UInt(padding.W)))} else {None}

    // 1 if (x * log2e).ex is larger than x.ex.
    val xexd = Output(UInt(1.W))
  })

  val (xsgn, xex0, xman0) = FloatChiselUtil.decompose(spec, io.x)

  val log2 = (a:Double) => {log(a) / log(2.0)}
  val xExOvfLimit = math.ceil(log2(maskL(exW)-exBias)).toLong // log2(255-127 = 128) = 7
  val xExUdfLimit = math.ceil(log2(abs(0 - exBias)  )).toLong // log2(|0-127| = 127) > 6

  val xIntW  = max(xExOvfLimit, xExUdfLimit).toInt
  val xFracW = manW + padding

  val extraMan = padding + xIntW

  // --------------------------------------------------------------------------
  // x = x * log2e

  val log2eSpec = new RealSpec(8, 0x7F, 37, false, false, true) // XXX determined empirically
  val log2e = new RealGeneric(log2eSpec, log2(E)) // ~ 1.4427

  val xprod = Cat(1.U(1.W), xman0) * log2e.manW1.toBigInt.U((log2eSpec.manW+1).W)
  val xprodW = xprod.getWidth
  assert(xprodW == (1+manW) + (1+log2eSpec.manW))
  val xprodMoreThan2 = xprod(xprodW - 1) === 1.U

  assert(log2eSpec.manW - extraMan > 0)

  val xprodbp        = spec.manW + log2eSpec.manW
  val xprodRoundBits = xprodbp - (spec.manW + extraMan)

  val xprodSticky  = (xprodMoreThan2 & xprod(xprodRoundBits-1)) |
                     xprod(xprodRoundBits-2, 0).orR
  val xprodRound   = Mux(xprodMoreThan2, xprod(xprodRoundBits), xprod(xprodRoundBits-1))
  val xprodShifted = Mux(xprodMoreThan2, xprod(xprodbp, xprodRoundBits+1),
                                         xprod(xprodbp-1, xprodRoundBits))
  val xprodLSB     = xprodShifted(0)

  val xprodInc = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
    xprodLSB, xprodRound, xprodSticky)

  val xprodRounded = xprodShifted +& xprodInc
  assert(xprodRounded.getWidth == manW+1 + extraMan)
  val xprodMoreThan2AfterRounded = xprodRounded(manW+extraMan)

  val xman = xprodRounded(manW+extraMan-1, 0)
  val xex  = xex0 + xprodMoreThan2AfterRounded + xprodMoreThan2

  assert(xprodMoreThan2AfterRounded +& xprodMoreThan2 =/= 2.U)

  val xexd0 = xprodMoreThan2AfterRounded + xprodMoreThan2
  io.xexd := ShiftRegister(xexd0, nStage)

//   printf("-------------------------\n")
//   printf("PreProc: x     = %b|%b(%d)|%b\n", xsgn, xex0, xex0, xman0)
//   printf("PreProc: xman0 = %b\n", xman0)
//   printf("PreProc: xman  = %b\n", xman)
//   printf("PreProc: xex0  = %b\n", xex0)
//   printf("PreProc: xex   = %b\n", xex)
//   printf("PreProc: xexd  = %b\n", xexd0)

  // --------------------------------------------------------------------------
  // do the same thing as pow2

  val xValW = xIntW + xFracW + 1
  val xVal  = Cat(1.U(1.W), xman)

  val xshift0 = (xIntW + exBias).U(exW.W) - xex
  val xshift  = xshift0(xIntW-1, 0)
  val xValShifted = xVal >> (xshift - 1.U) // XXX
  val xValRounded = (xValShifted >> 1.U) + xValShifted(0)
  // xVal becomes 0 if xint is large enough.
  // To check zinf correctly, we need to output (x*log2e).ex and pass it to OtherPath.
//   printf("PreProc: xVal = %b\n", xVal)
//   printf("PreProc: xRounded = %b\n", xValRounded)

  val xint0  = xValRounded(xIntW+xFracW-1, xFracW)
  val xfrac0 = xValRounded(xFracW-1, 0)
//   printf("PreProc: xint0 = %b\n", xint0)
//   printf("PreProc: xfrac0= %b\n", xfrac0)

  // if xsgn == 1, we negate xint to calculate exbias. but we later do that in
  // ExpOtherPath.
  // To avoid overflow, we need to extend the width here. Since xintW is smaller
  // than exW, this operation is safe. TODO: check if any parameter breaks this
  // relationship.
  val xint = xint0 +& (xsgn.asBool && (xfrac0 =/= 0.U)).asUInt
  io.xint := ShiftRegister(xint(exW-1, 0), nStage)

  val xfracNeg = (1<<xFracW).U((xFracW+1).W) - xfrac0
  val xfrac = Mux(xsgn === 0.U, xfrac0, xfracNeg(xFracW-1, 0))
//   printf("PreProc: xint  = %b\n", xint)
//   printf("PreProc: xfrac = %b\n", xfrac)

  val adr0 = xfrac(xFracW-1, (xFracW-1)-adrW+1)
  io.adr := ShiftRegister(adr0, nStage)

  if(order != 0) {
    val dx0  = Cat(~xfrac(xFracW-1-adrW), xfrac(xFracW-1-adrW-1, padding))
    io.dx.get := ShiftRegister(dx0, nStage)
  }

  if(padding != 0) {
    val xfracLSBs0 = xfrac(padding-1, 0)
    io.xfracLSBs.get := ShiftRegister(xfracLSBs0, nStage)
  }
}
