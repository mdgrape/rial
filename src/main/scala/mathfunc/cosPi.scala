//% @file sinPi.scala
//
// x -> sin(pi * x)
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import chisel3._
import chisel3.util._

import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.FloatChiselUtil
import rial.arith._

import rial.math.SinPiSim
import rial.mathfunc._

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class CosPiPreProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias
  val exAdrW = SinPiSim.calcExAdrW(spec)

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val io = IO(new Bundle {
    val x   = Input (UInt(spec.W.W))
    val adr = Output(UInt((exAdrW+adrW).W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
    val xConverted = new SinPiPreProcessOutput(spec)
  })

  val (xsgn, xex, xman) = FloatChiselUtil.decompose(spec,  io.x)

  // --------------------------------------------------------------------------
  // convert cos to sin

  val from0_ex0_neg = xman - (1 << (manW-1)).U
  val from0_ex0_pos = (1 << (manW-1)).U - xman
  val from0_ex0     = Mux(xman(manW-1), from0_ex0_neg, from0_ex0_pos)
  val shift_ex0     = PriorityEncoder(Reverse(from0_ex0)) + 1.U
  val norm_ex0      = from0_ex0 << shift_ex0

  val from0_ex1 = xman
  val shift_ex1 = PriorityEncoder(Reverse(from0_ex1)) + 1.U
  val norm_ex1  = from0_ex1 << shift_ex1

  val align_exmanW = ((1 << manW) - 1 + exBias).U - xex
  val from0_exmanW = ((1 << manW).U - (((1<<manW).U + xman) >> align_exmanW))(manW-1, 0)
  val shift_exmanW = PriorityEncoder(Reverse(from0_exmanW)) + 1.U
  val norm_exmanW  = from0_exmanW << shift_exmanW

  val yman = MuxCase(((1 << manW) - 1).U(manW.W), Seq(
    (xex === (exBias    ).U)                      -> norm_ex0(manW-1, 0),
    (xex === (exBias - 1).U)                      -> norm_ex1(manW-1, 0),
    ((exBias-manW).U < xex && xex < (exBias-1).U) -> norm_exmanW(manW-1, 0)
  ))
  val yex  = MuxCase((exBias-2).U(exW.W), Seq(
    (xex === (exBias    ).U)                      -> Mux(from0_ex0    === 0.U, 0.U, exBias.U     - shift_ex0   ),
    (xex === (exBias - 1).U)                      -> Mux(from0_ex1    === 0.U, 0.U, (exBias-1).U - shift_ex1   ),
    ((exBias-manW).U < xex && xex < (exBias-1).U) -> Mux(from0_exmanW === 0.U, 0.U, (exBias-1).U - shift_exmanW)
  ))
  assert(shift_ex0    <= exBias.U)
  assert(shift_ex1    <= (exBias-1).U)
  assert(shift_exmanW <= (exBias-1).U)

  val exOfs = (exBias-2).U - yex
  val exAdr = exOfs.asUInt()(exAdrW-1, 0)

  val adr0 = Cat(exAdr, yman(manW-1, manW-adrW)) // concat exAdr + man
  val dr0  = Cat(~yman(manW-adrW-1), yman(manW-adrW-2,0))

  io.adr   := ShiftRegister(adr0,  nStage)
  if(order != 0) {
    io.dx.get := ShiftRegister(dr0, nStage)
  }
  io.xConverted.xConvertedEx  := ShiftRegister(yex,  nStage)
  io.xConverted.xConvertedMan := ShiftRegister(yman, nStage)
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

// The same as sin.

// -------------------------------------------------------------------------
//                        _        _     _                    _   _
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_|
//                                               |_|
// -------------------------------------------------------------------------

// Some edge cases and signs are different from sinPi.

class CosPiOtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val io = IO(new Bundle {
    val x          = Input(UInt(spec.W.W))
    val xConverted = Flipped(new SinPiPreProcessOutput(spec))
    val zother     = new SinPiNonTableOutput(spec)
  })

  val (xsgn,  xex,  xman) = FloatChiselUtil.decompose(spec, io.x)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(spec, io.x)

  val yex  = io.xConverted.xConvertedEx
  val yman = io.xConverted.xConvertedMan

//   val out_of_bounds = ((exBias+1).U <= xex) && (xman =/= 0.U) // 2 < |x|
//   val znan  = xnan || xinf || out_of_bounds
//   val zzero = ((exBias.U <= xex)    && xman === 0.U) ||
//               (yex === 0.U          && yman === 0.U)
//   val zone  = (yex === (exBias-1).U && yman === 0.U)
//
//   val zSgn  = (!zzero && !znan) && ((xsgn === 1.U) || (xex === exBias.U))

  val out_of_bounds = ((exBias+1).U === xex && xman =/= 0.U) || ((exBias+1).U < xex)
  val znan  = xnan || xinf || out_of_bounds

  val zzero = (xex === (exBias-1).U && xman === 0.U)               || // x == 0.5
              (xex === (exBias  ).U && xman === (1 << (manW-1)).U) || // x == 1.5
              (yex === 0.U          && yman === 0.U)                  // y == 0
  val zone  =  xzero                                 || // x == 0
              (xex === (exBias  ).U && xman === 0.U) || // x == 1
              (xex === (exBias+1).U && xman === 0.U) || // x == 2
              (yex === (exBias-1).U && yman === 0.U)    // y == 0.5

  val zSgn = MuxCase(0.U(1.W), Seq(
    (xex === (exBias  ).U  && !zzero)             -> ~(xman(manW-1)),
    (xex === (exBias-1).U  && !zzero)             -> 1.U(1.W),
    ((exBias-manW).U < xex && xex < (exBias-1).U) -> 0.U(1.W)
  ))

  // --------------------------------------------------------------------------
  // linear approximation around zero

  val linearThreshold = SinPiSim.calcLinearThreshold(manW)
  val pi = new RealGeneric(spec, Pi)

  val isLinear = (yex < (linearThreshold + exBias).U)

  // y is in [0, 0.5) so it never overflows
  val linearProdEx        = (pi.ex-exBias).U(exW.W) + yex
  val linearProdMan       = (pi.man + (1<<manW)).toLong.U((manW+1).W) *
                            (yman   + (1<<manW).U((manW+1).W))
  val linearProdbp        = manW + manW
  val linearProdMoreThan2 = linearProdMan(linearProdbp+1)
  val linearProdRoundBits = linearProdbp - manW

  val linearProdSticky    = linearProdMan(linearProdRoundBits-2, 0).orR |
                           (linearProdMoreThan2 & linearProdMan(linearProdRoundBits-1))
  val linearProdRound     = Mux(linearProdMoreThan2, linearProdMan(linearProdRoundBits),
                                                     linearProdMan(linearProdRoundBits-1))
  val linearProdShift     = Mux(linearProdMoreThan2, linearProdMan(linearProdbp, linearProdRoundBits+1),
                                                     linearProdMan(linearProdbp-1, linearProdRoundBits))
  val linearProdLSB       = linearProdShift(0)

  val linearProdInc       = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
                                linearProdLSB, linearProdRound, linearProdSticky)
  val linearProdRounded   = linearProdShift +& linearProdInc
  val linearProdMoreThan2AfterRound = linearProdRounded(manW)

  val zManLinear = linearProdRounded(manW-1, 0)
  val zExLinear  = linearProdEx + (linearProdMoreThan2 | linearProdMoreThan2AfterRound)

  // TODO: cubic

  val zIsNonTable = znan || zzero || zone || isLinear
  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)

  val zEx = Mux(isLinear, zExLinear,
            Mux(znan,     Fill(exW, 1.U(1.W)),
            Mux(zzero,    0.U,
            Mux(zone,     exBias.U,
                          yex))))   // default

  val zeroFlush = znan || zzero || zone
  val zMan = Mux(zeroFlush, Cat(znan, 0.U((manW-1).W)), zManLinear)

  io.zother.zsgn := ShiftRegister(zSgn, nStage)
  io.zother.zex  := ShiftRegister(zEx,  nStage)
  io.zother.zman := ShiftRegister(zMan, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

// The same as sin.

