//% @file acos1.scala
//
// ACos function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._

import spire.math.SafeLong
import spire.math.Real
import spire.implicits._

import chisel3._
import chisel3.util._
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith._

// Calculate acos(x) if x is in range [-1, 1].
// otherwise, 0 or pi for positive and negative x, respectively.
//
// Here, we calculate acos(x) by using the following two phases.
//
// 1. x -> sqrt(1 - |x|)
// 2. x -> acos(1 - x^2)
//
// note that: `acos(1 - sqrt(1-|x|)^2) = acos(1-(1-|x|)) = acos(|x|)`.
//
// So, we have additional preprocess stage that calculates 1-|x|.
// Phase1 uses Phase1Pre/PostProcess and Sqrt table.
// Phase2 uses the specific table that calculates acos(1 - x^2)
//
// This file defines PreProcess and PostProcess of acos phase 1.

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

private[rial] object ACosSpecialValue {
  val width     = 2
  val xNaNInf   = 0.U(width.W) // x == nan/inf, acos(x) == nan
  val xZero     = 1.U(width.W) // x ==  0,      acos(x) == pi/2
  val xLargePos = 2.U(width.W) // x >=  1,      acos(x) == +0
  val xLargeNeg = 3.U(width.W) // x <= -1,      acos(x) == pi
}

private[rial] class ACosFlags extends Bundle {
  val isSpecial = Bool()
  val xsgn      = UInt(1.W)
  val flag      = UInt(ACosSpecialValue.width.W)
}

private[rial] class ACosPhase1PreProcess(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val dxW    = polySpec.dxW
  val order  = polySpec.order
  val exAdrW = 1 // sqrt extra table address

  val io = IO(new Bundle {
    val en  = Input (UInt(1.W))
    val x   = Flipped(new DecomposedRealOutput(spec))

    // 1-|x| for sqrt table
    val adr = Output(UInt((exAdrW+adrW).W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }

    // for sqrt non-table path
    val y = new DecomposedRealOutput(spec)
    // for postproc (encode special value flags for acos-2)
    val special = Output(new ACosFlags())
  })

  val xsgn = enableIf(io.en, io.x.sgn)
  val xex  = enableIf(io.en, io.x.ex )
  val xman = enableIf(io.en, io.x.man)

  val xLargerThan1 = xex >= exBias.U

  // -------------------------------------------------------------------------
  // calc 1 - |x|

  val xShiftMax = (1+manW+2).U
  val xShift0   = exBias.U(exW.W) - xex
  val xShift    = Mux(xShiftMax < xShift0, xShiftMax, xShift0)

  val xShifted = (Cat(1.U(1.W), xman, 0.U(2.W)) >> xShift)
  val xSubtracted = Mux(xLargerThan1, 0.U, ~(xShifted(manW+2-1, 0)) +& 1.U)

  val xNormalizeShift = PriorityEncoder(Reverse(xSubtracted))
  val xNormalized     = (xSubtracted << xNormalizeShift)(1+manW+2-1, 0)

//   printf("cir: xex            = %d\n", xex)
//   printf("cir: xman           = %b\n", xman)
//   printf("cir: xShift         = %d\n", xShift)
//   printf("cir: xShifted       = %b\n", xShifted)
//   printf("cir: xSub           = %b(W=%d)\n", xSubtracted, xSubtracted.getWidth.U)
//   printf("cir: xNormalizeShift= %d\n", xNormalizeShift)
//   printf("cir: xNormalized    = %b\n", xNormalized)
  assert(xNormalized(1+manW+2-1) === 1.U || xNormalized === 0.U) // normalized or zero

  val xRounded          = xNormalized(manW+1, 2) +& xNormalized(1)
  val xRoundedMoreThan2 = xRounded(manW)
  val xConvertedEx0  = exBias.U(exW.W) - xNormalizeShift + xRoundedMoreThan2
  val xConvertedMan0 = xRounded(manW-1, 0)

  val xSmallEnough = io.x.zero || xConvertedEx0 === exBias.U // -> special value

  val xConvertedEx  = Mux(xSmallEnough, exBias.U(exW.W), Mux(xLargerThan1, 0.U, xConvertedEx0))
  val xConvertedMan = Mux(xSmallEnough, 0.U(manW.W),     Mux(xLargerThan1, 0.U, xConvertedMan0))

  // -------------------------------------------------------------------------
  // do the same thing as sqrt

  val adr = enableIf(io.en, Cat(xConvertedEx(0), xConvertedMan(manW-1, dxW)))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx = enableIf(io.en, Cat(~xConvertedMan(dxW-1), xConvertedMan(dxW-2, 0)))
    io.dx.get := ShiftRegister(dx, nStage)
  }

//   printf("cir: xConvertedEx  = %d\n", xConvertedEx )
//   printf("cir: xConvertedMan = %b\n", xConvertedMan)

  // pass 1-|x| to sqrtOtherPath
  io.y.sgn  := ShiftRegister(enableIf(io.en, 0.U(1.W)),      nStage)
  io.y.ex   := ShiftRegister(enableIf(io.en, xConvertedEx),  nStage)
  io.y.man  := ShiftRegister(enableIf(io.en, xConvertedMan), nStage)
  io.y.zero := ShiftRegister(enableIf(io.en, xLargerThan1),  nStage)
  io.y.inf  := ShiftRegister(enableIf(io.en, io.x.inf),      nStage)
  io.y.nan  := ShiftRegister(enableIf(io.en, io.x.nan),      nStage)

  // -------------------------------------------------------------------------
  // check if special value

  val isSpecial = io.x.inf || io.x.nan || xSmallEnough || xLargerThan1
  val flag = MuxCase(0.U(ACosSpecialValue.width.W), Array(
    (io.x.inf || io.x.nan)         -> ACosSpecialValue.xNaNInf,
    (xSmallEnough)                 -> ACosSpecialValue.xZero,
    (xLargerThan1 && xsgn === 1.U) -> ACosSpecialValue.xLargeNeg,
    (xLargerThan1 && xsgn === 0.U) -> ACosSpecialValue.xLargePos
  ))

  io.special.isSpecial := ShiftRegister(enableIf(io.en, isSpecial), nStage)
  io.special.flag      := ShiftRegister(enableIf(io.en, flag),      nStage)
  io.special.xsgn      := ShiftRegister(enableIf(io.en, xsgn),      nStage)
}

// -------------------------------------------------------------------------
// use sqrt coeff table, sqrt non-table path

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

private[rial] class ACosPhase1PostProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val nStage = stage.total
  def getStage = nStage

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val order     = polySpec.order
  val extraBits = polySpec.extraBits

  val io = IO(new Bundle {
    val en     = Input(UInt(1.W))
    val zother = Input(new RoundingNonTableOutput(spec)) // from sqrt non-table
    val zres   = Input(UInt(fracW.W))                    // from table
    val flags  = Input(new ACosFlags)                    // from acos
    val z      = Output(UInt(spec.W.W))
  })

  val zsgn = WireDefault(0.U(1.W))
  val zex  = WireDefault(0.U(exW.W))
  val zman = WireDefault(0.U(manW.W))

  val z = enableIf(io.en, Cat(zsgn, zex, zman))
  io.z := ShiftRegister(z, nStage)

  // --------------------------------------------------------------------------
  // normal sqrt result

  val zexSqrt  = io.zother.zex

  val zmanSqrt = WireDefault(0.U(manW.W))
  if(extraBits == 0) {
    zmanSqrt := io.zres
  } else {
    val zman0 = dropLSB(extraBits, io.zres) +& io.zres(extraBits-1)
    val polynomialOvf = zman0(manW)
    zmanSqrt := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW-1,0))
  }

  // --------------------------------------------------------------------------
  // encode special values

  val zsgnSpecial = io.flags.xsgn
  val zexSpecial  = Fill(exW, 1.U(1.W)) // use nan-boxing
  val zmanSpecial = Mux(io.flags.flag === ACosSpecialValue.xNaNInf, Fill(manW, 1.U(1.W)), io.flags.flag)

  zsgn := zsgnSpecial // encode x sign (sqrt is always positive)
  zex  := Mux(io.flags.isSpecial,  zexSpecial,  zexSqrt)
  zman := Mux(io.flags.isSpecial, zmanSpecial, zmanSqrt)
}

