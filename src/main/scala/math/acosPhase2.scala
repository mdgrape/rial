//% @file acos.scala
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

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

private[rial] class ACosPhase2PreProcess(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
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
    val en  = Input (UInt(1.W))
    val x   = Flipped(new DecomposedRealOutput(spec))
    val adr = Output(UInt(adrW.W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }

    // decoded flags
    val flags = Output(new ACosFlags)
  })

  val xex  = enableIf(io.en, io.x.ex )
  val xman = enableIf(io.en, io.x.man)

  val xShift = exBias.U - xex
  val xAligned = Cat(1.U(1.W), xman) >> xShift

  val adr = enableIf(io.en, xAligned(manW-1, dxW))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx = enableIf(io.en, Cat(~xAligned(dxW-1), xAligned(dxW-2, 0)))
    io.dx.get := ShiftRegister(dx, nStage)
  }

  val flags = WireDefault(0.U.asTypeOf(new ACosFlags))
  flags.isSpecial := io.x.ex === Fill(exW, 1.U(1.W))
  flags.xsgn      := io.x.sgn
  flags.flag      := Mux(io.x.man === Fill(manW, 1.U(1.W)), ACosSpecialValue.xNaNInf, io.x.man)

  io.flags := ShiftRegister(flags, nStage)
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

private[rial] class ACosTableCoeff(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
) extends Module {

  val manW   = spec.manW
  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input(UInt((1+adrW).W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {
    val tableI = ACosPhase2Sim.acosTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    // sign mode 1 = 2's complement and no sign bit
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/1)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    assert(maxCbit(0) == fracW)
    assert(coeff(0).getWidth == fracW)

    io.cs.cs(0) := enableIf(io.en, coeff(0))

  } else {
    val tableI = ACosPhase2Sim.acosTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(0 <= diffWidth)

      if(diffWidth != 0) {
        val ci  = coeff(i)
        val msb = ci(cbit(i)-1)
        coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffs.cs(i) := coeff(i)
      }
    }
    io.cs := enableIf(io.en, coeffs)
  }
}
private[rial] object ACosTableCoeff {
  def getCBits(
    spec:     RealSpec,
    polySpec: PolynomialSpec
  ): Seq[Int] = {

    val order     = polySpec.order
    val adrW      = polySpec.adrW
    val extraBits = polySpec.extraBits
    val fracW     = polySpec.fracW

    if(order == 0) {
      return Seq(fracW)
    } else {
      return ACosPhase2Sim.acosTableGeneration( order, adrW, spec.manW, fracW ).cbit
    }
  }
  def getCalcW(
    spec:     RealSpec,
    polySpec: PolynomialSpec
  ): Seq[Int] = {

    val order     = polySpec.order
    val adrW      = polySpec.adrW
    val extraBits = polySpec.extraBits
    val fracW     = polySpec.fracW

    if(order == 0) {
      return Seq(fracW)
    } else {
      return ACosPhase2Sim.acosTableGeneration( order, adrW, spec.manW, fracW ).calcWidth
    }
  }
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

private[rial] class ACosPhase2PostProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val nStage = stage.total
  def getStage = nStage

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order
  val extraBits = polySpec.extraBits

  val realPi      = new RealGeneric(spec, Pi)
  val realHalfPi  = new RealGeneric(spec, Pi * 0.5)

  val io = IO(new Bundle {
    val en    = Input(Bool())
    val flags = Input(new ACosFlags())
    val zex0  = Input(UInt(exW.W))
    val zman0 = Input(UInt(manW.W))
    val exInc = Input(UInt(1.W))
    val z     = Output(UInt(spec.W.W))
  })

  // -----------------------------------------------------------------------
  // calc acos(|x|)

  val zMan0  = io.zman0(manW-1, 0)
  val zEx0   = io.zex0 + io.exInc // table result has ex == exBias

  // -----------------------------------------------------------------------
  // calc pi - z in case of x < 0

  val piW = 2+manW+2 // 2 for int part, manW+2 after the decimal point
  val pi = (realPi.manW1 << 3).toBigInt.U(piW.W)

  val zAligned = Cat(1.U(2.W), zMan0, 0.U(2.W)) >> (exBias.U - zEx0)
  val zNeg     = pi - zAligned
  val zNegMoreThan2 = zNeg(piW-1)
  val zNegShifted   = Mux(zNegMoreThan2 === 1.U, zNeg(piW-1, 3), zNeg(piW-2, 2))
  val zNegRounded   = zNegShifted + Mux(zNegMoreThan2 === 1.U, zNeg(2), zNeg(1))
  val zNegEx0  = exBias.U(exW.W) + zNegMoreThan2
  val zNegMan0 = zNegRounded(manW-1, 0)

  // -----------------------------------------------------------------------
  // check if z is special value

  val xNaNInf   = io.flags.isSpecial && io.flags.flag === ACosSpecialValue.xNaNInf
  val xZero     = io.flags.isSpecial && io.flags.flag === ACosSpecialValue.xZero
  val xLargePos = io.flags.isSpecial && io.flags.flag === ACosSpecialValue.xLargePos
  val xLargeNeg = io.flags.isSpecial && io.flags.flag === ACosSpecialValue.xLargeNeg

  // if x ==  0, return half pi
  // if x == +1, return 0
  // if x == -1, return pi
  // if x == nan/inf, return nan
  val zman = Mux(xNaNInf,   Cat(1.U(1.W), 0.U((manW-1).W)),
             Mux(xZero,     realHalfPi.man.toBigInt.U(manW.W),
             Mux(xLargePos, 0.U(manW.W),
             Mux(xLargeNeg, realPi.man.toBigInt.U(manW.W),
             Mux(io.flags.xsgn === 0.U, zMan0, zNegMan0)))))

  val zex  = Mux(xNaNInf,   Fill(exW, 1.U(1.W)),
             Mux(xZero,     realHalfPi.ex.U(exW.W),
             Mux(xLargePos, 0.U(exW.W),
             Mux(xLargeNeg, realPi.ex.U(exW.W),
             Mux(io.flags.xsgn === 0.U, zEx0, zNegEx0)))))

  val z = enableIf(io.en, Cat(0.U(1.W), zex, zman))
  io.z := ShiftRegister(z, nStage)
}
