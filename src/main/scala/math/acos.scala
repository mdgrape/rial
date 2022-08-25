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

import rial.math._

// Calculate acos(x) if x is in range [-1, 1].
// otherwise, 0 or pi for positive and negative x, respectively.
//
// Here, we calculate acos(x) by using the following two stages.
//
// 1. x -> sqrt(1 - |x|)
// 2. x -> acos(1 - x^2)
//
// acos(1 - sqrt(1-|x|)^2) = acos(1-(1-|x|)) = acos(|x|)
//
// So, we have additional preprocess stage that calculates 1-|x|.
//

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

object ACosSpecialValue {
  val W = 2
  val xZero    = 0.U(W.W)
  val xOne     = 1.U(W.W)
  val xNormal  = 2.U(W.W)
  val xNaNInf  = 3.U(W.W)
}

class ACosFlags extends Bundle {
  val xsgn    = UInt(1.W)
  val special = UInt(ACosSpecialValue.W.W)
}

class ACosPhase1PreProcess(
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
    val adr = Output(UInt((exAdrW+adrW).W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }

    // input to SqrtOtherPath
    val y   = new DecomposedRealOutput(spec)
    val special = Output(new ACosFlags())
  })

  val xsgn = enable(io.en, io.x.sgn)
  val xex  = enable(io.en, io.x.ex )
  val xman = enable(io.en, io.x.man)

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

  val xConvertedEx  = Mux(io.x.zero || xConvertedEx0 === exBias.U, exBias.U(exW.W), Mux(xLargerThan1, 0.U, xConvertedEx0))
  val xConvertedMan = Mux(io.x.zero || xConvertedEx0 === exBias.U, 0.U(manW.W),     Mux(xLargerThan1, 0.U, xConvertedMan0))

  // -------------------------------------------------------------------------
  // do the same thing as sqrt

  val adr = enable(io.en, Cat(xConvertedEx(0), xConvertedMan(manW-1, dxW)))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx = enable(io.en, Cat(~xConvertedMan(dxW-1), xConvertedMan(dxW-2, 0)))
    io.dx.get := ShiftRegister(dx, nStage)
  }

//   printf("cir: xConvertedEx  = %d\n", xConvertedEx )
//   printf("cir: xConvertedMan = %b\n", xConvertedMan)

  // pass 1-|x| to sqrtOtherPath
  io.y.sgn  := ShiftRegister(enable(io.en, 0.U(1.W)),       nStage)
  io.y.ex   := ShiftRegister(enable(io.en, xConvertedEx),  nStage)
  io.y.man  := ShiftRegister(enable(io.en, xConvertedMan), nStage)
  io.y.zero := ShiftRegister(enable(io.en, xLargerThan1),   nStage)
  io.y.inf  := ShiftRegister(enable(io.en, io.x.inf),       nStage)
  io.y.nan  := ShiftRegister(enable(io.en, io.x.nan),       nStage)

  // -------------------------------------------------------------------------
  // check if special value

  val xNaNInf = io.x.inf || io.x.nan
  val xZero   = io.x.zero || xConvertedEx0 === exBias.U
  val specialFlag = Mux(xNaNInf,      ACosSpecialValue.xNaNInf,
                    Mux(xZero,        ACosSpecialValue.xZero,
                    Mux(xLargerThan1, ACosSpecialValue.xOne,
                                      ACosSpecialValue.xNormal)))
  io.special.special := ShiftRegister(enable(io.en, specialFlag), nStage)
  io.special.xsgn := ShiftRegister(enable(io.en, xsgn), nStage)
}

// ============================================================================

class ACosPhase2PreProcess(
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
  })

  val xex  = enable(io.en, io.x.ex )
  val xman = enable(io.en, io.x.man)

  val xShift = exBias.U - xex
  val xAligned = Cat(1.U(1.W), xman) >> xShift

  val adr = enable(io.en, xAligned(manW-1, dxW))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx = enable(io.en, Cat(~xAligned(dxW-1), xAligned(dxW-2, 0)))
    io.dx.get := ShiftRegister(dx, nStage)
  }
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class ACosTableCoeff(
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

    io.cs.cs(0) := enable(io.en, coeff(0))

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
    io.cs := enable(io.en, coeffs)
  }
}
object ACosTableCoeff {
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

class ACosPostProcess(
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
    val x     = Flipped(new DecomposedRealOutput(spec))
    val flags = Input(new ACosFlags())
    val zres  = Input(UInt(fracW.W))
    val z     = Output(UInt(spec.W.W))
  })

  // -----------------------------------------------------------------------
  // calc acos(|x|)

  val lhs = Cat(1.U(1.W), io.zres)
  val rhs = Cat(1.U(1.W), io.x.man)

  val zProd = lhs * rhs
  val zMoreThan2 = zProd((fracW+1) + (manW+1)-1)
  val zShifted   = Mux(zMoreThan2, zProd((manW+1)+(fracW+1)-2, (fracW+1)  ),
                                   zProd((manW+1)+(fracW+1)-3, (fracW+1)-1))
  val zRounded   = zShifted +& Mux(zMoreThan2, zProd((fracW+1)-1), zProd((fracW+1)-2))
  val zMoreThan2AfterRound = zRounded(manW)
  val zExInc = zMoreThan2 | zMoreThan2AfterRound
  val zMan0  = zRounded(manW-1, 0)
  val zEx0   = io.x.ex + zExInc // table result has ex == exBias

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

  val special = io.flags.special
  val xNaNInf = special === ACosSpecialValue.xNaNInf
  val xZero   = special === ACosSpecialValue.xZero
  val xOne    = special === ACosSpecialValue.xOne
  // if x == 0, return half pi
  // if x == 1, return 0 or pi
  // if x == nan/inf, return nan
  val zman = Mux(xNaNInf, Cat(1.U(1.W), 0.U((manW-1).W)),
             Mux(xZero, realHalfPi.man.toBigInt.U(manW.W),
             Mux(xOne && io.flags.xsgn === 0.U, 0.U(manW.W),
             Mux(xOne && io.flags.xsgn === 1.U, realPi.man.toBigInt.U(manW.W),
             Mux(io.flags.xsgn === 0.U, zMan0, zNegMan0)))))

  val zex  = Mux(xNaNInf, Fill(exW, 1.U(1.W)),
             Mux(xZero, realHalfPi.ex.U(exW.W),
             Mux(xOne && io.flags.xsgn === 0.U, 0.U(exW.W),
             Mux(xOne && io.flags.xsgn === 1.U, realPi.ex.U(exW.W),
             Mux(io.flags.xsgn === 0.U, zEx0, zNegEx0)))))

  val z = enable(io.en, Cat(0.U(1.W), zex, zman))
  io.z := ShiftRegister(z, nStage)
}
