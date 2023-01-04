//% @file sigmoid.scala
//
// square root function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

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
import rial.arith.FloatChiselUtil

import rial.math._

// An implementation of the standard sigmoid function.

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class SigmoidPreProcessOutput(val spec: RealSpec) extends Bundle {
  val xsgn  = Output(UInt(1.W))
  val zone  = Output(Bool())
  val zzero = Output(Bool())
  val znan  = Output(Bool())
}

class SigmoidPreProcess(
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
    val out = new SigmoidPreProcessOutput(spec)
  })

  val rangeMaxLog2   = SigmoidSim.tableRangeMaxLog2(manW)
  val rangeMaxEx     = exBias + rangeMaxLog2 - 1
  val xExLargeEnough = rangeMaxEx.U < io.x.ex

  val znan  = io.x.nan
  val zone  = io.x.sgn === 0.U && xExLargeEnough
  val zzero = io.x.sgn === 1.U && xExLargeEnough

  val xScaleShift = rangeMaxEx.U(exW.W) - io.x.ex

  val xScaled = Cat(1.U, io.x.man) >> xScaleShift

  printf("sigmoid   : xScaled = %b\n", xScaled)

  // we use hidden bit, not only the mantissa
  val adr  = enable(io.en, xScaled(manW, dxW+1))
  io.adr := ShiftRegister(adr, nStage)
  printf("sigmoid   : adr = %b\n", adr)

  if(order != 0) {
    val dx   = enable(io.en, Cat(~xScaled(dxW), xScaled(dxW-1, 0)))
    io.dx.get := ShiftRegister(dx, nStage)
    printf("sigmoid   : dx = %b\n", dx)
  }

  io.out.xsgn  := ShiftRegister(io.x.sgn, nStage)
  io.out.znan  := ShiftRegister(znan , nStage)
  io.out.zone  := ShiftRegister(zone , nStage)
  io.out.zzero := ShiftRegister(zzero, nStage)
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class SigmoidTableCoeff(
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
    val adr = Input(UInt(adrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {

    val tableI = SigmoidSim.tableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    // sign mode 1 = 2's complement and no sign bit
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/1)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    assert(maxCbit(0) == fracW)
    assert(coeff(0).getWidth == fracW)

    io.cs.cs(0) := enable(io.en, coeff(0))

  } else {

    val tableI = SigmoidSim.tableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(0 <= diffWidth)

      if(diffWidth != 0) {
        val ci  = coeff(i)
//         val msb = ci(cbit(i)-1)
        val msb = ci.head(1)
        coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffs.cs(i) := coeff(i)
      }
      printf(f"sigmoid   : coeff(${i}) = %%d(%%b)\n", coeffs.cs(i).asSInt,  coeffs.cs(i))
    }
    io.cs := enable(io.en, coeffs)
  }
}
object SigmoidTableCoeff {
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
      return SigmoidSim.tableGeneration( order, adrW, spec.manW, fracW ).cbit
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
      return SigmoidSim.tableGeneration( order, adrW, spec.manW, fracW ).calcWidth
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

class SigmoidPostProcess(
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

  val calcW  = manW + extraBits

  val io = IO(new Bundle {
    val en     = Input(UInt(1.W))
    val preout = Flipped(new SigmoidPreProcessOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val z      = Output(UInt(spec.W.W))
  })

  printf("sigmoid   : zres = %b\n", io.zres)

  val zsgn  = 0.U(1.W)
  val znan  = io.preout.znan
  val zone  = io.preout.zone
  val zzero = io.preout.zzero

  val zPos = (BigInt(1) << (calcW+1)).U((calcW+2).W) - io.zres
  val zNeg = io.zres
  val zRes = Mux(io.preout.xsgn === 0.U, zPos, zNeg)

  printf("sigmoid   : zPos = %b : %b - %b\n", zPos, (BigInt(1) << (calcW+1)).U((calcW+2).W), io.zres)
  printf("sigmoid   : zNeg = %b\n", zNeg)
  printf("sigmoid   : zRes = %b\n", zRes)

  val zShift = Mux(zRes === 0.U, zRes.getWidth.U, PriorityEncoder(Reverse(zRes)))
  val zmanW1 = zRes << zShift
  val zmanW1Rounded = (zmanW1 >> (extraBits+1)) + zmanW1(extraBits)

  printf("sigmoid   : zShift = %d\n", zShift)
  printf("sigmoid   : zmanW1 = %b\n", zmanW1)
  printf("sigmoid   : zmanW1Rounded = %b\n", zmanW1Rounded)

  val zman = Mux(zmanW1Rounded > (BigInt(1) << (manW+1)).U, 0.U(manW.W), zmanW1Rounded(manW-1, 0))

  printf("sigmoid   : zman = %b\n", zman)

  val zex = exBias.U(exW.W) - zShift

  val z0 =
    Mux(znan,  Cat(zsgn, Fill(exW, 1.U(1.W)), 1.U(1.W), 0.U((manW-1).W)),
    Mux(zone,  Cat(zsgn, exBias.U, 0.U(manW.W)),
    Mux(zzero, 0.U(spec.W.W), Cat(zsgn, zex, zman))))

  val z = enable(io.en, z0)
  io.z := ShiftRegister(z, nStage)
}
