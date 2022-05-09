//% @file GenGaussian.scala
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.rng

import scala.language.reflectiveCalls
import scala.math._

import chisel3._
import chisel3.util._

import rial.arith._
import rial.math._
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._

//
// generate Gaussian using Box-Muller.
//
// z1 = sqrt(-2log(x)) cos(2pi*y)
// z2 = sqrt(-2log(x)) sin(2pi*y)
//

// -----------------------------------------------------------------------------
// Fixed -> FP sin/cos(2pi*y)
//
// sin/cos are used to determine the argument (not an absolute value), so y can
// be a fixed point.
// By doing this, the resolution (the difference between the current and the
// next value) around 2piy ~ 0, pi, 2pi will be almost the same. It is not a
// problem.
//

class BoxMullerSinCos2PiPreProcessOutput(val rndW: Int, val spec: RealSpec) extends Bundle {
  val zone     = Output(Bool())
  val zsgn     = Output(UInt(1.W))
  val zzero    = Output(Bool())
  val xShift   = Output(UInt(log2Up(rndW-2).W))
  val xAligned = Output(UInt((rndW-2).W))
}

object BoxMullerSinCos2PiTableCoeff {
  def genTable(
    polySpec: PolynomialSpec
  ): FuncTableInt = {

    // --------------------------------------------------------------------------
    // in a range x in [0, 1/4), 4 < sin(2pix)/x <= 2pi.
    // So 1/2 < sin(2pix)/8x <= pi/4 = 0.78.. < 1.
    val order  = polySpec.order
    val fracW  = polySpec.fracW
    val adrW   = polySpec.adrW

    val tableD = new FuncTableDouble( x0 => {
      val x = x0 / 4.0 // convert [0, 1) to the input range, [0, 1/4)
      sin(x * Pi * 2.0) / (8.0*x)
    }, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  def getCBits(
    polySpec: PolynomialSpec
  ): Seq[Int] = {
    val tableI = genTable(polySpec)
    tableI.cbit
  }
}

class BoxMullerSinCos2PiPreProc(
    rndW: Int,      // width of input fixedpoint
    spec: RealSpec, // output width
    polySpec: PolynomialSpec,
    stage: PipelineStageConfig
  ) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order  = polySpec.order
  val fracW  = polySpec.fracW
  val adrW   = polySpec.adrW
  val dxW    = polySpec.dxW

  assert(rndW > 2 + adrW + dxW) // manW + 2; satisfied by (i32, f32) or (i64, f64)

  val tableI = BoxMullerSinCos2PiTableCoeff.genTable(polySpec)
  val cbit = tableI.cbit

  // --------------------------------------------------------------------------

  val io = IO(new Bundle {
    val isSin = Input(Bool())
    val rnd   = Input(UInt(rndW.W))

    val cs    = Flipped(new TableCoeffInput(cbit))
    val dx    = if (order != 0) { Some(Output(UInt(dxW.W))) } else { None }
    val out   = new BoxMullerSinCos2PiPreProcessOutput(rndW, spec)
  })

  // detect special values
  val xExactQuot = io.rnd(rndW-3, 0) === 0.U
  val xHighIs0   = io.rnd(rndW-1, rndW-2) === 0.U
  val xHighIs1   = io.rnd(rndW-1, rndW-2) === 1.U
  val xHighIs2   = io.rnd(rndW-1, rndW-2) === 2.U
  val xHighIs3   = io.rnd(rndW-1, rndW-2) === 3.U

  // convert x in [0, 1) to [0, 1/4).
  // Since we convert x range, in case of x is exactly N/4, the resulting x
  // will become indistinguishable from 0. but sin(0) and sin(2pi/4) completely
  // differ to each other. We need to check if x is exactly 1/4, 2/4, or 3/4.
  val zone = Wire(Bool())
  val zsgn = Wire(UInt(1.W))
  val x = Wire(UInt((rndW-2).W))
  when(io.isSin) {

    x := Mux(io.rnd(rndW-2) === 1.U, ~(io.rnd(rndW-3, 0))+1.U, io.rnd(rndW-3, 0))

    // x is exactly 1/4 or 3/4
    zone := (xHighIs1 || xHighIs3) && xExactQuot
    // if x is exactly 0 or 1/2, then the result is exactly zero, not -0.
    zsgn := io.rnd(rndW-1) & ~((xHighIs0 || xHighIs2) && xExactQuot).asUInt

  } otherwise {

    x := Mux(io.rnd(rndW-2) === 0.U, ~(io.rnd(rndW-3, 0))+1.U, io.rnd(rndW-3, 0))

    // x is exactly 0 or 1/2
    zone := (xHighIs0 || xHighIs2) && xExactQuot
    // if x is exactly 1/4 or 3/4, then the result is exactly zero, not -0.
    zsgn := (io.rnd(rndW-1) ^ io.rnd(rndW-2)) &
           ~((xHighIs1 || xHighIs3) && xExactQuot).asUInt
  }

  if(order != 0) {
    val dx = Cat(~x(x.getWidth-adrW-1), x(x.getWidth-adrW-2, x.getWidth-adrW-dxW))
    io.dx.get := ShiftRegister(dx, nStage)
  }
  io.out.zsgn := ShiftRegister(zsgn, nStage)
  io.out.zone := ShiftRegister(zone, nStage)

  // --------------------------------------------------------------------------

  val adr = x(x.getWidth-1, x.getWidth-adrW)

  val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
  val coeff = getSlices(coeffTable(adr), coeffWidth)

  val coeffs = Wire(new TableCoeffInput(cbit))
  for (i <- 0 to order) {
    coeffs.cs(i) := coeff(i)
  }

  io.cs := ShiftRegister(coeffs, nStage)

  // --------------------------------------------------------------------------

  val zzero  = x === 0.U
  val xShift = Mux(zzero, 0.U, PriorityEncoder(Reverse(x)))
  val xAligned = (x << xShift)(x.getWidth-1, 0)
  assert(xAligned(x.getWidth-1) === 1.U || zzero)

  io.out.zzero    := ShiftRegister(zzero,    nStage)
  io.out.xShift   := ShiftRegister(xShift,   nStage)
  io.out.xAligned := ShiftRegister(xAligned, nStage)
}

class BoxMullerSinCos2PiPostProc(
    rndW: Int,      // width of input fixedpoint
    spec: RealSpec, // output width
    polySpec: PolynomialSpec,
    stage: PipelineStageConfig
  ) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order  = polySpec.order
  val fracW  = polySpec.fracW
  val adrW   = polySpec.adrW
  val dxW    = polySpec.dxW

  assert(rndW > 2 + adrW + dxW) // manW + 2; satisfied by (i32, f32) or (i64, f64)

  val io = IO(new Bundle {
    val pre  = Flipped(new BoxMullerSinCos2PiPreProcessOutput(rndW, spec))
    val zres = Input(UInt(fracW.W))
    val z    = Output(UInt(spec.W.W))
  })

  val zone     = io.pre.zone
  val zsgn     = io.pre.zsgn
  val zzero    = io.pre.zzero
  val xShift   = io.pre.xShift
  val xAligned = io.pre.xAligned

  val zProd      = xAligned * io.zres
  val zProdW     = xAligned.getWidth + fracW

  val zMoreThan2 = zProd(zProdW - 1)
  val zShifted   = Mux(zMoreThan2, zProd(zProdW-2, zProdW-manW-1),
                                   zProd(zProdW-3, zProdW-manW-2))

  val zRounded   = zShifted +& Mux(zMoreThan2, zProd(zProdW-manW-2),
                                               zProd(zProdW-manW-3))
  val zMoreThan2AfterRound = zRounded(manW)
  val zmanW1     = Mux(zMoreThan2AfterRound, Cat(1.U(1.W), 0.U(manW.W)),
                                             Cat(1.U(1.W), zRounded(manW-1, 0)))
  val zexInc     = zMoreThan2 + zMoreThan2AfterRound

  val zex  = Mux(zone, exBias.U(exW.W), Mux(zzero, 0.U, (exBias-1).U(exW.W) - xShift + zexInc))
  val zman = Mux(zone || zzero, 0.U, zmanW1(manW-1, 0))

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := ShiftRegister(z, nStage)
}


//
// z = sqrt(-2log(1-x)) at fixedpoint precision.
//
class Sqrt2LogX(
  rndW: Int,      // width of input fixedpoint
  spec: RealSpec, // output width
  polySpec: PolynomialSpec,
  stage: PipelineStageConfig
  ) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order  = polySpec.order
  val fracW  = polySpec.fracW
  val extraBits = fracW - manW
  val adrW   = polySpec.adrW
  val dxW    = polySpec.dxW

  assert(rndW >= adrW + dxW)

  val io = IO(new Bundle {
    val x = Input(UInt(rndW.W))
    val z = Output(UInt(spec.W.W))
  })

  val x = io.x

  // if rndW=32, xmin = 2^-32. then sqrt(-2 * log(2^-32)) ~ 6.66. intW = 3.
  val xmin = pow(2.0, -rndW)
  val zmax = sqrt(-2.0 * log(xmin))
  val intW = log2UpD(zmax)
  if(extraBits < intW) {
    println(f"WARNING: extraBits(${extraBits}) in polynomial is less than "+
            f"the recommended intW(${intW}).")
  }

  // sqrt(-2log(1-x)); to avoid inf, the x axis is reversed
  val tableD = new FuncTableDouble( x0 => {
    val scale = pow(2.0, -intW) // scale to make z < 1
    if(x0 < xmin) {
      zmax * scale
    } else {
      sqrt(-2*log(1.0-x0)) * scale
    }
  }, order )
  tableD.addRange(0.0, 1.0, 1<<adrW)
  val tableI = new FuncTableInt( tableD, fracW )
  val cbit = tableI.cbit
  val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)

  val polyEval = Module(new PolynomialEval(spec, polySpec, cbit, stage))

  val adr   = x(rndW-1, rndW-adrW)
  val coeff = getSlices(coeffTable(adr), coeffWidth)

  val coeffs = Wire(new TableCoeffInput(cbit))
  for (i <- 0 to order) {
    coeffs.cs(i) := coeff(i)
  }

  polyEval.io.coeffs := coeffs
  if(order != 0) {
    val dx = Cat(~x(x.getWidth-adrW-1), x(x.getWidth-adrW-2, x.getWidth-adrW-dxW))
    polyEval.io.dx.get := dx
  }

  val zres = polyEval.io.result // sin(2pix)/8x, in [0.5, 1). W = fracW

  // shift to normalize z

  val zShift = Mux(zres === 0.U, 0.U, PriorityEncoder(Reverse(zres)))
  val zAligned = (zres << zShift)(fracW-1, 0)
  assert(zAligned(fracW-1) === 1.U || zres === 0.U)

  assert(0 <= fracW-manW-2)
  val zman0 = zAligned(fracW-2, fracW-manW-1) +& zAligned(fracW-manW-2)
  val zmanMoreThan2AfterRound = zman0(manW)
  val zex0 = (exBias + intW - 1).U(exW.W) - zShift + zmanMoreThan2AfterRound

  val zman = zman0(manW-1, 0)
  val zex  = Mux(zres === 0.U, 0.U, zex0)
  val zsgn = 0.U(1.W)

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := z
}
