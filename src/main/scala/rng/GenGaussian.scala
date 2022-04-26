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

//
// sin/cos are used to determine the argument, so y can be a fixed point.
// By doing this, the resolution (the difference between the current and the
// next value) around 2piy ~ 0, pi, 2pi will be almost the same. It is not a
// problem.
//
class SinCos2Pi(
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
    val isSin = Input(Bool())
    val x     = Input(UInt(rndW.W))
    val z     = Output(UInt(spec.W.W))
  })

  // detect special values
  val xExactQuot = io.x(rndW-3, 0) === 0.U
  val xHighIs0   = io.x(rndW-1, rndW-2) === 0.U
  val xHighIs1   = io.x(rndW-1, rndW-2) === 1.U
  val xHighIs2   = io.x(rndW-1, rndW-2) === 2.U
  val xHighIs3   = io.x(rndW-1, rndW-2) === 3.U

  // convert x in [0, 1) to [0, 1/4).
  // Since we convert x range, in case of x is exactly N/4, the resulting x
  // will become indistinguishable from 0. but sin(0) and sin(2pi/4) completely
  // differ to each other. We need to check if x is exactly 1/4, 2/4, or 3/4.
  val zone = Wire(Bool())
  val zsgn = Wire(UInt(1.W))
  val x = Wire(UInt((rndW-2).W))
  when(io.isSin) {

    x := Mux(io.x(rndW-2) === 1.U, ~(io.x(rndW-3, 0))+1.U, io.x(rndW-3, 0))

    // x is exactly 1/4 or 3/4
    zone := (xHighIs1 || xHighIs3) && xExactQuot
    // if x is exactly 0 or 1/2, then the result is exactly zero, not -0.
    zsgn := io.x(rndW-1) & ~((xHighIs0 || xHighIs2) && xExactQuot).asUInt

  } otherwise {

    x := Mux(io.x(rndW-2) === 0.U, ~(io.x(rndW-3, 0))+1.U, io.x(rndW-3, 0))

    // x is exactly 0 or 1/2
    zone := (xHighIs0 || xHighIs2) && xExactQuot
    // if x is exactly 1/4 or 3/4, then the result is exactly zero, not -0.
    zsgn := (io.x(rndW-1) ^ io.x(rndW-2)) &
            ~((xHighIs1 || xHighIs3) && xExactQuot).asUInt
  }
  val adr = x(x.getWidth-1, x.getWidth-adrW)

  // in a range x in [0, 1/4), 4 < sin(2pix)/x <= 2pi.
  // So 1/2 < sin(2pix)/8x <= pi/4 = 0.78.. < 1.
  val tableD = new FuncTableDouble( x0 => {
    val x = x0 / 4.0 // convert [0, 1) to the input range, [0, 1/4)
    sin(x * Pi * 2.0) / (8.0*x)
  }, order )
  tableD.addRange(0.0, 1.0, 1<<adrW)
  val tableI = new FuncTableInt( tableD, fracW )
  val cbit = tableI.cbit

  val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
  val coeff = getSlices(coeffTable(adr), coeffWidth)

  val polyEval = Module(new PolynomialEval(spec, polySpec, cbit, stage))

  val coeffs = Wire(new TableCoeffInput(cbit))
  for (i <- 0 to order) {
    coeffs.cs(i) := coeff(i)
  }

  polyEval.io.coeffs := coeffs
  if(order != 0) {
    val dx = Cat(~x(x.getWidth-adrW-1), x(x.getWidth-adrW-2, x.getWidth-adrW-dxW))
    polyEval.io.dx.get := dx
  }

  val polyRes = polyEval.io.result // sin(2pix)/8x, in [0.5, 1). W = fracW

  // multiply x * sin(2pix)/x

  val xzero  = x === 0.U
  val xShift = Mux(xzero, 0.U, PriorityEncoder(Reverse(x)))
  val xAligned = (x << xShift)(x.getWidth-1, 0)
  assert(xAligned(x.getWidth-1) === 1.U || xzero)

  val zProd      = xAligned * polyRes
  val zProdW     = x.getWidth + fracW

  val zMoreThan2 = zProd(fracW+x.getWidth - 1)
  val zShifted   = Mux(zMoreThan2, zProd(zProdW-2, zProdW-manW-1),
                                   zProd(zProdW-3, zProdW-manW-2))

  val zRounded   = zShifted +& Mux(zMoreThan2, zProd(zProdW-manW-2),
                                               zProd(zProdW-manW-3))
  val zMoreThan2AfterRound = zRounded(manW)
  val zmanW1     = Mux(zMoreThan2AfterRound, Cat(1.U(1.W), 0.U(manW.W)),
                                             Cat(1.U(1.W), zRounded(manW-1, 0)))
  val zexInc     = zMoreThan2 + zMoreThan2AfterRound

  val zex  = Mux(zone, exBias.U(exW.W), Mux(xzero, 0.U, (exBias-1).U(exW.W) - xShift + zexInc))
  val zman = Mux(zone || xzero, 0.U, zmanW1(manW-1, 0))

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := z
}
