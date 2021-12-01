//% @file sqrt.scala
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

// sqrt(x): floating => floating
// - if x < 0, returns 0.
class SqrtPreProcess(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val manW = spec.manW

  val io = IO(new Bundle {
    val x   = Input (UInt(spec.W.W))
    val adr = Output(UInt((1+adrW).W))
    val dx  = Output(UInt((manW-adrW).W))
  })

  io.adr := io.x(manW, manW-adrW) // include LSB of x.ex
  io.dx  := Cat(~io.x(manW-adrW-1), io.x(manW-adrW-2,0))
}

class SqrtTableCoeff(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val manW  = spec.manW
  val calcW = manW + extraBits

  val tableI = SqrtSim.sqrtTableGeneration( nOrder, adrW, manW, calcW )
  val cbits  = tableI.cbit

  // TODO: in case of order == 0, the width and table generation strategy are
  //       different from order == 2.

  val io = IO(new Bundle {
    val adr = Input (UInt((1+adrW).W))
    val c0  = Output(UInt(cbits(0).W)) // TODO: order == 0 ?
    val c1  = Output(UInt(cbits(1).W)) // TODO
    val c2  = Output(UInt(cbits(2).W)) // TODO
  })

  val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
  val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

  io.c0 := coeff(0)
  io.c1 := coeff(1)
  io.c2 := coeff(2)
}

// No pathway other than table interpolation. just calculate ex and sgn.
class SqrtOtherPath(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val x     = Input (UInt(spec.W.W))
    val zex   = Output(UInt(spec.exW.W))
    val zsgn  = Output(UInt(1.W))
    val znan  = Output(Bool())
    val zinf  = Output(Bool())
    val zzero = Output(Bool())
  })

  val (xsgn,  xex,  xman) = FloatChiselUtil.decompose(spec, io.x)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(spec, io.x)

  val xneg = if(spec.disableSign) {false.B} else {xsgn === 1.U(1.W)}

  val znan  = xnan
  val zinf  = xinf
  val zzero = xzero || xneg

  io.znan  := znan
  io.zinf  := zinf
  io.zzero := zzero

  val xExNobias = xex - exBias.U
  val zex0 = xExNobias(exW-1) ## xExNobias(exW-1, 1) // arith right shift
  io.zex  := Mux(zinf || znan, maskU(exW),
             Mux(zzero,        0.U(exW.W),
                               zex0 + exBias.U))
  io.zsgn := 0.U(1.W) // always positive.
}

class SqrtPostProcess(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val nStage = stage.total
  def getStage() = nStage

  val io = IO(new Bundle {
    // ex and some conditions
    val zex   = Input(UInt(spec.exW.W))
    val zsgn  = Input(UInt(1.W))
    val znan  = Input(Bool())
    val zinf  = Input(Bool())
    val zzero = Input(Bool())
    // table interpolation results
    val zres  = Input(UInt((spec.manW + extraBits).W))
    // output
    val z     = Output(UInt(spec.W.W))
  })

  val zman0 = dropLSB(extraBits, io.zres) +& io.zres(extraBits-1)
  val polynomialOvf = zman0(manW)
  val zeroFlush     = io.zinf || io.zzero || io.znan
  val zman          = Mux(zeroFlush,     Cat(io.znan, 0.U((manW-1).W)),
                      Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman0(manW-1,0)))

  val z0 = Cat(io.zsgn, io.zex, zman)

  io.z   := ShiftRegister(z0, nStage)
}
