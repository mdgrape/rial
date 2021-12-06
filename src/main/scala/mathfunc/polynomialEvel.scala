//% @file tableEval.scala
//
// polynomial approx using coefficient table
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._

import rial.arith.RealSpec
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._

class TableCoeffInput( val cbit: Seq[Int] ) extends Bundle {
  val cs = Input(MixedVec(cbit.map{w => UInt(w.W)}))
}

class PolynomialEval(
  val spec:           RealSpec,
  val nOrder:         Int,
  val adrW:           Int,
  val extraBits:      Int,
  val cbit:           Seq[Int],
  val stage:          PipelineStageConfig,
  val enableRounding: Boolean = false
) extends Module {

  val nStage = stage.total

  val manW  = spec.manW
  val dxW   = manW - adrW
  val fracW = manW + extraBits

  val order = if(adrW == manW) {0} else {nOrder}

  val io = IO(new Bundle {
    val coeffs = new TableCoeffInput(cbit)
    val dx     = Input(UInt(dxW.W))
    val result = Output(UInt(fracW.W))
  })

  val res = Wire(UInt(fracW.W))

  if(order == 0) {

    res := io.coeffs.cs(0)

  } else {

    def hornerC( c: SInt, z: SInt, dx: SInt ) : SInt = {
      val zw   = z.getWidth
      val dxBp = dx.getWidth-1
      val dxw  = min(zw, dxBp)
      val dxl  = dx(dxBp, dxBp-dxw).asSInt
      val prod = dxl * z
      val prod_sft = dropLSB(dxw, prod)
      val sum = c +& prod_sft // Extend for safe
      sum
    }

    val coeffS = io.coeffs.cs.map( x => x.asSInt )

    val resS = coeffS.init.foldRight(coeffS.last)(
      (c,z) => hornerC( c, z, io.dx.asSInt )
    )
    res := resS.asUInt
  }
  io.result := ShiftRegister(res, nStage)
}

