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

class TableCoeffIO( val cbit: Seq[Int] ) extends Bundle {
  val cs = Input(MixedVec(cbit.map(w => UInt(w.W))))
}

class TableEval(
  val spec:           RealSpec,
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

  val order = cbit.length

  val io = IO(new Bundle {
    val coeffs = new TableCoeffIO(cbit)
    val dx     = UInt(dxW.W)
    val result = Output(UInt(fracW.W))
  })

  val res = UInt(fracW.W)

  if(order == 0) {

    res := coeffs(0)

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

    res := io.coeffs.init.foldRight(io.coeffs.last)(
      (c,z) => hornerC( c.asSInt, z, io.dx.asSInt )
    )
  }
  io.result := ShiftRegister(res, nStage)
}

