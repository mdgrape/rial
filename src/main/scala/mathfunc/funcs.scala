//% @file mathFunc.scala
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

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
import rial.mathfunc._

object SelectFunc {
  val selectFuncW = 4
  val selectSqrt  = 1.U(selectFuncW.W)
}

class MathFunctions(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val nStage = stage.total
  def getStage() = nStage

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order = if(adrW == manW) {0} else {nOrder}

  val maxAdrW  = adrW + 4        // TODO automatically calculate this from spec
  val maxCbit  = Seq(27, 21, 20) // TODO ditto
  val maxCalcW = Seq(27, 22, 20) // TODO ditto

  def getMaxAdrW()  = maxAdrW
  def getMaxCbit()  = maxCbit
  def getMaxCalcW() = maxCalcW

  val io = IO(new Bundle {
    val sel = Input(UInt(SelectFunc.selectFuncW.W))
    val x = Input (UInt(spec.W.W))
    val y = Input (UInt(spec.W.W))
    val z = Output(UInt(spec.W.W))
  })

  // .-------. .-----------------.   .-------------.
  // |       |-'  .------------. '---| polynomial  |   .---------.
  // |pre-   |====|(adr) table |-----| interpolate |---| post-   |
  // |process|=+  '------------'     '-------------' .-| process |
  // '-------' I  .--------------------------------. | '---------'
  //           +==|(x) non-table path (e.g. taylor)|-'
  //            ^ '--------------------------------'
  //            |
  //   we are here

  val sqrtPre   = Module(new SqrtPreProcess (spec, order, adrW, extraBits,                   stage, enableRangeCheck, enablePolynomialRounding))
  val sqrtTab   = Module(new SqrtTableCoeff (spec, order, adrW, extraBits, maxAdrW, maxCbit, stage, enableRangeCheck, enablePolynomialRounding))
  val sqrtOther = Module(new SqrtOtherPath  (spec, order, adrW, extraBits,                   stage, enableRangeCheck, enablePolynomialRounding))
  val sqrtPost  = Module(new SqrtPostProcess(spec, order, adrW, extraBits,                   stage, enableRangeCheck, enablePolynomialRounding))

  sqrtPre.io.x   := io.x
  sqrtTab.io.adr := sqrtPre.io.adr
  sqrtOther.io.x := io.x

  //                      we are here
  //                               |
  // .-------. +=================+ v .-------------.
  // |       |=+  .------------. +===|(dx) poly    |   .---------.
  // |pre-   |----| tableCoeff |=====|(cs) interp  |---| post-   |
  // |process|-.  '------------'     '-------------' .-| process |
  // '-------' |  .--------------------------------. | '---------'
  //           '--| non-table path (e.g. taylor)   |-'
  //              '--------------------------------'

  val polynomialEval = Module(new PolynomialEval(spec, order, adrW, extraBits, maxCbit, stage, enablePolynomialRounding))

  val nullTab = 0.U.asTypeOf(new TableCoeffInput(maxCbit))

  polynomialEval.io.dx     := MuxCase(0.U, Seq(
    (io.sel === SelectFunc.selectSqrt) -> sqrtPre.io.dx
  ))
  polynomialEval.io.coeffs.cs := MuxCase(nullTab.cs, Seq(
    (io.sel === SelectFunc.selectSqrt) -> sqrtTab.io.cs.cs
  ))

  //                                         we are here
  // .-------. .-----------------.   .-------------.  |
  // |       |-'  .------------. '---| polynomial  |  v .--------------.
  // |pre-   |----| tableCoeff |-----| interpolate |====|(zres)   post |
  // |process|-.  '------------'     '-------------'  +=|(zother) proc |
  // '-------' |  .--------------------------------.  I '--------------'
  //           '--| non-table path (e.g. taylor)   |==+
  //              '--------------------------------'

  sqrtPost.io.zother <> sqrtOther.io.zother
  sqrtPost.io.zres   := polynomialEval.io.result

  val z0 = MuxCase(0.U, Seq(
    (io.sel === SelectFunc.selectSqrt) -> sqrtPost.io.z
  ))

  io.z := ShiftRegister(z0, stage.total)
}
