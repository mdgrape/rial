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
  val selectNone       = 0.U(selectFuncW.W)
  val selectSqrt       = 1.U(selectFuncW.W)
  val selectInvSqrt    = 2.U(selectFuncW.W)
  val selectReciprocal = 3.U(selectFuncW.W)
  val selectSinPi      = 4.U(selectFuncW.W)
  val selectCosPi      = 5.U(selectFuncW.W)
  val selectACos       = 6.U(selectFuncW.W)
  val selectATan2      = 7.U(selectFuncW.W)
  val selectExp        = 8.U(selectFuncW.W)
  val selectLog        = 9.U(selectFuncW.W)
}

class MathFunctions(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  // TODO: stage config
  val nStage = stage.total
  def getStage() = nStage

  val polySpec = new PolynomialSpec(spec, nOrder, adrW, extraBits,
    enableRangeCheck, enablePolynomialRounding)

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order = polySpec.order

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

  val sqrtPre   = Module(new SqrtPreProcess (spec, polySpec, stage))
  val sqrtTab   = Module(new SqrtTableCoeff (spec, polySpec, maxAdrW, maxCbit, stage))
  val sqrtOther = Module(new SqrtOtherPath  (spec, polySpec, stage))
  val sqrtPost  = Module(new SqrtPostProcess(spec, polySpec, stage))

  sqrtPre.io.x   := io.x
  sqrtTab.io.adr := sqrtPre.io.adr
  sqrtOther.io.x := io.x

  val invsqrtTab   = Module(new InvSqrtTableCoeff (spec, polySpec, maxAdrW, maxCbit, stage))
  val invsqrtOther = Module(new InvSqrtOtherPath  (spec, polySpec, stage))
  val invsqrtPost  = Module(new InvSqrtPostProcess(spec, polySpec, stage))

  invsqrtTab.io.adr := sqrtPre.io.adr
  invsqrtOther.io.x := io.x

  val recPre   = Module(new ReciprocalPreProcess (spec, polySpec, stage))
  val recTab   = Module(new ReciprocalTableCoeff (spec, polySpec, maxAdrW, maxCbit, stage))
  val recOther = Module(new ReciprocalOtherPath  (spec, polySpec, stage))
  val recPost  = Module(new ReciprocalPostProcess(spec, polySpec, stage))

  recPre.io.x   := io.x
  recTab.io.adr := recPre.io.adr
  recOther.io.x := io.x

  val sinPiPre   = Module(new SinPiPreProcess (spec, polySpec, stage))
  val sinPiTab   = Module(new SinPiTableCoeff (spec, polySpec, maxAdrW, maxCbit, stage))
  val sinPiOther = Module(new SinPiOtherPath  (spec, polySpec, stage))
  val sinPiPost  = Module(new SinPiPostProcess(spec, polySpec, stage))

  val cosPiPre   = Module(new CosPiPreProcess (spec, polySpec, stage))
  val cosPiOther = Module(new CosPiOtherPath  (spec, polySpec, stage))

  sinPiPre.io.x            := io.x
  cosPiPre.io.x            := io.x

  sinPiTab.io.adr          := Mux(io.sel === SelectFunc.selectSinPi,
                                  sinPiPre.io.adr, cosPiPre.io.adr)

  sinPiOther.io.xConverted := sinPiPre.io.xConverted
  cosPiOther.io.xConverted := cosPiPre.io.xConverted

  sinPiOther.io.x          := io.x
  cosPiOther.io.x          := io.x

  //                      we are here
  //                               |
  // .-------. +=================+ v .-------------.
  // |       |=+  .------------. +===|(dx) poly    |   .---------.
  // |pre-   |----| tableCoeff |=====|(cs) interp  |---| post-   |
  // |process|-.  '------------'     '-------------' .-| process |
  // '-------' |  .--------------------------------. | '---------'
  //           '--| non-table path (e.g. taylor)   |-'
  //              '--------------------------------'

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, maxCbit, stage))

  val nullTab = 0.U.asTypeOf(new TableCoeffInput(maxCbit))

  if(order != 0) {
    polynomialEval.io.dx.get := MuxCase(0.U, Seq(
      (io.sel === SelectFunc.selectSqrt)       -> sqrtPre .io.dx.get,
      (io.sel === SelectFunc.selectInvSqrt)    -> sqrtPre .io.dx.get, // same as sqrt
      (io.sel === SelectFunc.selectReciprocal) -> recPre  .io.dx.get,
      (io.sel === SelectFunc.selectSinPi)      -> sinPiPre.io.dx.get,
      (io.sel === SelectFunc.selectCosPi)      -> cosPiPre.io.dx.get
    ))
  }
  polynomialEval.io.coeffs.cs <> MuxCase(nullTab.cs, Seq(
    (io.sel === SelectFunc.selectSqrt)       -> sqrtTab   .io.cs.cs,
    (io.sel === SelectFunc.selectInvSqrt)    -> invsqrtTab.io.cs.cs,
    (io.sel === SelectFunc.selectReciprocal) -> recTab    .io.cs.cs,
    (io.sel === SelectFunc.selectSinPi)      -> sinPiTab  .io.cs.cs,
    (io.sel === SelectFunc.selectCosPi)      -> sinPiTab  .io.cs.cs // same as sinPi
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

  invsqrtPost.io.zother <> invsqrtOther.io.zother
  invsqrtPost.io.zres   := polynomialEval.io.result

  recPost.io.zother <> recOther.io.zother
  recPost.io.zres   := polynomialEval.io.result

  sinPiPost.io.zother := Mux(io.sel === SelectFunc.selectSinPi,
    sinPiOther.io.zother, cosPiOther.io.zother)
  sinPiPost.io.zres   := polynomialEval.io.result

  val z0 = MuxCase(0.U, Seq(
    (io.sel === SelectFunc.selectSqrt)       -> sqrtPost.io.z,
    (io.sel === SelectFunc.selectInvSqrt)    -> invsqrtPost.io.z,
    (io.sel === SelectFunc.selectReciprocal) -> recPost.io.z,
    (io.sel === SelectFunc.selectSinPi)      -> sinPiPost.io.z,
    (io.sel === SelectFunc.selectCosPi)      -> sinPiPost.io.z // same as sinPi
  ))

  io.z := ShiftRegister(z0, stage.total)
}
