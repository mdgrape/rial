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
  val W          = 4
  val None       = 0.U(W.W)
  val Sqrt       = 1.U(W.W)
  val InvSqrt    = 2.U(W.W)
  val Reciprocal = 3.U(W.W)
  val SinPi      = 4.U(W.W)
  val CosPi      = 5.U(W.W)
  val ACos       = 6.U(W.W)
  val ATan2      = 7.U(W.W)
  val Exp        = 8.U(W.W)
  val Log        = 9.U(W.W)
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
    val sel = Input(UInt(SelectFunc.W.W))
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

  sinPiTab.io.adr          := Mux(io.sel === SelectFunc.SinPi,
                                  sinPiPre.io.adr, cosPiPre.io.adr)

  sinPiOther.io.xConverted := sinPiPre.io.xConverted
  cosPiOther.io.xConverted := cosPiPre.io.xConverted

  sinPiOther.io.x          := io.x
  cosPiOther.io.x          := io.x

  val acosPre   = Module(new ACosPreProcess (spec, polySpec, stage))
  val acosTab   = Module(new ACosTableCoeff (spec, polySpec, maxAdrW, maxCbit, stage))
  val acosOther = Module(new ACosOtherPath  (spec, polySpec, stage))
  val acosPost  = Module(new ACosPostProcess(spec, polySpec, stage))

  acosPre.io.x   := io.x
  acosTab.io.adr := acosPre.io.adr
  acosOther.io.x := io.x

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
      (io.sel === SelectFunc.Sqrt)       -> sqrtPre .io.dx.get,
      (io.sel === SelectFunc.InvSqrt)    -> sqrtPre .io.dx.get, // same as sqrt
      (io.sel === SelectFunc.Reciprocal) -> recPre  .io.dx.get,
      (io.sel === SelectFunc.SinPi)      -> sinPiPre.io.dx.get,
      (io.sel === SelectFunc.CosPi)      -> cosPiPre.io.dx.get,
      (io.sel === SelectFunc.ACos)       -> acosPre .io.dx.get
    ))
  }
  polynomialEval.io.coeffs.cs <> MuxCase(nullTab.cs, Seq(
    (io.sel === SelectFunc.Sqrt)       -> sqrtTab   .io.cs.cs,
    (io.sel === SelectFunc.InvSqrt)    -> invsqrtTab.io.cs.cs,
    (io.sel === SelectFunc.Reciprocal) -> recTab    .io.cs.cs,
    (io.sel === SelectFunc.SinPi)      -> sinPiTab  .io.cs.cs,
    (io.sel === SelectFunc.CosPi)      -> sinPiTab  .io.cs.cs,// same as sinPi
    (io.sel === SelectFunc.ACos)       -> acosTab   .io.cs.cs
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

  sinPiPost.io.zother := Mux(io.sel === SelectFunc.SinPi,
    sinPiOther.io.zother, cosPiOther.io.zother)
  sinPiPost.io.zres   := polynomialEval.io.result

  acosPost.io.zother <> acosOther.io.zother
  acosPost.io.zres   := polynomialEval.io.result

  val z0 = MuxCase(0.U, Seq(
    (io.sel === SelectFunc.Sqrt)       -> sqrtPost.io.z,
    (io.sel === SelectFunc.InvSqrt)    -> invsqrtPost.io.z,
    (io.sel === SelectFunc.Reciprocal) -> recPost.io.z,
    (io.sel === SelectFunc.SinPi)      -> sinPiPost.io.z,
    (io.sel === SelectFunc.CosPi)      -> sinPiPost.io.z, // same as sinPi
    (io.sel === SelectFunc.ACos)       -> acosPost.io.z
  ))

  io.z := ShiftRegister(z0, stage.total)
}
