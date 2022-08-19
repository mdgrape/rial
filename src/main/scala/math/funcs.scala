//% @file mathFunc.scala
//
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

object SelectFunc {
  val W           =  4
  val None        =  0.U(W.W)
  val Sqrt        =  1.U(W.W)
  val InvSqrt     =  2.U(W.W)
  val Reciprocal  =  3.U(W.W)
  val Sin         =  4.U(W.W)
  val Cos         =  5.U(W.W)
  val ACos1       =  6.U(W.W)
  val ACos2       =  7.U(W.W)
  val ATan2Stage1 =  8.U(W.W)
  val ATan2Stage2 =  9.U(W.W)
  val Pow2        = 10.U(W.W)
  val Exp         = 11.U(W.W)
  val Log2        = 12.U(W.W)
  val Log         = 13.U(W.W)
}

class DecomposedRealOutput(val spec: RealSpec) extends Bundle {
  val sgn  = Output(UInt(1.W))
  val ex   = Output(UInt(spec.exW.W))
  val man  = Output(UInt(spec.manW.W))
  val zero = Output(Bool())
  val inf  = Output(Bool())
  val nan  = Output(Bool())
}
class DecomposeReal(val spec: RealSpec) extends Module {
  val io = IO(new Bundle {
    val real   = Input(UInt(spec.W.W))
    val decomp = new DecomposedRealOutput(spec)
  })
  val (sgn,  ex,  man) = FloatChiselUtil.decompose (spec, io.real)
  val (zero, inf, nan) = FloatChiselUtil.checkValue(spec, io.real)
  io.decomp.sgn  := sgn
  io.decomp.ex   := ex
  io.decomp.man  := man
  io.decomp.zero := zero
  io.decomp.inf  := inf
  io.decomp.nan  := nan
}

// # Overview
//
//             .--preStage    .--tableCalcGap .-- calcPostGap
//             |      .--preCalcGap   .--calcStage    .--postStage
//         ____|____  |       |  _____|_____  |  _____|_____
//        '         ' '       ' '           ' ' '           '
//       .----.-.----.-.-----.-.-----.-.-----.-.-----.-.-----.
//       |    |v|    |v|table|v|Calc1|v|Calc2|v|     |v|     |
// in -> |Pre1| |Pre2| :-----'-'-----'-'-----: |Post1| |Post2| -> out
//       |    | |    | |      non-table      | |     | |     |
//       '----'-'----'-'---------------------'-'-----'-'-----'
//
// TODO: consider setting nStage for each function. (like, sqrt does not need
//       multiple cycles in its preprocess, but sincos may need.)
//
class MathFuncPipelineConfig(
  val preStage:     PipelineStageConfig, // clock cycles in preprocess
  val calcStage:    PipelineStageConfig, // clock cycles in table/polynomial and non-table path
  val postStage:    PipelineStageConfig, // clock cycles in postprocess
  val preCalcGap:   Boolean, // if true, add register between preprocess and calculation stage
  val tableCalcGap: Boolean, // if true, add register between table and calculation stage (+1 to calcStage for OtherPath)
  val calcPostGap:  Boolean, // if true, add register between calculation and postprocess stage
  ) {

  def total = {
    preStage.total +
    calcStage.total +
    postStage.total +
    (if(preCalcGap)   {1} else {0}) +
    (if(tableCalcGap) {1} else {0}) +
    (if(calcPostGap)  {1} else {0})
  }
  def getString = {
    f"pre: ${preStage.getString}, " +
    f"calc: ${calcStage.getString}, " +
    f"post: ${postStage.getString}, " +
    f"pre-calc: ${preCalcGap}, table-calc: ${tableCalcGap}, calc-post: ${calcPostGap}"
  }
}

object MathFuncPipelineConfig {
  def none = {
    new MathFuncPipelineConfig(
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      false,
      false,
      false)
  }
}

// # Overview
//
//                                     table/Calc
//                                          |
//                .------.  _  .---------.  _  .------------.
// sel -----------|      |-|v|-'.-------.'-|v|-| chebyshev  |    _  .-------.
//    .---------. | pre- |-| |--| table |--| |-| polynomial |---|v|-| post- |-> z
// x -|decompose|-| proc | | |  '-------'  '-' '------------' .-| |-| proc  |
// y -|         |-|      |-| |-. .--------------------------. | '_' '-------'
//    '---------' '------' '-' '-|  non-table (taylor etc)  |-'  .
//                          .    '--------------------------'    |
//                          |                                    |
//                   Preproc/Calc                          Calc/Postproc
//

class MathFunctions(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : MathFuncPipelineConfig,
  val dxW0 : Option[Int] = None,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  assert(!spec.disableSign) // assuming the float type allows negative values

  // later we use this in ShiftRegister
  val pcGap = if(stage.preCalcGap ) {1} else {0}
  val tcGap = if(stage.tableCalcGap){1} else {0}
  val cpGap = if(stage.calcPostGap) {1} else {0}

  val nPreStage  = stage.preStage.total
  val nCalcStage = stage.calcStage.total
  val nPostStage = stage.postStage.total

  val nOtherStage = nCalcStage + tcGap
  val otherStage = PipelineStageConfig.atOut(nOtherStage)

  println(f"nPreStage  = ${nPreStage }")
  println(f"nCalcStage = ${nCalcStage}")
  println(f"nPostStage = ${nPostStage}")

  val nStage = stage.total
  def getStage = nStage

  val polySpec = new PolynomialSpec(spec, nOrder, adrW, extraBits, dxW0,
    enableRangeCheck, enablePolynomialRounding)
  val order = polySpec.order

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val maxCbit  = Seq(
    ACosTableCoeff.getCBits(spec, polySpec),
    SqrtTableCoeff.getCBits(spec, polySpec),
    InvSqrtTableCoeff.getCBits(spec, polySpec),
    ReciprocalTableCoeff.getCBits(spec, polySpec),
    SinCosTableCoeff.getCBits(spec, polySpec),
    ATan2Stage2TableCoeff.getCBits(spec, polySpec),
    ExpTableCoeff.getCBits(spec, polySpec),
    LogTableCoeff.getCBits(spec, polySpec)
    ).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  val maxCalcW = Seq(
    ACosTableCoeff.getCalcW(spec, polySpec),
    SqrtTableCoeff.getCalcW(spec, polySpec),
    InvSqrtTableCoeff.getCalcW(spec, polySpec),
    ReciprocalTableCoeff.getCalcW(spec, polySpec),
    SinCosTableCoeff.getCalcW(spec, polySpec),
    ATan2Stage2TableCoeff.getCalcW(spec, polySpec),
    ExpTableCoeff.getCalcW(spec, polySpec),
    LogTableCoeff.getCalcW(spec, polySpec)
    ).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  def getMaxCbit  = maxCbit
  def getMaxCalcW = maxCalcW

  println(f"acos        cbits = ${ACosTableCoeff       .getCBits(spec, polySpec)} calcW = ${ACosTableCoeff       .getCalcW(spec, polySpec)}")
  println(f"sqrt        cbits = ${SqrtTableCoeff       .getCBits(spec, polySpec)} calcW = ${SqrtTableCoeff       .getCalcW(spec, polySpec)}")
  println(f"invsqrt     cbits = ${InvSqrtTableCoeff    .getCBits(spec, polySpec)} calcW = ${InvSqrtTableCoeff    .getCalcW(spec, polySpec)}")
  println(f"rec         cbits = ${ReciprocalTableCoeff .getCBits(spec, polySpec)} calcW = ${ReciprocalTableCoeff .getCalcW(spec, polySpec)}")
  println(f"sincos      cbits = ${SinCosTableCoeff     .getCBits(spec, polySpec)} calcW = ${SinCosTableCoeff     .getCalcW(spec, polySpec)}")
  println(f"atan2Stage2 cbits = ${ATan2Stage2TableCoeff.getCBits(spec, polySpec)} calcW = ${ATan2Stage2TableCoeff.getCalcW(spec, polySpec)}")
  println(f"pow2        cbits = ${ExpTableCoeff        .getCBits(spec, polySpec)} calcW = ${ExpTableCoeff        .getCalcW(spec, polySpec)}")
  println(f"log2        cbits = ${LogTableCoeff        .getCBits(spec, polySpec)} calcW = ${LogTableCoeff        .getCalcW(spec, polySpec)}")
  println(f"maximum     cbits = ${maxCbit} calcW = ${maxCalcW}")

  val io = IO(new Bundle {
    val sel = Input(UInt(SelectFunc.W.W))
    val x = Input (UInt(spec.W.W))
    val y = Input (UInt(spec.W.W))
    val z = Output(UInt(spec.W.W))
  })

  val xdecomp = Module(new DecomposeReal(spec))
  val ydecomp = Module(new DecomposeReal(spec))
  xdecomp.io.real := io.x
  ydecomp.io.real := io.y

  val yIsLarger = io.x(spec.W-2, 0) < io.y(spec.W-2, 0) // cmp without sign bit

  // These registers delay sel, xdec, ydec by nPreStage and nCalcStage cycles.
  // The values of PC/CP registers can be used as if those are the output of
  // Preprocess/Calculation stage, respectively.
  // If nStages are zero, then it becomes a wire.
  val selPCReg  = ShiftRegister(io.sel,            nPreStage)
  val xdecPCReg = ShiftRegister(xdecomp.io.decomp, nPreStage)
  val ydecPCReg = ShiftRegister(ydecomp.io.decomp, nPreStage)
  val selPCGapReg  = ShiftRegister(selPCReg,  pcGap)
  val xdecPCGapReg = ShiftRegister(xdecPCReg, pcGap)
  val ydecPCGapReg = ShiftRegister(ydecPCReg, pcGap)

  val selCPReg  = ShiftRegister(selPCGapReg,  tcGap + nCalcStage)
  val xdecCPReg = ShiftRegister(xdecPCGapReg, tcGap + nCalcStage)
  val ydecCPReg = ShiftRegister(ydecPCGapReg, tcGap + nCalcStage)
  val selCPGapReg  = ShiftRegister(selCPReg,  cpGap)
  val xdecCPGapReg = ShiftRegister(xdecCPReg, cpGap)
  val ydecCPGapReg = ShiftRegister(ydecCPReg, cpGap)

  val yIsLargerPCReg    = ShiftRegister(yIsLarger, nPreStage)
  val yIsLargerPCGapReg = ShiftRegister(yIsLargerPCReg, pcGap)
  val yIsLargerCPReg    = ShiftRegister(yIsLargerPCGapReg, tcGap + nCalcStage)
  val yIsLargerCPGapReg = ShiftRegister(yIsLargerCPReg, cpGap)

  // ==========================================================================
  //
  // .-------.  _  .----------------.   .-------------.
  // |       |-|v|-' .------------. '---| polynomial  |    _  .---------.
  // |pre-   |=| |===|(adr) table |-----| interpolate |---|v|-| post-   |
  // |process|=| |=+ '------------'     '-------------' .-| |-| process |
  // '-------' '-' I .--------------------------------. | '_' '---------'
  //               +=|(x) non-table path (e.g. taylor)|-'
  //               ^ '--------------------------------'
  //               |
  // now we are here

  // --------------------------------------------------------------------------
  // acos

  val acosFlagReg = Reg(new ACosFlags())

  val acos1Pre  = Module(new ACosStage1PreProcess(spec, polySpec, stage.preStage))
  val acos2Pre  = Module(new ACosStage2PreProcess(spec, polySpec, stage.preStage))
  val acosTab   = Module(new ACosTableCoeff (spec, polySpec, maxCbit))
  val acosPost  = Module(new ACosPostProcess(spec, polySpec, stage.postStage))

  val acos1PreAdrPCGapReg = ShiftRegister(acos1Pre.io.adr, pcGap)

  acos1Pre.io.en        := (io.sel === SelectFunc.ACos1)
  acos1Pre.io.x         := xdecomp.io.decomp
  acos2Pre.io.en        := (io.sel === SelectFunc.ACos2)
  acos2Pre.io.x         := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  acosTab.io.en         := selPCGapReg === SelectFunc.ACos2
  acosTab.io.adr        := ShiftRegister(acos2Pre.io.adr, pcGap)

  when(selPCReg === SelectFunc.ACos1) {
    acosFlagReg := acos1Pre.io.special
  }

  // after preprocess
  // XXX Since `PCReg`s are delayed by nPreStage, the timing is the same as acosPre output.
  when(selPCReg =/= SelectFunc.ACos1) {
    assert(acos1Pre.io.adr === 0.U)
    assert(acos1Pre.io.dx.getOrElse(0.U) === 0.U)
  }
  when(selPCReg =/= SelectFunc.ACos2) {
    assert(acos2Pre.io.adr === 0.U)
    assert(acos2Pre.io.dx.getOrElse(0.U) === 0.U)
  }

  when(selPCGapReg =/= SelectFunc.ACos2) {
    assert(acosTab.io.cs.asUInt === 0.U)
  }

  // --------------------------------------------------------------------------
  // sqrt

  val sqrtPre   = Module(new SqrtPreProcess (spec, polySpec, stage.preStage))
  val sqrtTab   = Module(new SqrtTableCoeff (spec, polySpec, maxCbit))
  val sqrtOther = Module(new SqrtOtherPath  (spec, polySpec, otherStage))
  val sqrtPost  = Module(new SqrtPostProcess(spec, polySpec, stage.postStage))

  val sqrtPreAdrPCGapReg = ShiftRegister(sqrtPre.io.adr, pcGap)

  sqrtPre.io.en  := (io.sel === SelectFunc.Sqrt || io.sel === SelectFunc.InvSqrt)
  sqrtPre.io.x   := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  sqrtTab.io.en  := (selPCGapReg === SelectFunc.Sqrt || selPCGapReg === SelectFunc.ACos1)
  sqrtTab.io.adr := sqrtPreAdrPCGapReg | acos1PreAdrPCGapReg
  sqrtOther.io.x := Mux(selPCGapReg === SelectFunc.Sqrt, xdecPCGapReg, ShiftRegister(acos1Pre.io.y, pcGap))

  // after preprocess
  when(selPCReg =/= SelectFunc.Sqrt && selPCReg =/= SelectFunc.InvSqrt) {
    assert(sqrtPre.io.adr === 0.U)
    assert(sqrtPre.io.dx.getOrElse(0.U) === 0.U)
  }
  when(selPCGapReg =/= SelectFunc.Sqrt) {
    assert(sqrtTab.io.cs.asUInt === 0.U || selPCGapReg === SelectFunc.ACos1)
  }

  // --------------------------------------------------------------------------
  // invsqrt

  val invsqrtTab   = Module(new InvSqrtTableCoeff (spec, polySpec, maxCbit))
  val invsqrtOther = Module(new InvSqrtOtherPath  (spec, polySpec, otherStage))
  val invsqrtPost  = Module(new InvSqrtPostProcess(spec, polySpec, stage.postStage))

  // (preprocess is same as sqrt)
  // ------ Preprocess-Calculate ------
  invsqrtTab.io.en  := selPCGapReg === SelectFunc.InvSqrt
  invsqrtTab.io.adr := sqrtPreAdrPCGapReg
  invsqrtOther.io.x := xdecPCGapReg

  when(selPCGapReg =/= SelectFunc.InvSqrt) {
    assert(invsqrtTab.io.cs.asUInt === 0.U)
  }

  // --------------------------------------------------------------------------
  // reciprocal

  val recPre   = Module(new ReciprocalPreProcess (spec, polySpec, stage.preStage))
  val recTab   = Module(new ReciprocalTableCoeff (spec, polySpec, maxCbit))
  val recOther = Module(new ReciprocalOtherPath  (spec, polySpec, otherStage))
  val recPost  = Module(new ReciprocalPostProcess(spec, polySpec, stage.postStage))

  // atan2 uses reciprocal 1/max(x,y) to calculate min(x,y)/max(x,y).
  val recUseY = (io.sel === SelectFunc.ATan2Stage1) && yIsLarger

  recPre.io.en  := (io.sel === SelectFunc.Reciprocal) || (io.sel === SelectFunc.ATan2Stage1)
  recPre.io.x   := Mux(recUseY, ydecomp.io.decomp, xdecomp.io.decomp)
  // ------ Preprocess-Calculate ------
  recTab.io.en  := (selPCGapReg === SelectFunc.Reciprocal) ||
                   (selPCGapReg === SelectFunc.ATan2Stage1)
  recTab.io.adr := ShiftRegister(recPre.io.adr, pcGap)
  recOther.io.x := xdecPCGapReg

  when(selPCReg =/= SelectFunc.Reciprocal && selPCReg =/= SelectFunc.ATan2Stage1) {
    assert(recPre.io.adr === 0.U)
    assert(recPre.io.dx.getOrElse(0.U) === 0.U)
  }
  when(selPCGapReg =/= SelectFunc.Reciprocal && selPCGapReg =/= SelectFunc.ATan2Stage1) {
    assert(recTab.io.cs.asUInt === 0.U)
  }

  // --------------------------------------------------------------------------
  // sincos

  val sincosPre   = Module(new SinCosPreProcess (spec, polySpec, stage.preStage))
  val sincosTab   = Module(new SinCosTableCoeff (spec, polySpec, maxCbit))
  val sincosPost  = Module(new SinCosPostProcess(spec, polySpec, stage.postStage))

  sincosPre.io.en    := (io.sel === SelectFunc.Sin) || (io.sel === SelectFunc.Cos)
  sincosPre.io.isSin := (io.sel === SelectFunc.Sin)
  sincosPre.io.x     := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  sincosTab.io.en    := (selPCGapReg === SelectFunc.Sin) ||
                        (selPCGapReg === SelectFunc.Cos)
  sincosTab.io.adr   := ShiftRegister(sincosPre.io.adr, pcGap)

  when(selPCReg =/= SelectFunc.Sin && selPCReg =/= SelectFunc.Cos) {
    assert(sincosPre.io.adr === 0.U)
    assert(sincosPre.io.dx.getOrElse(0.U) === 0.U)
  }

  when(selPCGapReg =/= SelectFunc.Sin && selPCGapReg =/= SelectFunc.Cos) {
    assert(sincosTab.io.cs.asUInt === 0.U)
  }

  // --------------------------------------------------------------------------
  // atan2

  val atan2Stage1Pre   = Module(new ATan2Stage1PreProcess (spec, polySpec, stage.preStage))
  val atan2Stage1Other = Module(new ATan2Stage1OtherPath  (spec, polySpec, otherStage))
  val atan2Stage1Post  = Module(new ATan2Stage1PostProcess(spec, polySpec, stage.postStage))

  // atan2Stage1Pre checks if x and y are special values.
  // for calculation, reciprocal is re-used.
  atan2Stage1Pre.io.en := (io.sel === SelectFunc.ATan2Stage1)
  atan2Stage1Pre.io.x  := xdecomp.io.decomp
  atan2Stage1Pre.io.y  := ydecomp.io.decomp
  atan2Stage1Pre.io.yIsLarger := yIsLarger
  // ------ Preprocess-Calculate ------
  atan2Stage1Other.io.x  := xdecPCGapReg
  atan2Stage1Other.io.y  := ydecPCGapReg
  atan2Stage1Other.io.yIsLarger := yIsLargerPCGapReg

  val atan2Stage2Pre   = Module(new ATan2Stage2PreProcess (spec, polySpec, stage.preStage))
  val atan2Stage2Tab   = Module(new ATan2Stage2TableCoeff (spec, polySpec, maxCbit))
  val atan2Stage2Other = Module(new ATan2Stage2OtherPath  (spec, polySpec, otherStage))
  val atan2Stage2Post  = Module(new ATan2Stage2PostProcess(spec, polySpec, stage.postStage))
  atan2Stage2Pre.io.en  := (io.sel === SelectFunc.ATan2Stage2)
  atan2Stage2Pre.io.x   := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  atan2Stage2Tab.io.en  := (selPCGapReg === SelectFunc.ATan2Stage2)
  atan2Stage2Tab.io.adr := ShiftRegister(atan2Stage2Pre.io.adr, pcGap)
  atan2Stage2Other.io.x := xdecPCGapReg

  when(selPCReg =/= SelectFunc.ATan2Stage2) {
    assert(atan2Stage2Pre.io.adr === 0.U)
    assert(atan2Stage2Pre.io.dx.getOrElse(0.U) === 0.U)
  }
  when(selPCGapReg =/= SelectFunc.ATan2Stage2) {
    assert(atan2Stage2Tab.io.cs.asUInt === 0.U)
  }

  // ------------------------------------------------------------------------
  // atan related status register.
  // atan2 stage1 must save some values until atan2 stage2. So, after the pre-
  // process, it saves some flags to register.

  val atan2FlagReg = Reg(new ATan2Flags())
  // the timing is at the cycle when atan2Stage1Pre completes
  when(selPCReg === SelectFunc.ATan2Stage1) {
    atan2FlagReg.status  := Cat(yIsLargerPCReg, xdecPCReg.sgn)
    atan2FlagReg.special := atan2Stage1Pre.io.special
    atan2FlagReg.ysgn    := ydecPCReg.sgn
  }
  // The register is updated only when atan2stage1 is executed. That means that
  // we don't need to care about the timing here. ATan2Stage2 can only be
  // executed after the Stage1 because it uses the result of stage1 as its input.
  atan2Stage2Other.io.flags := atan2FlagReg
  atan2Stage2Post.io.flags  := atan2FlagReg

  // --------------------------------------------------------------------------
  // pow2/exp

  val expPre   = Module(new ExpPreProcess (spec, polySpec, stage.preStage, false))
  val expTab   = Module(new ExpTableCoeff (spec, polySpec, maxCbit))
  val expOther = Module(new ExpOtherPath  (spec, polySpec, otherStage))
  val expPost  = Module(new ExpPostProcess(spec, polySpec, stage.postStage))

  expPre.io.en     := (io.sel === SelectFunc.Exp) || (io.sel === SelectFunc.Pow2)
  expPre.io.isexp.get := (io.sel === SelectFunc.Exp)
  expPre.io.x      := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  expTab.io.en     := (selPCGapReg === SelectFunc.Exp) ||
                      (selPCGapReg === SelectFunc.Pow2)
  expTab.io.adr    := ShiftRegister(expPre.io.adr, pcGap)
  expOther.io.x    := xdecPCGapReg
  expOther.io.xint := ShiftRegister(expPre.io.xint, pcGap)
  expOther.io.xexd := ShiftRegister(expPre.io.xexd, pcGap)

  if(expPre.io.xfracLSBs.isDefined) {
    expOther.io.xfracLSBs.get := ShiftRegister(expPre.io.xfracLSBs.get, pcGap)
  }

  when(selPCReg =/= SelectFunc.Exp && selPCReg =/= SelectFunc.Pow2) {
    assert(expPre.io.adr === 0.U)
    assert(expPre.io.dx.getOrElse(0.U) === 0.U)
  }

  when(selPCGapReg =/= SelectFunc.Exp && selPCGapReg =/= SelectFunc.Pow2) {
    assert(expTab.io.cs.asUInt === 0.U)
  }

  // --------------------------------------------------------------------------
  // log2/ln

  val logPre   = Module(new LogPreProcess (spec, polySpec, stage.preStage))
  val logTab   = Module(new LogTableCoeff (spec, polySpec, maxCbit))
  val logOther = Module(new LogOtherPath  (spec, polySpec, otherStage))
  val logPost  = Module(new LogPostProcess(spec, polySpec, stage.postStage, false)) // can be both log2 and log

  logPre.io.en      := (io.sel === SelectFunc.Log) || (io.sel === SelectFunc.Log2)
  logPre.io.x       := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  val logPreExAdr = logPre.io.adr(logPre.io.adr.getWidth-1, logPre.io.adr.getWidth-2)
  logTab.io.en      := (selPCGapReg === SelectFunc.Log) ||
                       (selPCGapReg === SelectFunc.Log2)
  logTab.io.adr     := ShiftRegister(logPre.io.adr, pcGap)
  logOther.io.x     := xdecPCGapReg

  when(selPCReg =/= SelectFunc.Log && selPCReg =/= SelectFunc.Log2) {
    assert(logPre.io.adr === 0.U)
    assert(logPre.io.dx.getOrElse(0.U) === 0.U)
  }
  when(selPCGapReg =/= SelectFunc.Log && selPCGapReg =/= SelectFunc.Log2) {
    assert(logTab.io.cs.asUInt === 0.U)
  }

  // ------------------------------------------------------------------------
  //                  now we are here
  //                               |
  // .-------. +=================+ v .-------------.
  // |       |=+  .------------. +===|(dx) poly    |   .---------.
  // |pre-   |----| tableCoeff |=====|(cs) interp  |---| post-   |
  // |process|-.  '------------'     '-------------' .-| process |
  // '-------' |  .--------------------------------. | '---------'
  //           '--| non-table path (e.g. taylor)   |-'
  //              '--------------------------------'

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, maxCbit, stage.calcStage))

  if(order != 0) {
    val polynomialDx = ShiftRegister(sqrtPre.io.dx.get |
                                      recPre.io.dx.get |
                                   sincosPre.io.dx.get |
                                    acos1Pre.io.dx.get |
                                    acos2Pre.io.dx.get |
                                      recPre.io.dx.get |
                              atan2Stage2Pre.io.dx.get |
                                      expPre.io.dx.get |
                                      logPre.io.dx.get, pcGap + tcGap)
//     printf("polynomialEval.io.dx.get    = %b\n", polynomialDx)
    polynomialEval.io.dx.get := polynomialDx
  }

  val polynomialCoef = ShiftRegister(sqrtTab.io.cs.asUInt |
                     invsqrtTab.io.cs.asUInt |
                         recTab.io.cs.asUInt |
                      sincosTab.io.cs.asUInt |
                        acosTab.io.cs.asUInt |
                 atan2Stage2Tab.io.cs.asUInt |
                         expTab.io.cs.asUInt |
                         logTab.io.cs.asUInt, tcGap)

//   printf("polynomialEval.io.coeffs.cs = %b\n", polynomialCoef.asUInt)
  polynomialEval.io.coeffs := polynomialCoef.asTypeOf(new TableCoeffInput(maxCbit))

  // ------------------------------------------------------------------------
  //                                         we are here
  // .-------. .-----------------.   .-------------.  |
  // |       |-'  .------------. '---| polynomial  |  v .--------------.
  // |pre-   |----| tableCoeff |-----| interpolate |====|(zres)   post |
  // |process|-.  '------------'     '-------------'  +=|(zother) proc |
  // '-------' |  .--------------------------------.  I '--------------'
  //           '--| non-table path (e.g. taylor)   |==+
  //              '--------------------------------'

  val polynomialResultCPGapReg = ShiftRegister(polynomialEval.io.result, cpGap)

  acosPost.io.en     := (selCPGapReg === SelectFunc.ACos2)
  acosPost.io.x      := xdecCPGapReg
  acosPost.io.flags  := acosFlagReg
  acosPost.io.zres   := polynomialResultCPGapReg

  sqrtPost.io.en     := (selCPGapReg === SelectFunc.Sqrt || selCPGapReg === SelectFunc.ACos1)
  sqrtPost.io.zother := ShiftRegister(sqrtOther.io.zother, cpGap)
  sqrtPost.io.zres   := polynomialResultCPGapReg

  invsqrtPost.io.en     := (selCPGapReg === SelectFunc.InvSqrt)
  invsqrtPost.io.zother := ShiftRegister(invsqrtOther.io.zother, cpGap)
  invsqrtPost.io.zres   := polynomialResultCPGapReg

  recPost.io.en     := selCPGapReg === SelectFunc.Reciprocal
  recPost.io.zother := ShiftRegister(recOther.io.zother, cpGap)
  recPost.io.zres   := polynomialResultCPGapReg

  sincosPost.io.en   := (selCPGapReg === SelectFunc.Sin) ||
                        (selCPGapReg === SelectFunc.Cos)
  sincosPost.io.pre  := ShiftRegister(sincosPre.io.out, pcGap + tcGap + nCalcStage + cpGap)
  sincosPost.io.zres := polynomialResultCPGapReg

  atan2Stage1Post.io.en     := (selCPGapReg === SelectFunc.ATan2Stage1)
  atan2Stage1Post.io.zother := ShiftRegister(atan2Stage1Other.io.zother, cpGap)
  atan2Stage1Post.io.zres   := polynomialResultCPGapReg
  atan2Stage1Post.io.minxy  := Mux(yIsLargerCPGapReg, xdecCPGapReg, ydecCPGapReg)

  atan2Stage2Post.io.en     := (selCPGapReg === SelectFunc.ATan2Stage2)
  atan2Stage2Post.io.zother := ShiftRegister(atan2Stage2Other.io.zother, cpGap)
  atan2Stage2Post.io.zres   := polynomialResultCPGapReg
  atan2Stage2Post.io.x      := xdecCPGapReg

  if(expPre.io.xfracLSBs.isDefined) {
    expPost.io.zCorrCoef.get := ShiftRegister(expOther.io.zCorrCoef.get, cpGap)
  }
  expPost.io.en     := selCPGapReg === SelectFunc.Pow2 ||
                       selCPGapReg === SelectFunc.Exp
  expPost.io.zother := ShiftRegister(expOther.io.zother, cpGap)
  expPost.io.zres   := polynomialResultCPGapReg

  logPost.io.en     := selCPGapReg === SelectFunc.Log ||
                       selCPGapReg === SelectFunc.Log2
  logPost.io.zother := ShiftRegister(logOther.io.zother, cpGap)
  logPost.io.zres   := polynomialResultCPGapReg
  logPost.io.isln.get := selCPGapReg === SelectFunc.Log

  val z0 = sqrtPost.io.z        |
           invsqrtPost.io.z     |
           recPost.io.z         |
           sincosPost.io.z      |
           acosPost.io.z        |
           atan2Stage1Post.io.z |
           atan2Stage2Post.io.z |
           expPost.io.z         |
           logPost.io.z

  io.z := z0
}

class MathFuncUnitFP32( stage : MathFuncPipelineConfig )
    extends MathFunctions( RealSpec.Float32Spec, 2, 8, 3, stage) {
}

object MathFuncUnitFP32_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new MathFuncUnitFP32(MathFuncPipelineConfig.none)) ) )
}

class MathFuncUnitBF16( stage : MathFuncPipelineConfig )
    extends MathFunctions( RealSpec.BFloat16Spec, 0, 7, 1, stage) {
}

object MathFuncUnitBF16_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new MathFuncUnitBF16(MathFuncPipelineConfig.none)) ) )
}

