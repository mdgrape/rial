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
  val W           =  4
  val None        =  0.U(W.W)
  val Sqrt        =  1.U(W.W)
  val InvSqrt     =  2.U(W.W)
  val Reciprocal  =  3.U(W.W)
  val Sin         =  4.U(W.W)
  val Cos         =  5.U(W.W)
  val ACos        =  6.U(W.W)
  val ATan2Stage1 =  7.U(W.W)
  val ATan2Stage2 =  8.U(W.W)
  val Pow2        =  9.U(W.W)
  val Exp         = 10.U(W.W)
  val Log2        = 11.U(W.W)
  val Log         = 12.U(W.W)
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

class MathFunctions(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  assert(!spec.disableSign)

  // TODO: stage config
  val nStage = stage.total
  def getStage() = nStage

  val polySpec = new PolynomialSpec(spec, nOrder, adrW, extraBits,
    enableRangeCheck, enablePolynomialRounding)

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order = polySpec.order

  println(f"acos        cbits = ${ACosTableCoeff.getCBits(spec, polySpec)}")
  println(f"sqrt        cbits = ${SqrtTableCoeff.getCBits(spec, polySpec)}")
  println(f"invsqrt     cbits = ${InvSqrtTableCoeff.getCBits(spec, polySpec)}")
  println(f"rec         cbits = ${ReciprocalTableCoeff.getCBits(spec, polySpec)}")
  println(f"sincos      cbits = ${SinCosTableCoeff.getCBits(spec, polySpec)}")
  println(f"atan2Stage2 cbits = ${ATan2Stage2TableCoeff.getCBits(spec, polySpec)}")
  println(f"pow2        cbits = ${Pow2TableCoeff.getCBits(spec, polySpec)}")
  println(f"log2        cbits = ${Log2TableCoeff.getCBits(spec, polySpec)}")

  println(f"acos        calcW = ${ACosTableCoeff.getCalcW(spec, polySpec)}")
  println(f"sqrt        calcW = ${SqrtTableCoeff.getCalcW(spec, polySpec)}")
  println(f"invsqrt     calcW = ${InvSqrtTableCoeff.getCalcW(spec, polySpec)}")
  println(f"rec         calcW = ${ReciprocalTableCoeff.getCalcW(spec, polySpec)}")
  println(f"sincos      calcW = ${SinCosTableCoeff.getCalcW(spec, polySpec)}")
  println(f"atan2Stage2 calcW = ${ATan2Stage2TableCoeff.getCalcW(spec, polySpec)}")
  println(f"pow2        calcW = ${Pow2TableCoeff.getCalcW(spec, polySpec)}")
  println(f"log2        calcW = ${Log2TableCoeff.getCalcW(spec, polySpec)}")

  val maxCbit  = Seq(
    ACosTableCoeff.getCBits(spec, polySpec),
    SqrtTableCoeff.getCBits(spec, polySpec),
    InvSqrtTableCoeff.getCBits(spec, polySpec),
    ReciprocalTableCoeff.getCBits(spec, polySpec),
    SinCosTableCoeff.getCBits(spec, polySpec),
    ATan2Stage2TableCoeff.getCBits(spec, polySpec),
    Pow2TableCoeff.getCBits(spec, polySpec),
    Log2TableCoeff.getCBits(spec, polySpec)
    ).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  val maxCalcW = Seq(
    ACosTableCoeff.getCalcW(spec, polySpec),
    SqrtTableCoeff.getCalcW(spec, polySpec),
    InvSqrtTableCoeff.getCalcW(spec, polySpec),
    ReciprocalTableCoeff.getCalcW(spec, polySpec),
    SinCosTableCoeff.getCalcW(spec, polySpec),
    ATan2Stage2TableCoeff.getCalcW(spec, polySpec),
    Pow2TableCoeff.getCalcW(spec, polySpec),
    Log2TableCoeff.getCalcW(spec, polySpec)
    ).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  def getMaxCbit()  = maxCbit
  def getMaxCalcW() = maxCalcW

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

  val yIsLarger = io.x(spec.W-2, 0) < io.y(spec.W-2, 0) // without sign bit

  // .-------. .-----------------.   .-------------.
  // |       |-'  .------------. '---| polynomial  |   .---------.
  // |pre-   |====|(adr) table |-----| interpolate |---| post-   |
  // |process|=+  '------------'     '-------------' .-| process |
  // '-------' I  .--------------------------------. | '---------'
  //           +==|(x) non-table path (e.g. taylor)|-'
  //            ^ '--------------------------------'
  //            |
  // now we are here

  val acosPre   = Module(new ACosPreProcess (spec, polySpec, stage))
  val acosTab   = Module(new ACosTableCoeff (spec, polySpec, maxCbit, stage))
  val acosOther = Module(new ACosOtherPath  (spec, polySpec, stage))
  val acosPost  = Module(new ACosPostProcess(spec, polySpec, stage))

  acosPre.io.en  := (io.sel === SelectFunc.ACos)
  acosPre.io.x   := io.x
  acosTab.io.en  := (io.sel === SelectFunc.ACos) && (!acosPre.io.useSqrt)
  acosTab.io.adr := acosPre.io.adr
  acosOther.io.x       := xdecomp.io.decomp
  acosOther.io.useSqrt := acosPre.io.useSqrt
  acosOther.io.yex     := acosPre.io.yex
  acosOther.io.yman    := acosPre.io.yman

  when(io.sel =/= SelectFunc.ACos) {
    assert(acosPre.io.adr === 0.U)
    if(acosPre.io.dx.isDefined) {
      assert(acosPre.io.dx.get  === 0.U)
    }
  }
  when(io.sel =/= SelectFunc.ACos || acosPre.io.useSqrt) {
    assert(acosTab.io.cs.asUInt === 0.U)
  }

  val sqrtPre   = Module(new SqrtPreProcess (spec, polySpec, stage))
  val sqrtTab   = Module(new SqrtTableCoeff (spec, polySpec, maxCbit, stage))
  val sqrtOther = Module(new SqrtOtherPath  (spec, polySpec, stage))
  val sqrtPost  = Module(new SqrtPostProcess(spec, polySpec, stage))

  sqrtPre.io.en  := (io.sel === SelectFunc.Sqrt || io.sel === SelectFunc.InvSqrt)
  sqrtPre.io.x   := io.x
  // later we need to add a register here to make it pipe
  sqrtTab.io.en  := (io.sel === SelectFunc.Sqrt || acosPre.io.useSqrt)
  sqrtTab.io.adr := sqrtPre.io.adr | acosPre.io.adr
  sqrtOther.io.x := xdecomp.io.decomp

  when(!(io.sel === SelectFunc.Sqrt || io.sel === SelectFunc.InvSqrt)) {
    assert(sqrtPre.io.adr === 0.U)
    if(sqrtPre.io.dx.isDefined) {
      assert(sqrtPre.io.dx.get === 0.U)
    }
  }
  when(io.sel =/= SelectFunc.Sqrt) {
    assert(sqrtTab.io.cs.asUInt === 0.U || (io.sel === SelectFunc.ACos && acosPre.io.useSqrt))
  }

  val invsqrtTab   = Module(new InvSqrtTableCoeff (spec, polySpec, maxCbit, stage))
  val invsqrtOther = Module(new InvSqrtOtherPath  (spec, polySpec, stage))
  val invsqrtPost  = Module(new InvSqrtPostProcess(spec, polySpec, stage))

  invsqrtTab.io.en  := io.sel === SelectFunc.InvSqrt
  invsqrtTab.io.adr := sqrtPre.io.adr
  invsqrtOther.io.x := xdecomp.io.decomp

  when(io.sel =/= SelectFunc.InvSqrt) {
    assert(invsqrtTab.io.cs.asUInt === 0.U)
  }

  val recPre   = Module(new ReciprocalPreProcess (spec, polySpec, stage))
  val recTab   = Module(new ReciprocalTableCoeff (spec, polySpec, maxCbit, stage))
  val recOther = Module(new ReciprocalOtherPath  (spec, polySpec, stage))
  val recPost  = Module(new ReciprocalPostProcess(spec, polySpec, stage))

  // atan2 uses reciprocal 1/max(x,y) to calculate min(x,y)/max(x,y).
  val recUseY = (io.sel === SelectFunc.ATan2Stage1) && yIsLarger

  recPre.io.en  := (io.sel === SelectFunc.Reciprocal) || (io.sel === SelectFunc.ATan2Stage1)
  recPre.io.x   := Mux(recUseY, io.y, io.x)

  recTab.io.en  := (io.sel === SelectFunc.Reciprocal) || (io.sel === SelectFunc.ATan2Stage1)
  recTab.io.adr := recPre.io.adr
  recOther.io.x := xdecomp.io.decomp

  when(!((io.sel === SelectFunc.Reciprocal) || (io.sel === SelectFunc.ATan2Stage1))) {
    assert(recPre.io.adr === 0.U)
    if(recPre.io.dx.isDefined) {
      assert(recPre.io.dx.get  === 0.U)
    }
  }

  when(!((io.sel === SelectFunc.Reciprocal) || (io.sel === SelectFunc.ATan2Stage1))) {
    assert(recTab.io.cs.asUInt === 0.U)
  }

  val sincosPre   = Module(new SinCosPreProcess (spec, polySpec, stage))
  val sincosTab   = Module(new SinCosTableCoeff (spec, polySpec, maxCbit, stage))
  val sincosOther = Module(new SinCosOtherPath  (spec, polySpec, stage))
  val sincosPost  = Module(new SinCosPostProcess(spec, polySpec, stage))

  sincosPre.io.en    := (io.sel === SelectFunc.Sin) || (io.sel === SelectFunc.Cos)
  sincosPre.io.isSin := (io.sel === SelectFunc.Sin)
  sincosPre.io.x     := io.x

  sincosTab.io.en    := (io.sel === SelectFunc.Sin) || (io.sel === SelectFunc.Cos)
  sincosTab.io.adr   := sincosPre.io.adr

  sincosOther.io.xConverted := sincosPre.io.xConverted
  sincosOther.io.x          := xdecomp.io.decomp

  when(io.sel =/= SelectFunc.Sin && io.sel =/= SelectFunc.Cos) {
    assert(sincosPre.io.adr === 0.U)
    if(sincosPre.io.dx.isDefined) {
      assert(sincosPre.io.dx.get  === 0.U)
    }
  }

  when(io.sel =/= SelectFunc.Sin && io.sel =/= SelectFunc.Cos) {
    assert(sincosTab.io.cs.asUInt === 0.U)
  }

  val atan2Stage1Pre   = Module(new ATan2Stage1PreProcess (spec, polySpec, stage))
  val atan2Stage1Other = Module(new ATan2Stage1OtherPath  (spec, polySpec, stage))
  val atan2Stage1Post  = Module(new ATan2Stage1PostProcess(spec, polySpec, stage))

  // atan2Stage1Pre checks if x and y are special values.
  // for calculation, reciprocal is re-used.
  atan2Stage1Pre.io.en:= (io.sel === SelectFunc.ATan2Stage1)
  atan2Stage1Pre.io.x := xdecomp.io.decomp
  atan2Stage1Pre.io.y := ydecomp.io.decomp

  atan2Stage1Other.io.x := xdecomp.io.decomp
  atan2Stage1Other.io.y := ydecomp.io.decomp
  atan2Stage1Other.io.yIsLarger := yIsLarger

  val atan2Stage2Pre   = Module(new ATan2Stage2PreProcess (spec, polySpec, stage))
  val atan2Stage2Tab   = Module(new ATan2Stage2TableCoeff (spec, polySpec, maxCbit, stage))
  val atan2Stage2Other = Module(new ATan2Stage2OtherPath  (spec, polySpec, stage))
  val atan2Stage2Post  = Module(new ATan2Stage2PostProcess(spec, polySpec, stage))
  atan2Stage2Pre.io.en  := (io.sel === SelectFunc.ATan2Stage2)
  atan2Stage2Pre.io.x   := io.x
  atan2Stage2Tab.io.en  := (io.sel === SelectFunc.ATan2Stage2)
  atan2Stage2Tab.io.adr := atan2Stage2Pre.io.adr
  atan2Stage2Other.io.x := xdecomp.io.decomp

  when(!(io.sel === SelectFunc.ATan2Stage2)) {
    assert(atan2Stage2Pre.io.adr === 0.U)
    if(atan2Stage2Pre.io.dx.isDefined) {
      assert(atan2Stage2Pre.io.dx.get  === 0.U)
    }
  }
  when(io.sel =/= SelectFunc.ATan2Stage2) {
    assert(atan2Stage2Tab.io.cs.asUInt === 0.U)
  }


  // ------------------------------------------------------------------------
  // atan related status register

  val atan2FlagReg = Reg(new ATan2Flags())
  when(io.sel === SelectFunc.ATan2Stage1) {
    // check special values ... TODO: need to consider the delay in sel and atan2Stage1PreProcess
    atan2FlagReg.status  := Cat(yIsLarger, io.x(spec.W-1))
    atan2FlagReg.special := atan2Stage1Pre.io.special
    atan2FlagReg.ysgn    :=     io.y(spec.W-1)
  }
  atan2Stage2Other.io.flags := atan2FlagReg

  val pow2Pre   = Module(new Pow2PreProcess (spec, polySpec, stage))
  val pow2Tab   = Module(new Pow2TableCoeff (spec, polySpec, maxCbit, stage))
  val pow2Other = Module(new Pow2OtherPath  (spec, polySpec, stage))
  val pow2Post  = Module(new Pow2PostProcess(spec, polySpec, stage))

  pow2Pre.io.en     := (io.sel === SelectFunc.Exp) || (io.sel === SelectFunc.Pow2)
  pow2Pre.io.isexp  := (io.sel === SelectFunc.Exp)
  pow2Pre.io.x      := io.x
  pow2Tab.io.en     := (io.sel === SelectFunc.Exp) || (io.sel === SelectFunc.Pow2)
  pow2Tab.io.adr    := pow2Pre.io.adr
  pow2Other.io.x    := xdecomp.io.decomp
  pow2Other.io.xint := pow2Pre.io.xint
  pow2Other.io.xexd := pow2Pre.io.xexd

  if(pow2Pre.io.xfracLSBs.isDefined) {
    pow2Other.io.xfracLSBs.get := pow2Pre.io.xfracLSBs.get
  }

  when(!((io.sel === SelectFunc.Exp) || (io.sel === SelectFunc.Pow2))) {
    assert(pow2Pre.io.adr === 0.U)
    if(pow2Pre.io.dx.isDefined) {
      assert(pow2Pre.io.dx.get  === 0.U)
    }
  }

  when(!((io.sel === SelectFunc.Exp) || (io.sel === SelectFunc.Pow2))) {
    assert(pow2Tab.io.cs.asUInt === 0.U)
  }

  val log2Pre   = Module(new Log2PreProcess (spec, polySpec, stage))
  val log2Tab   = Module(new Log2TableCoeff (spec, polySpec, maxCbit, stage))
  val log2Other = Module(new Log2OtherPath  (spec, polySpec, stage))
  val log2Post  = Module(new Log2PostProcess(spec, polySpec, stage))

  log2Pre.io.en      := (io.sel === SelectFunc.Log) || (io.sel === SelectFunc.Log2)
  log2Pre.io.x       := io.x
  log2Tab.io.en      := (io.sel === SelectFunc.Log) || (io.sel === SelectFunc.Log2)
  log2Tab.io.adr     := log2Pre.io.adr
  log2Other.io.x     := xdecomp.io.decomp
  log2Other.io.exadr := log2Pre.io.adr(log2Pre.io.adr.getWidth-1, log2Pre.io.adr.getWidth-2)
  log2Post.io.x      := xdecomp.io.decomp
  log2Post.io.exadr  := log2Pre.io.adr(log2Pre.io.adr.getWidth-1, log2Pre.io.adr.getWidth-2)
  log2Post.io.xmanbp := log2Other.io.xmanbp

  when(!((io.sel === SelectFunc.Log) || (io.sel === SelectFunc.Log2))) {
    assert(log2Pre.io.adr === 0.U)
    if(log2Pre.io.dx.isDefined) {
      assert(log2Pre.io.dx.get  === 0.U)
    }
  }
  when(!((io.sel === SelectFunc.Log) || (io.sel === SelectFunc.Log2))) {
    assert(log2Tab.io.cs.asUInt === 0.U)
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

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, maxCbit, stage))

  if(order != 0) {
    polynomialEval.io.dx.get := sqrtPre.io.dx.get |
                                 recPre.io.dx.get |
                              sincosPre.io.dx.get |
                                acosPre.io.dx.get |
                                 recPre.io.dx.get |
                         atan2Stage2Pre.io.dx.get |
                                pow2Pre.io.dx.get |
                                log2Pre.io.dx.get
  }

  polynomialEval.io.coeffs.cs := (sqrtTab.io.cs.cs.asUInt |
                               invsqrtTab.io.cs.cs.asUInt |
                                   recTab.io.cs.cs.asUInt |
                                sincosTab.io.cs.cs.asUInt |
                                  acosTab.io.cs.cs.asUInt |
                           atan2Stage2Tab.io.cs.cs.asUInt |
                                  pow2Tab.io.cs.cs.asUInt |
                                  log2Tab.io.cs.cs.asUInt
                                  ).asTypeOf(new MixedVec(maxCbit.map{w => UInt(w.W)}))

  //                                         we are here
  // .-------. .-----------------.   .-------------.  |
  // |       |-'  .------------. '---| polynomial  |  v .--------------.
  // |pre-   |----| tableCoeff |-----| interpolate |====|(zres)   post |
  // |process|-.  '------------'     '-------------'  +=|(zother) proc |
  // '-------' |  .--------------------------------.  I '--------------'
  //           '--| non-table path (e.g. taylor)   |==+
  //              '--------------------------------'

  sqrtPost.io.en     := (io.sel === SelectFunc.Sqrt)
  sqrtPost.io.zother := sqrtOther.io.zother
  sqrtPost.io.zres   := polynomialEval.io.result

  invsqrtPost.io.en     := (io.sel === SelectFunc.InvSqrt)
  invsqrtPost.io.zother := invsqrtOther.io.zother
  invsqrtPost.io.zres   := polynomialEval.io.result

  recPost.io.en     := io.sel === SelectFunc.Reciprocal
  recPost.io.zother := recOther.io.zother
  recPost.io.zres   := polynomialEval.io.result

  sincosPost.io.en     := (io.sel === SelectFunc.Sin || io.sel === SelectFunc.Cos)
  sincosPost.io.zother := sincosOther.io.zother
  sincosPost.io.zres   := polynomialEval.io.result

  acosPost.io.en     := (io.sel === SelectFunc.ACos)
  acosPost.io.zother := acosOther.io.zother
  acosPost.io.zres   := polynomialEval.io.result

  atan2Stage1Post.io.en     := (io.sel === SelectFunc.ATan2Stage1)
  atan2Stage1Post.io.zother := atan2Stage1Other.io.zother
  atan2Stage1Post.io.zres   := polynomialEval.io.result
  atan2Stage1Post.io.minxy  := Mux(yIsLarger, xdecomp.io.decomp, ydecomp.io.decomp)

  atan2Stage2Post.io.en     := (io.sel === SelectFunc.ATan2Stage2)
  atan2Stage2Post.io.zother := atan2Stage2Other.io.zother
  atan2Stage2Post.io.zres   := polynomialEval.io.result
  atan2Stage2Post.io.flags  := atan2FlagReg // TODO keep the value in frag reg

  if(pow2Pre.io.xfracLSBs.isDefined) {
    pow2Post.io.zCorrCoef.get := pow2Other.io.zCorrCoef.get
  }
  pow2Post.io.en     := (io.sel === SelectFunc.Pow2 || io.sel === SelectFunc.Exp)
  pow2Post.io.zother := pow2Other.io.zother
  pow2Post.io.zres   := polynomialEval.io.result

  log2Post.io.en     := (io.sel === SelectFunc.Log || io.sel === SelectFunc.Log2)
  log2Post.io.zother := log2Other.io.zother
  log2Post.io.zres   := polynomialEval.io.result
  log2Post.io.isln   := io.sel === SelectFunc.Log

  val z0 = sqrtPost.io.z        |
           invsqrtPost.io.z     |
           recPost.io.z         |
           sincosPost.io.z      |
           acosPost.io.z        |
           atan2Stage1Post.io.z |
           atan2Stage2Post.io.z |
           pow2Post.io.z        |
           log2Post.io.z

  io.z := ShiftRegister(z0, stage.total)
}

class MathFuncUnit( stage : PipelineStageConfig )
    extends MathFunctions( RealSpec.Float32Spec, 2, 8, 3, stage) {
}

object MathFuncUnit_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new MathFuncUnit(PipelineStageConfig.none)) ) )
}

