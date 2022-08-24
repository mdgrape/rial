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


/* Enumerator to represent which function is used in [[rial.math.MathFunctions]]. */
object FuncKind extends Enumeration {
  type FuncKind = Value
  val Sqrt, InvSqrt, Reciprocal, Sin, Cos, ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2, Exp, Log = Value
}


/** A Config class for [[rial.math.MathFunctions]] Module.
 *
 *  @constructor create a new MathFuncConfig.
 *  @param funcs the list of functions that should be supported.
 */
class MathFuncConfig(
  val funcs: Seq[FuncKind.FuncKind]
) {
  assert(funcs.length > 0, "At least one function should be supported")
  import FuncKind._

  /** Checks if a function is supported.
   *
   *  @param fn An enumerator of the function.
   *  @return true if the function is supported. false if not.
   */
  def has(fn: FuncKind): Boolean = {
    funcs.exists(_==fn)
  }

  if(has(ACosPhase1) || has(ACosPhase2)) {
    assert(has(ACosPhase1) && has(ACosPhase2), "ACos requires both phase 1 and 2.")

    // since acosPhase1 uses sqrt table, we auto-generate sqrt.
    // and invsqrt shares the preproc with sqrt, it also auto-generate invsqrt.
    // we can remove this assumption in principle, but the implementation will
    // be painful and redundant.
    assert(has(Sqrt) && has(InvSqrt), """TODO: remove this assumption later.
      | Since acosPhase1 uses sqrt table, we auto-generate sqrt.
      | And invsqrt shares the preproc with sqrt, it also auto-generate invsqrt.
      | In principle, we can remove this assumption, but the implementation will
      | be painful and redundant.
      """.stripMargin)
  }
  if(has(ATan2Phase1) || has(ATan2Phase2)) {
    assert(has(ATan2Phase1) && has(ATan2Phase2), "ATan2 requires both phase 1 and 2.")

    assert(has(Reciprocal), """TODO: remove this assumption later.
      | Since atan2Phase1 uses reciprocal table, we auto-generate reciprocal.
      | In principle, we can remove this assumption, but the implementation will
      | be painful and redundant.
      """.stripMargin)
  }

  /** The width of UInt to represent function select signal.
   */
  val signalW = log2Up(funcs.length)

  /** Returns function select signal. It starts from 1.
   *  If the passed function is not supported, returns 0.
   *
   *  @param fn An enumerator of the function
   *  @return the signal that corresponds to the function
   */
  def signal(fn: FuncKind): UInt = {
    (funcs.indexWhere(_==fn) + 1).U(signalW.W)
  }

  /** Returns function select signal that runs no function.
   *
   *  @return the signal that corresponds to no function.
   */
  def signalNone(): UInt = {
    0.U(signalW.W)
  }
}

/** Factory for [[rial.math.MathFuncConfig]].
 */
object MathFuncConfig {
  import FuncKind._

  /* All functions are supported.
   */
  val all = new MathFuncConfig(Seq(Sqrt, InvSqrt, Reciprocal, Sin, Cos,
    ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2, Exp, Log))
}

/** A Bundle that is returned from [[rial.math.DecomposeReal]].
 *
 * It decomposes a UInt that represents real value using [[rial.arith.RealSpec]]
 * and check if the value is zero, inf, or nan.
 *
 * @param spec A [[rial.arith.RealSpec]] corresponding to the input floating point number.
 */
class DecomposedRealOutput(val spec: RealSpec) extends Bundle {
  val sgn  = Output(UInt(1.W))
  val ex   = Output(UInt(spec.exW.W))
  val man  = Output(UInt(spec.manW.W))
  val zero = Output(Bool())
  val inf  = Output(Bool())
  val nan  = Output(Bool())
}
/** A Module that decomposes a UInt that represents Real.
 *
 * @param spec A [[rial.arith.RealSpec]] corresponding to the input floating point number.
 */
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

/** Config class to set pipeline stages in a [[rial.math.MathFunctions]] module.
 *
 * Overview:
 *
 * {{{
 * //             .--preStage    .--tableCalcGap .-- calcPostGap
 * //             |      .--preCalcGap   .--calcStage    .--postStage
 * //         ____|____  |       |  _____|_____  |  _____|_____
 * //        '         ' '       ' '           ' ' '           '
 * //       .----.-.----.-.-----.-.-----.-.-----.-.-----.-.-----.
 * //       |    |v|    |v|table|v|Calc1|v|Calc2|v|     |v|     |
 * // in -> |Pre1| |Pre2| :-----'-'-----'-'-----: |Post1| |Post2| -> out
 * //       |    | |    | |      non-table      | |     | |     |
 * //       '----'-'----'-'---------------------'-'-----'-'-----'
 * }}}
 *
 * TODO: consider setting nStage for each function. (like, sqrt does not need
 *       multiple cycles in its preprocess, but sincos may need.)
 *
 * @constructor create a new MathFuncConfig.
 * @param preStage     clock cycles in preprocess.
 * @param calcStage    clock cycles in table/polynomial and non-table path.
 * @param postStage    clock cycles in postprocess.
 * @param preCalcGap   if true, add register between preprocess and calculation stage
 * @param tableCalcGap if true, add register between table and calculation stage (+1 to calcStage for OtherPath)
 * @param calcPostGap  if true, add register between calculation and postprocess stage
 *
 */
class MathFuncPipelineConfig(
  val preStage:     PipelineStageConfig,
  val calcStage:    PipelineStageConfig,
  val postStage:    PipelineStageConfig,
  val preCalcGap:   Boolean,
  val tableCalcGap: Boolean,
  val calcPostGap:  Boolean,
  ) {

  /** Calculates the total latency in clock cycles.
   */
  def total = {
    preStage.total +
    calcStage.total +
    postStage.total +
    (if(preCalcGap)   {1} else {0}) +
    (if(tableCalcGap) {1} else {0}) +
    (if(calcPostGap)  {1} else {0})
  }
  /** Generates a string that represents the current config.
   */
  def getString = {
    f"pre: ${preStage.getString}, " +
    f"calc: ${calcStage.getString}, " +
    f"post: ${postStage.getString}, " +
    f"pre-calc: ${preCalcGap}, table-calc: ${tableCalcGap}, calc-post: ${calcPostGap}"
  }
}

/** Factory for [[rial.math.MathFuncPipelineConfig]].
 */
object MathFuncPipelineConfig {

  /** constructs a config corresponding to the single cycle mathfunc.
   */
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

/** A module that calculates math functions.
 *
 * # Overview
 *
 * It takes two floating point numbers, x and y, and returns one floating point
 * number, z. Only ATan2 uses both x and y. Others just ignores y.
 *
 * TODO: make io.y optional and disable when ATan2 is not required.
 *
 * {{{
 * //                                     table/Calc
 * //                                          |
 * //                .------.  _  .---------.  _  .------------.
 * // sel -----------|      |-|v|-'.-------.'-|v|-| chebyshev  |    _  .-------.
 * //    .---------. | pre- |-| |--| table |--| |-| polynomial |---|v|-| post- |-> z
 * // x -|decompose|-| proc | | |  '-------'  '-' '------------' .-| |-| proc  |
 * // y -|         |-|      |-| |-. .--------------------------. | '_' '-------'
 * //    '---------' '------' '-' '-|  non-table (taylor etc)  |-'  .
 * //                          .    '--------------------------'    |
 * //                          |                                    |
 * //                   Preproc/Calc                          Calc/Postproc
 * }}}
 *
 *
 * @constructor     Constructs a new MathFunctions module.
 * @param fncfg     determines which functions are supported.
 * @param spec      determines the spec of input/output floating point number.
 * @param nOrder    determines the order of polynomial.
 * @param adrW      determines the default bit width of the tables.
 * @param extraBits determines the extra bits of polynomial results.
 * @param stage     determines the pipeline stages of this module.
 * @param dxW0      determines the default bit width of the polynomial input.
 *
 * table address and input bit width of polynomial depend on a function.
 * For detail, see the corresponding math function module.
 *
 */
class MathFunctions(
  val fncfg: MathFuncConfig,
  val spec : RealSpec,
  val nOrder: Int, val adrW : Int, val extraBits : Int,
  val stage : MathFuncPipelineConfig,
  val dxW0 : Option[Int] = None,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  import FuncKind._;

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

  /** Returns the bit width of the coefficients of polynomial used by the function.
   *  Used to determine maximum bit width of [[rial.math.PolynomialEval]] module
   *  embedded in a [[rial.math.MathFunctions]]
   *
   *  @return bit width of coefficients used in the polynomial of the function.
   */
  def getCBit(fn: FuncKind.FuncKind): Seq[Int] = {
    import FuncKind._
    fn match {
      case Sqrt        => SqrtTableCoeff.getCBits(spec, polySpec)
      case InvSqrt     => InvSqrtTableCoeff.getCBits(spec, polySpec)
      case Reciprocal  => ReciprocalTableCoeff.getCBits(spec, polySpec)
      case Sin         => SinCosTableCoeff.getCBits(spec, polySpec)
      case Cos         => SinCosTableCoeff.getCBits(spec, polySpec)
      case ACosPhase1  => SqrtTableCoeff.getCBits(spec, polySpec)
      case ACosPhase2  => ACosTableCoeff.getCBits(spec, polySpec)
      case ATan2Phase1 => ReciprocalTableCoeff.getCBits(spec, polySpec)
      case ATan2Phase2 => ATan2Stage2TableCoeff.getCBits(spec, polySpec)
      case Exp         => ExpTableCoeff.getCBits(spec, polySpec)
      case Log         => LogTableCoeff.getCBits(spec, polySpec)
    }
  }

  /** Returns the bit width of the temporary of polynomial used by the function.
   *  Used to determine maximum bit width of [[rial.math.PolynomialEval]] module
   *  embedded in a [[rial.math.MathFunctions]]
   *
   *  @return bit width of temporary used in the polynomial of the function.
   */
  def getCalcW(fn: FuncKind.FuncKind): Seq[Int] = {
    import FuncKind._
    fn match {
      case Sqrt        => SqrtTableCoeff.getCalcW(spec, polySpec)
      case InvSqrt     => InvSqrtTableCoeff.getCalcW(spec, polySpec)
      case Reciprocal  => ReciprocalTableCoeff.getCalcW(spec, polySpec)
      case Sin         => SinCosTableCoeff.getCalcW(spec, polySpec)
      case Cos         => SinCosTableCoeff.getCalcW(spec, polySpec)
      case ACosPhase1  => SqrtTableCoeff.getCalcW(spec, polySpec)
      case ACosPhase2  => ACosTableCoeff.getCalcW(spec, polySpec)
      case ATan2Phase1 => ReciprocalTableCoeff.getCalcW(spec, polySpec)
      case ATan2Phase2 => ATan2Stage2TableCoeff.getCalcW(spec, polySpec)
      case Exp         => ExpTableCoeff.getCalcW(spec, polySpec)
      case Log         => LogTableCoeff.getCalcW(spec, polySpec)
    }
  }

  val maxCbit  = fncfg.funcs.map(f => getCBit(f)).
    reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )
  val maxCalcW = fncfg.funcs.map(f => getCalcW(f)).
    reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  def getMaxCbit  = maxCbit
  def getMaxCalcW = maxCalcW

  println(f"[${if(fncfg.has(ACosPhase2         )){"x"}else{" "}}] acos    cbits = ${ACosTableCoeff       .getCBits(spec, polySpec)} calcW = ${ACosTableCoeff       .getCalcW(spec, polySpec)}")
  println(f"[${if(fncfg.has(Sqrt               )){"x"}else{" "}}] sqrt    cbits = ${SqrtTableCoeff       .getCBits(spec, polySpec)} calcW = ${SqrtTableCoeff       .getCalcW(spec, polySpec)}")
  println(f"[${if(fncfg.has(InvSqrt            )){"x"}else{" "}}] invsqrt cbits = ${InvSqrtTableCoeff    .getCBits(spec, polySpec)} calcW = ${InvSqrtTableCoeff    .getCalcW(spec, polySpec)}")
  println(f"[${if(fncfg.has(Reciprocal         )){"x"}else{" "}}] rec     cbits = ${ReciprocalTableCoeff .getCBits(spec, polySpec)} calcW = ${ReciprocalTableCoeff .getCalcW(spec, polySpec)}")
  println(f"[${if(fncfg.has(Sin)||fncfg.has(Cos)){"x"}else{" "}}] sincos  cbits = ${SinCosTableCoeff     .getCBits(spec, polySpec)} calcW = ${SinCosTableCoeff     .getCalcW(spec, polySpec)}")
  println(f"[${if(fncfg.has(ATan2Phase2        )){"x"}else{" "}}] atan2   cbits = ${ATan2Stage2TableCoeff.getCBits(spec, polySpec)} calcW = ${ATan2Stage2TableCoeff.getCalcW(spec, polySpec)}")
  println(f"[${if(fncfg.has(Exp                )){"x"}else{" "}}] exp     cbits = ${ExpTableCoeff        .getCBits(spec, polySpec)} calcW = ${ExpTableCoeff        .getCalcW(spec, polySpec)}")
  println(f"[${if(fncfg.has(Log                )){"x"}else{" "}}] log     cbits = ${LogTableCoeff        .getCBits(spec, polySpec)} calcW = ${LogTableCoeff        .getCalcW(spec, polySpec)}")
  println(f"maximum cbits = ${maxCbit} calcW = ${maxCalcW}")

  val io = IO(new Bundle {
    val sel = Input(UInt(fncfg.signalW.W))
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
  // Output float.
  //
  // later we will insert a value to the map.

  val zs = fncfg.funcs.map( fn => { fn -> Wire(UInt(spec.W.W)) }).toMap
  zs.values.foreach(_ := 0.U) // init
  val z0 = zs.values.reduce(_|_)
  io.z := z0

  // ==========================================================================
  // Polynomial Evaluator

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, maxCbit, stage.calcStage))

  val tableCoeffWidth = maxCbit.reduce(_+_)
  val polynomialCoefs = fncfg.funcs.map( fn => {
    fn -> Wire(UInt(tableCoeffWidth.W))
  }).toMap
  polynomialCoefs.values.foreach(_:=0.U)

  val polynomialCoef = ShiftRegister(polynomialCoefs.values.reduce(_|_), tcGap)
  polynomialEval.io.coeffs := polynomialCoef.asTypeOf(new TableCoeffInput(maxCbit))

  val polynomialDxs = if (order == 0) { None } else { Some(fncfg.funcs.map(fn => {
      fn -> Wire(UInt(polySpec.dxW.W))
    }).toMap)
  }
  if(order != 0) {
    polynomialCoefs.values.foreach(_:=0.U) // init

    val polynomialDx = ShiftRegister(polynomialDxs.get.values.reduce(_|_), pcGap + tcGap)
//     printf("polynomialEval.io.dx.get    = %b\n", polynomialDx)
    polynomialEval.io.dx.get := polynomialDx
  }

  val polynomialResultCPGapReg = ShiftRegister(polynomialEval.io.result, cpGap)


  // ==========================================================================
  // ACos

  assert(fncfg.signal(ACosPhase1) =/= 0.U)
  assert(fncfg.signal(ACosPhase2) =/= 0.U)

  val acosFlagReg = Reg(new ACosFlags())

  val acos1Pre  = Module(new ACosStage1PreProcess(spec, polySpec, stage.preStage))
  val acos2Pre  = Module(new ACosStage2PreProcess(spec, polySpec, stage.preStage))
  val acosTab   = Module(new ACosTableCoeff (spec, polySpec, maxCbit))
  val acosPost  = Module(new ACosPostProcess(spec, polySpec, stage.postStage))

  val acos1PreAdrPCGapReg = ShiftRegister(acos1Pre.io.adr, pcGap)

  acos1Pre.io.en := io.sel === fncfg.signal(ACosPhase1)
  acos1Pre.io.x  := xdecomp.io.decomp
  acos2Pre.io.en := io.sel === fncfg.signal(ACosPhase2)
  acos2Pre.io.x  := xdecomp.io.decomp

  if(order != 0) {
    polynomialDxs.get(ACosPhase1)  := acos1Pre.io.dx.get
    polynomialDxs.get(ACosPhase2)  := acos2Pre.io.dx.get
  }

  // ------ Preprocess-Calculate ------

  acosTab.io.en  := selPCGapReg === fncfg.signal(ACosPhase2)
  acosTab.io.adr := ShiftRegister(acos2Pre.io.adr, pcGap)

  polynomialCoefs(ACosPhase2)  := acosTab.io.cs.asUInt

  acosPost.io.en     := (selCPGapReg === fncfg.signal(ACosPhase2))
  acosPost.io.x      := xdecCPGapReg
  acosPost.io.flags  := acosFlagReg
  acosPost.io.zres   := polynomialResultCPGapReg

  zs(ACosPhase2)  := acosPost.io.z

  when(selPCReg === fncfg.signal(ACosPhase1)) {
    acosFlagReg := acos1Pre.io.special
  }

  // after preprocess
  // XXX Since `PCReg`s are delayed by nPreStage, the timing is the same as acosPre output.
  when(selPCReg =/= fncfg.signal(ACosPhase1)) {
    assert(acos1Pre.io.adr === 0.U)
    assert(acos1Pre.io.dx.getOrElse(0.U) === 0.U)
  }
  when(selPCReg =/= fncfg.signal(ACosPhase2)) {
    assert(acos2Pre.io.adr === 0.U)
    assert(acos2Pre.io.dx.getOrElse(0.U) === 0.U)
  }

  when(selPCGapReg =/= fncfg.signal(ACosPhase2)) {
    assert(acosTab.io.cs.asUInt === 0.U)
  }

  // ==========================================================================
  // sqrt

  val sqrtPre   = Module(new SqrtPreProcess (spec, polySpec, stage.preStage))
  val sqrtTab   = Module(new SqrtTableCoeff (spec, polySpec, maxCbit))
  val sqrtOther = Module(new SqrtOtherPath  (spec, polySpec, otherStage))
  val sqrtPost  = Module(new SqrtPostProcess(spec, polySpec, stage.postStage))

  val sqrtPreAdrPCGapReg = ShiftRegister(sqrtPre.io.adr, pcGap)

  sqrtPre.io.en  := io.sel === fncfg.signal(Sqrt) || io.sel === fncfg.signal(InvSqrt)
  sqrtPre.io.x   := xdecomp.io.decomp
  if(order != 0) {
    polynomialDxs.get(Sqrt) := sqrtPre.io.dx.get
  }

  // ------ Preprocess-Calculate ------

  sqrtTab.io.en  := selPCGapReg === fncfg.signal(Sqrt) || selPCGapReg === fncfg.signal(ACosPhase1)
  sqrtTab.io.adr := sqrtPreAdrPCGapReg | acos1PreAdrPCGapReg

  // redundant...
  polynomialCoefs(Sqrt)        := sqrtTab.io.cs.asUInt
  polynomialCoefs(ACosPhase1)  := sqrtTab.io.cs.asUInt

  sqrtOther.io.x := Mux(selPCGapReg === fncfg.signal(Sqrt), xdecPCGapReg, ShiftRegister(acos1Pre.io.y, pcGap))

  sqrtPost.io.en     := (selCPGapReg === fncfg.signal(Sqrt) || selCPGapReg === fncfg.signal(ACosPhase1))
  sqrtPost.io.zother := ShiftRegister(sqrtOther.io.zother, cpGap)
  sqrtPost.io.zres   := polynomialResultCPGapReg

  // redundant...
  zs(ACosPhase1)  := sqrtPost.io.z
  zs(Sqrt) := sqrtPost.io.z

  // after preprocess
  when(selPCReg =/= fncfg.signal(Sqrt) && selPCReg =/= fncfg.signal(InvSqrt)) {
    assert(sqrtPre.io.adr === 0.U)
    assert(sqrtPre.io.dx.getOrElse(0.U) === 0.U)
  }
  when(selPCGapReg =/= fncfg.signal(Sqrt)) {
    assert(sqrtTab.io.cs.asUInt === 0.U || selPCGapReg === fncfg.signal(ACosPhase1))
  }

  // ==========================================================================
  // invsqrt
  val invsqrtTab   = Module(new InvSqrtTableCoeff (spec, polySpec, maxCbit))
  val invsqrtOther = Module(new InvSqrtOtherPath  (spec, polySpec, otherStage))
  val invsqrtPost  = Module(new InvSqrtPostProcess(spec, polySpec, stage.postStage))

  // (preprocess is same as sqrt)

  if(order != 0) {
    // redundant...
    polynomialDxs.get(InvSqrt) := sqrtPre.io.dx.get
  }

  invsqrtTab.io.en  := selPCGapReg === fncfg.signal(InvSqrt)
  invsqrtTab.io.adr := sqrtPreAdrPCGapReg

  polynomialCoefs(InvSqrt) := invsqrtTab.io.cs.asUInt

  invsqrtOther.io.x := xdecPCGapReg

  invsqrtPost.io.en     := (selCPGapReg === fncfg.signal(InvSqrt))
  invsqrtPost.io.zother := ShiftRegister(invsqrtOther.io.zother, cpGap)
  invsqrtPost.io.zres   := polynomialResultCPGapReg

  zs(InvSqrt) := invsqrtPost.io.z

  when(selPCGapReg =/= fncfg.signal(InvSqrt)) {
    assert(invsqrtTab.io.cs.asUInt === 0.U)
  }

  // ==========================================================================
  // reciprocal

  val recPre   = Module(new ReciprocalPreProcess (spec, polySpec, stage.preStage))
  val recTab   = Module(new ReciprocalTableCoeff (spec, polySpec, maxCbit))
  val recOther = Module(new ReciprocalOtherPath  (spec, polySpec, otherStage))
  val recPost  = Module(new ReciprocalPostProcess(spec, polySpec, stage.postStage))

  // atan2 uses reciprocal 1/max(x,y) to calculate min(x,y)/max(x,y).
  val recUseY = (io.sel === fncfg.signal(ATan2Phase1)) && yIsLarger

  recPre.io.en  := (io.sel === fncfg.signal(Reciprocal)) || (io.sel === fncfg.signal(ATan2Phase1))
  recPre.io.x   := Mux(recUseY, ydecomp.io.decomp, xdecomp.io.decomp)

  if(order != 0) {
    polynomialDxs.get(Reciprocal)  := recPre.io.dx.get
  }

  // ------ Preprocess-Calculate ------
  recTab.io.en  := (selPCGapReg === fncfg.signal(Reciprocal)) ||
                   (selPCGapReg === fncfg.signal(ATan2Phase1))
  recTab.io.adr := ShiftRegister(recPre.io.adr, pcGap)

  polynomialCoefs(Reciprocal) := recTab.io.cs.asUInt

  recOther.io.x := xdecPCGapReg

  recPost.io.en     := selCPGapReg === fncfg.signal(Reciprocal)
  recPost.io.zother := ShiftRegister(recOther.io.zother, cpGap)
  recPost.io.zres   := polynomialResultCPGapReg

  zs(Reciprocal) := recPost.io.z

  when(selPCReg =/= fncfg.signal(Reciprocal) && selPCReg =/= fncfg.signal(ATan2Phase1)) {
    assert(recPre.io.adr === 0.U)
    assert(recPre.io.dx.getOrElse(0.U) === 0.U)
  }
  when(selPCGapReg =/= fncfg.signal(Reciprocal) && selPCGapReg =/= fncfg.signal(ATan2Phase1)) {
    assert(recTab.io.cs.asUInt === 0.U)
  }

  // ==========================================================================
  // atan2

  val atan2Stage1Pre   = Module(new ATan2Stage1PreProcess (spec, polySpec, stage.preStage))
  val atan2Stage1Other = Module(new ATan2Stage1OtherPath  (spec, polySpec, otherStage))
  val atan2Stage1Post  = Module(new ATan2Stage1PostProcess(spec, polySpec, stage.postStage))

  // atan2Stage1Pre checks if x and y are special values.
  // for calculation, reciprocal is re-used.
  atan2Stage1Pre.io.en := (io.sel === fncfg.signal(ATan2Phase1))
  atan2Stage1Pre.io.x  := xdecomp.io.decomp
  atan2Stage1Pre.io.y  := ydecomp.io.decomp
  atan2Stage1Pre.io.yIsLarger := yIsLarger

  if(order != 0) {
    polynomialDxs.get(ATan2Phase1) := recPre.io.dx.get
  }

  // ------ Preprocess-Calculate ------
  atan2Stage1Other.io.x  := xdecPCGapReg
  atan2Stage1Other.io.y  := ydecPCGapReg
  atan2Stage1Other.io.yIsLarger := yIsLargerPCGapReg

  polynomialCoefs(ATan2Phase1) := recTab.io.cs.asUInt

  atan2Stage1Post.io.en     := (selCPGapReg === fncfg.signal(ATan2Phase1))
  atan2Stage1Post.io.zother := ShiftRegister(atan2Stage1Other.io.zother, cpGap)
  atan2Stage1Post.io.zres   := polynomialResultCPGapReg
  atan2Stage1Post.io.minxy  := Mux(yIsLargerCPGapReg, xdecCPGapReg, ydecCPGapReg)

  zs(ATan2Phase1) := atan2Stage1Post.io.z

  val atan2Stage2Pre   = Module(new ATan2Stage2PreProcess (spec, polySpec, stage.preStage))
  val atan2Stage2Tab   = Module(new ATan2Stage2TableCoeff (spec, polySpec, maxCbit))
  val atan2Stage2Other = Module(new ATan2Stage2OtherPath  (spec, polySpec, otherStage))
  val atan2Stage2Post  = Module(new ATan2Stage2PostProcess(spec, polySpec, stage.postStage))
  atan2Stage2Pre.io.en  := (io.sel === fncfg.signal(ATan2Phase2))
  atan2Stage2Pre.io.x   := xdecomp.io.decomp

  if(order != 0) {
    polynomialDxs.get(ATan2Phase2) := atan2Stage2Pre.io.dx.get
  }

  // ------ Preprocess-Calculate ------
  atan2Stage2Tab.io.en  := (selPCGapReg === fncfg.signal(ATan2Phase2))
  atan2Stage2Tab.io.adr := ShiftRegister(atan2Stage2Pre.io.adr, pcGap)
  atan2Stage2Other.io.x := xdecPCGapReg

  polynomialCoefs(ATan2Phase2) := atan2Stage2Tab.io.cs.asUInt

  atan2Stage2Post.io.en     := (selCPGapReg === fncfg.signal(ATan2Phase2))
  atan2Stage2Post.io.zother := ShiftRegister(atan2Stage2Other.io.zother, cpGap)
  atan2Stage2Post.io.zres   := polynomialResultCPGapReg
  atan2Stage2Post.io.x      := xdecCPGapReg

  zs(ATan2Phase2) := atan2Stage2Post.io.z

  when(selPCReg =/= fncfg.signal(ATan2Phase2)) {
    assert(atan2Stage2Pre.io.adr === 0.U)
    assert(atan2Stage2Pre.io.dx.getOrElse(0.U) === 0.U)
  }
  when(selPCGapReg =/= fncfg.signal(ATan2Phase2)) {
    assert(atan2Stage2Tab.io.cs.asUInt === 0.U)
  }

  // ------------------------------------------------------------------------
  // atan related status register.
  // atan2 stage1 must save some values until atan2 stage2. So, after the pre-
  // process, it saves some flags to register.

  val atan2FlagReg = Reg(new ATan2Flags())
  // the timing is at the cycle when atan2Stage1Pre completes
  when(selPCReg === fncfg.signal(ATan2Phase1)) {
    atan2FlagReg.status  := Cat(yIsLargerPCReg, xdecPCReg.sgn)
    atan2FlagReg.special := atan2Stage1Pre.io.special
    atan2FlagReg.ysgn    := ydecPCReg.sgn
  }
  // The register is updated only when atan2stage1 is executed. That means that
  // we don't need to care about the timing here. ATan2Stage2 can only be
  // executed after the Stage1 because it uses the result of stage1 as its input.
  atan2Stage2Other.io.flags := atan2FlagReg
  atan2Stage2Post.io.flags  := atan2FlagReg

  // ==========================================================================
  // sincos

  if(fncfg.has(Sin) || fncfg.has(Cos)) {
    val sincosPre   = Module(new SinCosPreProcess (spec, polySpec, stage.preStage))
    val sincosTab   = Module(new SinCosTableCoeff (spec, polySpec, maxCbit))
    val sincosPost  = Module(new SinCosPostProcess(spec, polySpec, stage.postStage))

    if(fncfg.has(Sin) && fncfg.has(Cos)) {
      sincosPre.io.en    := (io.sel === fncfg.signal(Sin)) || (io.sel === fncfg.signal(Cos))
      sincosPre.io.isSin := (io.sel === fncfg.signal(Sin))
      sincosPre.io.x     := xdecomp.io.decomp

      sincosTab.io.en  := (selPCGapReg === fncfg.signal(Sin)) ||
                          (selPCGapReg === fncfg.signal(Cos))
      sincosTab.io.adr := ShiftRegister(sincosPre.io.adr, pcGap)

      sincosPost.io.en   := (selCPGapReg === fncfg.signal(Sin)) ||
                            (selCPGapReg === fncfg.signal(Cos))
      sincosPost.io.pre  := ShiftRegister(sincosPre.io.out, pcGap + tcGap + nCalcStage + cpGap)
      sincosPost.io.zres := polynomialResultCPGapReg

      if(order != 0) {
        polynomialDxs.get(Sin) := sincosPre.io.dx.get
        polynomialDxs.get(Cos) := sincosPre.io.dx.get
      }
      polynomialCoefs(Sin) := sincosTab.io.cs.asUInt
      polynomialCoefs(Cos) := sincosTab.io.cs.asUInt
      zs(Sin) := sincosPost.io.z
      zs(Cos) := sincosPost.io.z

      when(selPCReg =/= fncfg.signal(Sin) && selPCReg =/= fncfg.signal(Cos)) {
        assert(sincosPre.io.adr === 0.U)
        assert(sincosPre.io.dx.getOrElse(0.U) === 0.U)
      }
      when(selPCGapReg =/= fncfg.signal(Sin) && selPCGapReg =/= fncfg.signal(Cos)) {
        assert(sincosTab.io.cs.asUInt === 0.U)
      }

    } else if (fncfg.has(Sin)) {

      sincosPre.io.en    := io.sel === fncfg.signal(Sin)
      sincosPre.io.isSin := true.B
      sincosPre.io.x     := xdecomp.io.decomp

      sincosTab.io.en  := selPCGapReg === fncfg.signal(Sin)
      sincosTab.io.adr := ShiftRegister(sincosPre.io.adr, pcGap)

      sincosPost.io.en   := selCPGapReg === fncfg.signal(Sin)
      sincosPost.io.pre  := ShiftRegister(sincosPre.io.out, pcGap + tcGap + nCalcStage + cpGap)
      sincosPost.io.zres := polynomialResultCPGapReg

      if(order != 0) {
        polynomialDxs.get(Sin) := sincosPre.io.dx.get
      }
      polynomialCoefs(Sin) := sincosTab.io.cs.asUInt
      zs(Sin) := sincosPost.io.z

      when(selPCReg =/= fncfg.signal(Sin)) {
        assert(sincosPre.io.adr === 0.U)
        assert(sincosPre.io.dx.getOrElse(0.U) === 0.U)
      }
      when(selPCGapReg =/= fncfg.signal(Sin)) {
        assert(sincosTab.io.cs.asUInt === 0.U)
      }

    } else {

      sincosPre.io.en    := io.sel === fncfg.signal(Cos)
      sincosPre.io.isSin := false.B
      sincosPre.io.x     := xdecomp.io.decomp

      sincosTab.io.en  := selPCGapReg === fncfg.signal(Cos)
      sincosTab.io.adr := ShiftRegister(sincosPre.io.adr, pcGap)

      sincosPost.io.en   := selCPGapReg === fncfg.signal(Cos)
      sincosPost.io.pre  := ShiftRegister(sincosPre.io.out, pcGap + tcGap + nCalcStage + cpGap)
      sincosPost.io.zres := polynomialResultCPGapReg

      if(order != 0) {
        polynomialDxs.get(Cos) := sincosPre.io.dx.get
      }
      polynomialCoefs(Cos) := sincosTab.io.cs.asUInt
      zs(Cos) := sincosPost.io.z

      when(selPCReg =/= fncfg.signal(Cos)) {
        assert(sincosPre.io.adr === 0.U)
        assert(sincosPre.io.dx.getOrElse(0.U) === 0.U)
      }
      when(selPCGapReg =/= fncfg.signal(Cos)) {
        assert(sincosTab.io.cs.asUInt === 0.U)
      }
    }
  }

  // ==========================================================================
  // exp

  if(fncfg.has(Exp)) {
    val expPre   = Module(new ExpPreProcess (spec, polySpec, stage.preStage))
    val expTab   = Module(new ExpTableCoeff (spec, polySpec, maxCbit))
    val expOther = Module(new ExpOtherPath  (spec, polySpec, otherStage))
    val expPost  = Module(new ExpPostProcess(spec, polySpec, stage.postStage))

    expPre.io.en     := (io.sel === fncfg.signal(Exp))
    expPre.io.x      := xdecomp.io.decomp

    if(order != 0) {
      polynomialDxs.get(Exp) := expPre.io.dx.get
    }

    // ------ Preprocess-Calculate ------
    expTab.io.en     := (selPCGapReg === fncfg.signal(Exp))
    expTab.io.adr    := ShiftRegister(expPre.io.adr, pcGap)

    polynomialCoefs(Exp) := expTab.io.cs.asUInt

    expOther.io.x    := xdecPCGapReg
    expOther.io.xint := ShiftRegister(expPre.io.xint, pcGap)
    expOther.io.xexd := ShiftRegister(expPre.io.xexd, pcGap)

    if(expPre.io.xfracLSBs.isDefined) {
      expOther.io.xfracLSBs.get := ShiftRegister(expPre.io.xfracLSBs.get, pcGap)
      expPost.io.zCorrCoef.get  := ShiftRegister(expOther.io.zCorrCoef.get, cpGap)
    }
    expPost.io.en     := selCPGapReg === fncfg.signal(Exp)
    expPost.io.zother := ShiftRegister(expOther.io.zother, cpGap)
    expPost.io.zres   := polynomialResultCPGapReg

    zs(Exp) := expPost.io.z

    when(selPCReg =/= fncfg.signal(Exp)) {
      assert(expPre.io.adr === 0.U)
      assert(expPre.io.dx.getOrElse(0.U) === 0.U)
    }

    when(selPCGapReg =/= fncfg.signal(Exp)) {
      assert(expTab.io.cs.asUInt === 0.U)
    }
  }

  // ==========================================================================
  // log2/ln

  if(fncfg.has(Log)) {
    val logPre   = Module(new LogPreProcess (spec, polySpec, stage.preStage))
    val logTab   = Module(new LogTableCoeff (spec, polySpec, maxCbit))
    val logOther = Module(new LogOtherPath  (spec, polySpec, otherStage))
    val logPost  = Module(new LogPostProcess(spec, polySpec, stage.postStage))

    logPre.io.en      := (io.sel === fncfg.signal(Log))
    logPre.io.x       := xdecomp.io.decomp
    if(order != 0) {
      polynomialDxs.get(Log) := logPre.io.dx.get
    }
    // ------ Preprocess-Calculate ------
    val logPreExAdr = logPre.io.adr(logPre.io.adr.getWidth-1, logPre.io.adr.getWidth-2)
    logTab.io.en      := (selPCGapReg === fncfg.signal(Log))
    logTab.io.adr     := ShiftRegister(logPre.io.adr, pcGap)

    polynomialCoefs(Log) := logTab.io.cs.asUInt

    logOther.io.x     := xdecPCGapReg

    logPost.io.en     := selCPGapReg === fncfg.signal(Log)
    logPost.io.zother := ShiftRegister(logOther.io.zother, cpGap)
    logPost.io.zres   := polynomialResultCPGapReg

    zs(Log) := logPost.io.z

    when(selPCReg =/= fncfg.signal(Log)) {
      assert(logPre.io.adr === 0.U)
      assert(logPre.io.dx.getOrElse(0.U) === 0.U)
    }
    when(selPCGapReg =/= fncfg.signal(Log)) {
      assert(logTab.io.cs.asUInt === 0.U)
    }
  }
}

class MathFuncUnitFP32( stage : MathFuncPipelineConfig )
    extends MathFunctions(MathFuncConfig.all, RealSpec.Float32Spec, 2, 8, 3, stage) {
}

object MathFuncUnitFP32_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new MathFuncUnitFP32(MathFuncPipelineConfig.none)) ) )
}

class MathFuncUnitBF16( stage : MathFuncPipelineConfig )
    extends MathFunctions(MathFuncConfig.all, RealSpec.BFloat16Spec, 0, 7, 1, stage) {
}

object MathFuncUnitBF16_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new MathFuncUnitBF16(MathFuncPipelineConfig.none)) ) )
}

