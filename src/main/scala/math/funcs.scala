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
import rial.arith.RoundSpec
import rial.arith.FloatChiselUtil
import rial.math._


/* Enumerator to represent which function is used in [[rial.math.MathFunctions]]. */
object FuncKind extends Enumeration {
  type FuncKind = Value
  val Sqrt, InvSqrt, Reciprocal, Sin, Cos, ACosPhase1, ACosPhase2, ATan2Phase1, ATan2Phase2, Exp, Log = Value

  def getString(fn: FuncKind): String = {
    fn match {
      case Sqrt => "Sqrt"
      case InvSqrt => "InvSqrt"
      case Reciprocal => "Reciprocal"
      case Sin => "Sin"
      case Cos => "Cos"
      case ACosPhase1 => "ACosPhase1"
      case ACosPhase2 => "ACosPhase2"
      case ATan2Phase1 => "ATan2Phase1"
      case ATan2Phase2 => "ATan2Phase2"
      case Exp => "Exp"
      case Log => "Log"
    }
  }
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

    assert(has(Sqrt), """
      | Since acosPhase1 uses sqrt table, we auto-generate sqrt.
      | In principle, we can remove this assumption, but the implementation will
      | be painful and redundant.
      """.stripMargin)
  }
  if(has(ATan2Phase1) || has(ATan2Phase2)) {
    assert(has(ATan2Phase1) && has(ATan2Phase2), "ATan2 requires both phase 1 and 2.")

    assert(has(Reciprocal), """
      | Since atan2Phase1 uses reciprocal table, we auto-generate reciprocal.
      | In principle, we can remove this assumption, but the implementation will
      | be painful and redundant.
      """.stripMargin)
  }

  /** The width of UInt to represent function select signal.
   */
  val signalW = log2Up(1 + funcs.length) // None + [funcs..]

  /** Returns function select signal. It starts from 1.
   *  If the passed function is not supported, fails.
   *
   *  @param fn An enumerator of the function
   *  @return the signal that corresponds to the function
   */
  def signal(fn: FuncKind): UInt = {
    assert(has(fn), f"function `${FuncKind.getString(fn)}` is currently not supported")
    (funcs.indexWhere(_==fn) + 1).U(signalW.W)
  }

  /** Returns function select signal that runs no function.
   *
   *  @return the signal that corresponds to no function.
   */
  def signalNone(): UInt = {
    0.U(signalW.W)
  }

  def getString: String = {
    if(funcs == MathFuncConfig.all.funcs) {
      "[All functions]"
    } else {
      funcs.map(fn => FuncKind.getString(fn)).mkString("[", ", ", "]")
    }
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

/** A module to multiply polynomial result and a coefficient.
 *  (fracW+1) * (manW+1) -> manW
 *
 * Some of the functions require multiplication at postprocess, like:
 * - sincos    approximates sin(x)/x      in polynomial, so z = polynomial * x
 * - acos      approximates acos(1-x^2)/x in polynomial, so z = polynomial * x
 * - atan2-1   approximates 1/max(x,y),   in polynomial, so z = polynomial * min(x,y)
 * - atan2-2   approximates atan(x)/x     in polynomial, so z = polynomial * x
 * - log(x~1)  approximates log(x)/(1-x)  in polynomial, so z = polynomial * (1-x)
 * - log(x!=1) approximates log2(x)       in polynomial, so z = polynomial * ln2
 *
 * To reduce area of polynomial stage, those share one multiplier.
 *
 */
class PostProcMultiplier(
  val realSpec: RealSpec,
  val roundSpec: RoundSpec,
  val polySpec: PolynomialSpec,
  val stage: PipelineStageConfig,
) extends Module {
  assert(polySpec.manW == realSpec.manW)

  val fracW = polySpec.fracW
  val manW  = realSpec.manW
  val nStage = stage.total

  val io = IO(new Bundle {
    val en    = Input(Bool())
    val lhs   = Input(UInt((1+fracW).W))
    val rhs   = Input(UInt((1+manW).W))
    val out   = Output(UInt(manW.W))
    val exInc = Output(UInt(1.W))
  })

//   printf("cir: en        = %b\n", io.en)
//   printf("cir: lhs(frac) = %b\n", io.lhs)
//   printf("cir: rhs(man)  = %b\n", io.rhs)

  val prod      = enable(io.en, io.lhs) * enable(io.en, io.rhs)
  val moreThan2 = prod((manW+1)+(fracW+1)-1)

  val roundBits = fracW
  val sticky    = prod(roundBits-2, 0).orR | (moreThan2 & prod(roundBits-1))
  val round     = Mux(moreThan2, prod(roundBits), prod(roundBits-1))
  val shifted   = Mux(moreThan2, prod((manW+1)+(fracW+1)-2, (fracW+1)  ),
                                 prod((manW+1)+(fracW+1)-3, (fracW+1)-1))
  val lsb       = shifted(0)
  val roundInc  = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky)
  val rounded   = shifted +& roundInc
  val moreThan2AfterRound = rounded(manW)

  io.exInc := ShiftRegister(moreThan2 | moreThan2AfterRound, nStage)
  io.out   := ShiftRegister(rounded(manW-1, 0),              nStage)
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
      case ATan2Phase2 => ATan2Phase2TableCoeff.getCBits(spec, polySpec)
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
      case ATan2Phase2 => ATan2Phase2TableCoeff.getCalcW(spec, polySpec)
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

  fncfg.funcs.foreach(fn => {
    println(f"${FuncKind.getString(fn)}%-11s: " +
      f"cbits = ${getCBit(fn) .map(b => f"${b}%3s").mkString("[", ",", "]")} " +
      f"calcW = ${getCalcW(fn).map(b => f"${b}%3s").mkString("[", ",", "]")}")
  })
  println(f"Maximum    : " +
    f"cbits = ${maxCbit .map(b => f"${b}%3s").mkString("[", ",", "]")} " +
    f"calcW = ${maxCalcW.map(b => f"${b}%3s").mkString("[", ",", "]")}")

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
  // PostProc multiplier

  def usePostProcMultiplier(fn: FuncKind.FuncKind): Boolean = {
//     fn == ACosPhase2 || fn == ATan2Phase1 || fn == ATan2Phase2 || fn == Sin || fn == Cos || fn == Log
    fn == ACosPhase2 || fn == ATan2Phase1 || fn == Log || fn == Sin || fn == Cos
  }

  val hasPostProcMultiplier = fncfg.funcs.exists(fn => usePostProcMultiplier(fn))

  val multStage = PipelineStageConfig.none // TODO
  val postProcMultiplier = if(hasPostProcMultiplier) {
    Some(Module(new PostProcMultiplier(spec, RoundSpec.roundToEven, polySpec, multStage)))
  } else {None}

  // regsiter?
  val postProcMultEn  = fncfg.funcs.filter(fn => usePostProcMultiplier(fn)).map( fn => { fn -> Wire(Bool()) }).toMap
  val postProcMultLhs = fncfg.funcs.filter(fn => usePostProcMultiplier(fn)).map( fn => { fn -> Wire(UInt((1+polySpec.fracW).W)) }).toMap
  val postProcMultRhs = fncfg.funcs.filter(fn => usePostProcMultiplier(fn)).map( fn => { fn -> Wire(UInt((1+manW).W)) }).toMap

  postProcMultEn .values.foreach(v => v := false.B)
  postProcMultLhs.values.foreach(v => v := 0.U)
  postProcMultRhs.values.foreach(v => v := 0.U)

//   postProcMultEn .foreach(kv => printf("cir: %d: en  = %b\n", fncfg.signal(kv._1), kv._2))
//   postProcMultLhs.foreach(kv => printf("cir: %d: lhs = %b\n", fncfg.signal(kv._1), kv._2))
//   postProcMultRhs.foreach(kv => printf("cir: %d: lhs = %b\n", fncfg.signal(kv._1), kv._2))

  if(hasPostProcMultiplier) {
    postProcMultiplier.get.io.en  := postProcMultEn.values.reduce(_|_)
    postProcMultiplier.get.io.lhs := postProcMultLhs.values.reduce(_|_)
    postProcMultiplier.get.io.rhs := postProcMultRhs.values.reduce(_|_)
  }

  // ==========================================================================
  // Output float.
  //
  // later we will insert a value to the map.

  val zs = fncfg.funcs.map( fn => { fn -> Wire(UInt(spec.W.W)) }).toMap
  zs.values.foreach(_ := 0.U) // init
  val z0 = zs.values.reduce(_|_)
  io.z := z0

  // ==========================================================================
  // ACos

  if(fncfg.has(ACosPhase1) || fncfg.has(ACosPhase2)) {
    assert(fncfg.signal(ACosPhase1) =/= 0.U)
    assert(fncfg.signal(ACosPhase2) =/= 0.U)

    val acosFlagReg = Reg(new ACosFlags())

    val acos1Pre  = Module(new ACosPhase1PreProcess(spec, polySpec, stage.preStage))
    val acos2Pre  = Module(new ACosPhase2PreProcess(spec, polySpec, stage.preStage))
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

    assert(hasPostProcMultiplier, "ACos requires post-proc multiplier")

    val isACos2 = selCPGapReg === fncfg.signal(ACosPhase2)
    postProcMultEn(ACosPhase2)  := isACos2
    postProcMultLhs(ACosPhase2) := enable(isACos2, Cat(1.U(1.W), polynomialResultCPGapReg))
    postProcMultRhs(ACosPhase2) := enable(isACos2, Cat(1.U(1.W), xdecCPGapReg.man))

    acosPost.io.en     := (selCPGapReg === fncfg.signal(ACosPhase2))
    acosPost.io.flags  := acosFlagReg
    acosPost.io.zex0   := xdecCPGapReg.ex
    acosPost.io.zman0  := postProcMultiplier.get.io.out   // TODO consider nStage
    acosPost.io.exInc  := postProcMultiplier.get.io.exInc

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

    // --------------------------------------------------------------------------
    // sqrt

    val sqrtPre   = Module(new SqrtPreProcess (spec, polySpec, stage.preStage))
    val sqrtTab   = Module(new SqrtTableCoeff (spec, polySpec, maxCbit))
    val sqrtOther = Module(new SqrtOtherPath  (spec, polySpec, otherStage))
    val sqrtPost  = Module(new SqrtPostProcess(spec, polySpec, stage.postStage))

    val sqrtPreAdrPCGapReg = ShiftRegister(sqrtPre.io.adr, pcGap)

    if(fncfg.has(InvSqrt)) {
      sqrtPre.io.en := io.sel === fncfg.signal(Sqrt) || io.sel === fncfg.signal(InvSqrt)
    } else {
      sqrtPre.io.en := io.sel === fncfg.signal(Sqrt)
    }
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
    if(fncfg.has(InvSqrt)) {
      when(selPCReg =/= fncfg.signal(Sqrt) && selPCReg =/= fncfg.signal(InvSqrt)) {
        assert(sqrtPre.io.adr === 0.U)
        assert(sqrtPre.io.dx.getOrElse(0.U) === 0.U)
      }
    } else {
      when(selPCReg =/= fncfg.signal(Sqrt)) {
        assert(sqrtPre.io.adr === 0.U)
        assert(sqrtPre.io.dx.getOrElse(0.U) === 0.U)
      }
    }
    when(selPCGapReg =/= fncfg.signal(Sqrt)) {
      assert(sqrtTab.io.cs.asUInt === 0.U || selPCGapReg === fncfg.signal(ACosPhase1))
    }

    if(fncfg.has(InvSqrt)) {
      // ----------------------------------------------------------------------
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
    }
  }

  // ==========================================================================
  // sqrt w/o acos

  if( ! fncfg.has(ACosPhase1) && ! fncfg.has(ACosPhase2) && fncfg.has(Sqrt)) {
    // --------------------------------------------------------------------------
    // sqrt

    val sqrtPre   = Module(new SqrtPreProcess (spec, polySpec, stage.preStage))
    val sqrtTab   = Module(new SqrtTableCoeff (spec, polySpec, maxCbit))
    val sqrtOther = Module(new SqrtOtherPath  (spec, polySpec, otherStage))
    val sqrtPost  = Module(new SqrtPostProcess(spec, polySpec, stage.postStage))

    val sqrtPreAdrPCGapReg = ShiftRegister(sqrtPre.io.adr, pcGap)

    if(fncfg.has(InvSqrt)) {
      sqrtPre.io.en  := io.sel === fncfg.signal(Sqrt) || io.sel === fncfg.signal(InvSqrt)
    } else {
      sqrtPre.io.en  := io.sel === fncfg.signal(Sqrt)
    }

    sqrtPre.io.x   := xdecomp.io.decomp
    if(order != 0) {
      polynomialDxs.get(Sqrt) := sqrtPre.io.dx.get
    }

    // ------ Preprocess-Calculate ------

    sqrtTab.io.en  := selPCGapReg === fncfg.signal(Sqrt)
    sqrtTab.io.adr := sqrtPreAdrPCGapReg

    // redundant...
    polynomialCoefs(Sqrt)        := sqrtTab.io.cs.asUInt

    sqrtOther.io.x := xdecPCGapReg

    sqrtPost.io.en     := selCPGapReg === fncfg.signal(Sqrt)
    sqrtPost.io.zother := ShiftRegister(sqrtOther.io.zother, cpGap)
    sqrtPost.io.zres   := polynomialResultCPGapReg

    // redundant...
    zs(Sqrt) := sqrtPost.io.z

    if(fncfg.has(InvSqrt)) {
      when(selPCReg =/= fncfg.signal(Sqrt) && selPCReg =/= fncfg.signal(InvSqrt)) {
        assert(sqrtPre.io.adr === 0.U)
        assert(sqrtPre.io.dx.getOrElse(0.U) === 0.U)
      }
    } else {
      when(selPCReg =/= fncfg.signal(Sqrt)) {
        assert(sqrtPre.io.adr === 0.U)
        assert(sqrtPre.io.dx.getOrElse(0.U) === 0.U)
      }
    }
    when(selPCGapReg =/= fncfg.signal(Sqrt)) {
      assert(sqrtTab.io.cs.asUInt === 0.U)
    }

    if(fncfg.has(InvSqrt)) {
      // --------------------------------------------------------------------------
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
    }
  }

  // invsqrt w/o sqrt
  if( ! fncfg.has(ACosPhase1) && ! fncfg.has(ACosPhase2) && ! fncfg.has(Sqrt) && fncfg.has(InvSqrt)) {
    // --------------------------------------------------------------------------
    // invsqrt
    val sqrtPre      = Module(new SqrtPreProcess (spec, polySpec, stage.preStage))
    val invsqrtTab   = Module(new InvSqrtTableCoeff (spec, polySpec, maxCbit))
    val invsqrtOther = Module(new InvSqrtOtherPath  (spec, polySpec, otherStage))
    val invsqrtPost  = Module(new InvSqrtPostProcess(spec, polySpec, stage.postStage))

    val sqrtPreAdrPCGapReg = ShiftRegister(sqrtPre.io.adr, pcGap)

    sqrtPre.io.en  := io.sel === fncfg.signal(InvSqrt)
    sqrtPre.io.x   := xdecomp.io.decomp

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

    when(selPCReg =/= fncfg.signal(InvSqrt)) {
      assert(sqrtPre.io.adr === 0.U)
      assert(sqrtPre.io.dx.getOrElse(0.U) === 0.U)
    }
    when(selPCGapReg =/= fncfg.signal(InvSqrt)) {
      assert(invsqrtTab.io.cs.asUInt === 0.U)
    }
  }

  // ==========================================================================
  // atan2

  if(fncfg.has(ATan2Phase1) || fncfg.has(ATan2Phase2)) {
    assert(fncfg.has(ATan2Phase1) && fncfg.has(ATan2Phase2) && fncfg.has(Reciprocal))

    // --------------------------------------------------------------------------
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

    // --------------------------------------------------------------------------
    // atan2

    val atan2Phase1Pre   = Module(new ATan2Phase1PreProcess (spec, polySpec, stage.preStage))
    val atan2Phase1Other = Module(new ATan2Phase1OtherPath  (spec, polySpec, otherStage))
    val atan2Phase1Post  = Module(new ATan2Phase1PostProcess(spec, polySpec, stage.postStage))

    // atan2Phase1Pre checks if x and y are special values.
    // for calculation, reciprocal is re-used.
    atan2Phase1Pre.io.en := (io.sel === fncfg.signal(ATan2Phase1))
    atan2Phase1Pre.io.x  := xdecomp.io.decomp
    atan2Phase1Pre.io.y  := ydecomp.io.decomp
    atan2Phase1Pre.io.yIsLarger := yIsLarger

    if(order != 0) {
      polynomialDxs.get(ATan2Phase1) := recPre.io.dx.get
    }

    // ------ Preprocess-Calculate ------
    atan2Phase1Other.io.x  := xdecPCGapReg
    atan2Phase1Other.io.y  := ydecPCGapReg
    atan2Phase1Other.io.yIsLarger := yIsLargerPCGapReg

    polynomialCoefs(ATan2Phase1) := recTab.io.cs.asUInt

    // TODO:
    // consider adding PostProcess-1 and -2, meaning before and after PostProcMult.

    val isATan21 = selCPGapReg === fncfg.signal(ATan2Phase1)
    val minXY    = Mux(yIsLargerCPGapReg, xdecCPGapReg, ydecCPGapReg)
    val zother   = ShiftRegister(atan2Phase1Other.io.zother, cpGap)
    val zFrac    = Mux(zother.maxXYMan0, 0.U, polynomialResultCPGapReg)

    postProcMultEn(ATan2Phase1)  := isATan21
    postProcMultLhs(ATan2Phase1) := enable(isATan21, Cat(1.U(1.W), zFrac))
    postProcMultRhs(ATan2Phase1) := enable(isATan21, Cat(1.U(1.W), minXY.man))

    atan2Phase1Post.io.en     := (selCPGapReg === fncfg.signal(ATan2Phase1))
    atan2Phase1Post.io.zother := zother
    atan2Phase1Post.io.zman0  := postProcMultiplier.get.io.out
    atan2Phase1Post.io.minxy  := minXY

    zs(ATan2Phase1) := atan2Phase1Post.io.z

    val atan2Phase2Pre   = Module(new ATan2Phase2PreProcess (spec, polySpec, stage.preStage))
    val atan2Phase2Tab   = Module(new ATan2Phase2TableCoeff (spec, polySpec, maxCbit))
    val atan2Phase2Other = Module(new ATan2Phase2OtherPath  (spec, polySpec, otherStage))
    val atan2Phase2Post  = Module(new ATan2Phase2PostProcess(spec, polySpec, stage.postStage))
    atan2Phase2Pre.io.en  := (io.sel === fncfg.signal(ATan2Phase2))
    atan2Phase2Pre.io.x   := xdecomp.io.decomp

    if(order != 0) {
      polynomialDxs.get(ATan2Phase2) := atan2Phase2Pre.io.dx.get
    }

    // ------ Preprocess-Calculate ------
    atan2Phase2Tab.io.en  := (selPCGapReg === fncfg.signal(ATan2Phase2))
    atan2Phase2Tab.io.adr := ShiftRegister(atan2Phase2Pre.io.adr, pcGap)
    atan2Phase2Other.io.x := xdecPCGapReg

    polynomialCoefs(ATan2Phase2) := atan2Phase2Tab.io.cs.asUInt

    atan2Phase2Post.io.en     := (selCPGapReg === fncfg.signal(ATan2Phase2))
    atan2Phase2Post.io.zother := ShiftRegister(atan2Phase2Other.io.zother, cpGap)
    atan2Phase2Post.io.zres   := polynomialResultCPGapReg
    atan2Phase2Post.io.x      := xdecCPGapReg

    zs(ATan2Phase2) := atan2Phase2Post.io.z

    when(selPCReg =/= fncfg.signal(ATan2Phase2)) {
      assert(atan2Phase2Pre.io.adr === 0.U)
      assert(atan2Phase2Pre.io.dx.getOrElse(0.U) === 0.U)
    }
    when(selPCGapReg =/= fncfg.signal(ATan2Phase2)) {
      assert(atan2Phase2Tab.io.cs.asUInt === 0.U)
    }

    // ------------------------------------------------------------------------
    // atan related status register.
    // atan2 stage1 must save some values until atan2 stage2. So, after the pre-
    // process, it saves some flags to register.

    val atan2FlagReg = Reg(new ATan2Flags())
    // the timing is at the cycle when atan2Phase1Pre completes
    when(selPCReg === fncfg.signal(ATan2Phase1)) {
      atan2FlagReg.status  := Cat(yIsLargerPCReg, xdecPCReg.sgn)
      atan2FlagReg.special := atan2Phase1Pre.io.special
      atan2FlagReg.ysgn    := ydecPCReg.sgn
    }
    // The register is updated only when atan2stage1 is executed. That means that
    // we don't need to care about the timing here. ATan2Phase2 can only be
    // executed after the Phase1 because it uses the result of stage1 as its input.
    atan2Phase2Other.io.flags := atan2FlagReg
    atan2Phase2Post.io.flags  := atan2FlagReg
  }

  // ==========================================================================
  // reciprocal w/o ATan2

  if( ! fncfg.has(ATan2Phase1) && ! fncfg.has(ATan2Phase2) && fncfg.has(Reciprocal)) {

    val recPre   = Module(new ReciprocalPreProcess (spec, polySpec, stage.preStage))
    val recTab   = Module(new ReciprocalTableCoeff (spec, polySpec, maxCbit))
    val recOther = Module(new ReciprocalOtherPath  (spec, polySpec, otherStage))
    val recPost  = Module(new ReciprocalPostProcess(spec, polySpec, stage.postStage))

    recPre.io.en  := io.sel === fncfg.signal(Reciprocal)
    recPre.io.x   := xdecomp.io.decomp

    if(order != 0) {
      polynomialDxs.get(Reciprocal)  := recPre.io.dx.get
    }

    // ------ Preprocess-Calculate ------
    recTab.io.en  := (selPCGapReg === fncfg.signal(Reciprocal))
    recTab.io.adr := ShiftRegister(recPre.io.adr, pcGap)

    polynomialCoefs(Reciprocal) := recTab.io.cs.asUInt

    recOther.io.x := xdecPCGapReg

    recPost.io.en     := selCPGapReg === fncfg.signal(Reciprocal)
    recPost.io.zother := ShiftRegister(recOther.io.zother, cpGap)
    recPost.io.zres   := polynomialResultCPGapReg

    zs(Reciprocal) := recPost.io.z

    when(selPCReg =/= fncfg.signal(Reciprocal)) {
      assert(recPre.io.adr === 0.U)
      assert(recPre.io.dx.getOrElse(0.U) === 0.U)
    }
    when(selPCGapReg =/= fncfg.signal(Reciprocal)) {
      assert(recTab.io.cs.asUInt === 0.U)
    }
  }

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

      if(order != 0) {
        polynomialDxs.get(Sin) := sincosPre.io.dx.get
        polynomialDxs.get(Cos) := sincosPre.io.dx.get
      }
      polynomialCoefs(Sin) := sincosTab.io.cs.asUInt
      polynomialCoefs(Cos) := sincosTab.io.cs.asUInt

      // TODO consider nStage
      val preOut = ShiftRegister(sincosPre.io.out, pcGap + tcGap + nCalcStage + cpGap)

      assert(hasPostProcMultiplier, "Sin/Cos requires post-proc multiplier")

      val isSin = (selCPGapReg === fncfg.signal(Sin))
      postProcMultEn(Sin)  := isSin
      postProcMultLhs(Sin) := enable(isSin, Cat(1.U(1.W), polynomialResultCPGapReg))
      postProcMultRhs(Sin) := enable(isSin, Cat(1.U(1.W), preOut.yman))

      val isCos = (selCPGapReg === fncfg.signal(Cos))
      postProcMultEn(Cos)  := isCos
      postProcMultLhs(Cos) := enable(isCos, Cat(1.U(1.W), polynomialResultCPGapReg))
      postProcMultRhs(Cos) := enable(isCos, Cat(1.U(1.W), preOut.yman))

      sincosPost.io.en   := (selCPGapReg === fncfg.signal(Sin)) ||
                            (selCPGapReg === fncfg.signal(Cos))
      sincosPost.io.pre  := preOut
      sincosPost.io.zman0  := postProcMultiplier.get.io.out
      sincosPost.io.zexInc := postProcMultiplier.get.io.exInc

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

      if(order != 0) {
        polynomialDxs.get(Sin) := sincosPre.io.dx.get
      }
      polynomialCoefs(Sin) := sincosTab.io.cs.asUInt

     // TODO consider nStage
      val preOut = ShiftRegister(sincosPre.io.out, pcGap + tcGap + nCalcStage + cpGap)

      assert(hasPostProcMultiplier, "Sin/Cos requires post-proc multiplier")

      val isSin = (selCPGapReg === fncfg.signal(Sin))
      postProcMultEn(Sin)  := isSin
      postProcMultLhs(Sin) := enable(isSin, Cat(1.U(1.W), polynomialResultCPGapReg))
      postProcMultRhs(Sin) := enable(isSin, Cat(1.U(1.W), preOut.yman))

      sincosPost.io.en   := (selCPGapReg === fncfg.signal(Sin))
      sincosPost.io.pre  := preOut
      sincosPost.io.zman0  := postProcMultiplier.get.io.out
      sincosPost.io.zexInc := postProcMultiplier.get.io.exInc

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

      if(order != 0) {
        polynomialDxs.get(Cos) := sincosPre.io.dx.get
      }
      polynomialCoefs(Cos) := sincosTab.io.cs.asUInt

      // TODO consider nStage
      val preOut = ShiftRegister(sincosPre.io.out, pcGap + tcGap + nCalcStage + cpGap)

      assert(hasPostProcMultiplier, "Sin/Cos requires post-proc multiplier")

      val isCos = (selCPGapReg === fncfg.signal(Cos))
      postProcMultEn(Cos)  := isCos
      postProcMultLhs(Cos) := enable(isCos, Cat(1.U(1.W), polynomialResultCPGapReg))
      postProcMultRhs(Cos) := enable(isCos, Cat(1.U(1.W), preOut.yman))

      sincosPost.io.en   := (selCPGapReg === fncfg.signal(Cos))
      sincosPost.io.pre  := preOut
      sincosPost.io.zman0  := postProcMultiplier.get.io.out
      sincosPost.io.zexInc := postProcMultiplier.get.io.exInc

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
    logTab.io.en  := (selPCGapReg === fncfg.signal(Log))
    logTab.io.adr := ShiftRegister(logPre.io.adr, pcGap)

    polynomialCoefs(Log) := logTab.io.cs.asUInt

    logOther.io.x := xdecPCGapReg
    val nonTableOut = ShiftRegister(logOther.io.zother, cpGap)

    assert(hasPostProcMultiplier, "Log requires post-proc multiplier")

    // XXX Right place to put the following logic?
    // - [x] put here
    // - [ ] split logPostProc
    // - [ ] add io ports to logPostProc

    val fracW = polySpec.fracW

    val x0_5to1_0  = nonTableOut.x0_5to1_0
    val x1_0to2_0  = nonTableOut.x1_0to2_0
    val xOtherwise = nonTableOut.xOtherwise

    val zsgn = nonTableOut.zsgn

    // case 1: x < 0.5 or 2 <= x. calc (ex + log2(1.man)) * ln2.
    val zfrac0 = Mux(zsgn === 0.U, polynomialResultCPGapReg, // means 0 <= xexNobias
                 Mux(polynomialResultCPGapReg === 0.U,
                   Fill(fracW, 1.U(1.W)),
                   ~polynomialResultCPGapReg + 1.U))

    val zfull0 = Cat(nonTableOut.zint, zfrac0)
    val zfullShifted = (zfull0 << nonTableOut.zintShift)(exW+fracW-1, 0)
    assert(zfullShifted(exW + fracW-1) === 1.U || !xOtherwise)

    // fracW+1 width
    val ymanW1 = Mux(x0_5to1_0, Cat(1.U(1.W), polynomialResultCPGapReg), // case 2
                 Mux(x1_0to2_0, Cat(polynomialResultCPGapReg, 0.U(1.W)), // case 3
                                zfullShifted(exW+fracW-1, exW-1)))       // case 1

    val isLog = selCPGapReg === fncfg.signal(Log)
    postProcMultEn(Log)  := isLog
    postProcMultLhs(Log) := enable(isLog, ymanW1)
    postProcMultRhs(Log) := enable(isLog, Cat(1.U(1.W), nonTableOut.constant))

    logPost.io.en     := selCPGapReg === fncfg.signal(Log)
    logPost.io.zother := nonTableOut
    logPost.io.zman0  := postProcMultiplier.get.io.out
    logPost.io.zexInc := postProcMultiplier.get.io.exInc

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

