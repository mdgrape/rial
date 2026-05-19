//% @file MathFunctions.scala
//
// Math functions module (simplified for sqrt/invsqrt/reciprocal)
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import chisel3._
import chisel3.util._

import rial.approx._
import rial.arith.RealSpec
import rial.util.PipelineStageConfig

import scala.collection.mutable.ArrayBuffer
import scala.reflect.ClassTag

/** A module that calculates math functions
 *
 * @constructor Constructs a new MathFunctions module.
 * @param cfg       Configuration that determines which functions are supported.
 * @param spec      RealSpec for input/output floating point numbers.
 * @param preStage  Pipeline stages for preprocess.
 * @param calcStage Pipeline stages for calculation (table lookup + polynomial eval).
 * @param postStage Pipeline stages for postprocess.
 */
class MathFunctions(
  val cfg:  MathFuncConfig,
  val spec: RealSpec,
  val stage: MathFuncPipelineConfig
) extends Module {

  import FuncKind._

  assert(!spec.disableSign, "assuming the float type allows negative values")

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val bitW   = cfg.bitW
  val coeffs = cfg.coeffs

  val io = IO(new Bundle {
    val sel = Input (UInt(cfg.signalW.W))
    val x   = Input (UInt(spec.W.W))
    val z   = Output(UInt(spec.W.W))
  })

  val selIn   = ShiftRegister(io.sel,                 1, 0.U, true.B)
  val selTab  = ShiftRegister(selIn,   stage.preStage  , 0.U, true.B)
  val selCalc = ShiftRegister(selTab,  stage.tableStage, 0.U, true.B)
  val selPost = ShiftRegister(selCalc, stage.calcStage , 0.U, true.B)
  val selOut  = ShiftRegister(selPost, stage.postStage , 0.U, true.B)

  // ==========================================================================
  // Output selection
  // ==========================================================================

  val zs = cfg.funcs.map(fn => {
    fn -> WireDefault(0.U.asTypeOf(Valid(UInt(spec.W.W))))
  }).toMap

  io.z := MuxCase(0.U, cfg.funcs.map(fn => {
    (selOut === cfg.signal(fn) && zs(fn).valid) -> zs(fn).bits
  }))

  // ==========================================================================
  // Decompose input and create delayed versions
  // ==========================================================================

  val xdecomp = Module(new DecomposeReal(spec))
  xdecomp.io.real.valid := cfg.isValid(io.sel)
  xdecomp.io.real.bits  := io.x

  val xIn   = Pipe(xdecomp.io.decomp, 1)
  val xTab  = Pipe(xIn,   stage.preStage  )
  val xCalc = Pipe(xTab,  stage.tableStage)
  val xPost = Pipe(xCalc, stage.calcStage )
  val xOut  = Pipe(xPost, stage.postStage )

  // ==========================================================================
  // Polynomial Evaluator (shared)
  // ==========================================================================

  val polynomialEval = Module(new PolynomialEval(bitW, stage.calcStage))

  polynomialEval.io.in.valid := false.B

  // These will be OR-ed together from different function modules
  val polynomialCoefW = bitW.coefBits.reduce(_+_)
  val polynomialCoefs = cfg.funcs.map(fn => {
    fn -> WireDefault(0.U.asTypeOf(Valid(UInt(polynomialCoefW.W))))
  }).toMap

  val polynomialDxs = if (bitW.nOrder == 0) { None } else { Some(
    cfg.funcs.map(fn => {
      fn -> WireDefault(0.U.asTypeOf(Valid(UInt(bitW.inBitW.W))))
    }).toMap
  )}

  val polynomialCoef = MuxCase(0.U.asTypeOf(Valid(UInt(polynomialCoefW.W))),
    cfg.funcs.map(fn => {
      (selCalc === cfg.signal(fn)) -> polynomialCoefs(fn)
    })
  )
  val polynomialDx = if (bitW.nOrder == 0) { None } else { Some(
    MuxCase(0.U.asTypeOf(Valid(UInt(bitW.inBitW.W))),
      cfg.funcs.map(fn => {
        (selCalc === cfg.signal(fn)) -> polynomialDxs.get(fn)
      })
    )
  )}

  if(bitW.nOrder == 0) { // w/o dx, table only
    polynomialEval.io.in.valid       := polynomialCoef.valid
    polynomialEval.io.in.bits.coeffs := polynomialCoef.bits.asTypeOf(
      MixedVec(bitW.coefBits.map(w => SInt(w.W)))
    )
  } else { // w/ dx

    polynomialEval.io.in.valid       := polynomialCoef.valid && polynomialDx.get.valid
    polynomialEval.io.in.bits.coeffs := polynomialCoef.bits.asTypeOf(
      MixedVec(bitW.coefBits.map(w => SInt(w.W)))
    )
    polynomialEval.io.in.bits.x.get  := polynomialDx.get.bits.asSInt
  }

  // ==========================================================================
  // actual coeffs and pre/post process
  // ==========================================================================

  val pres: Map[FuncKind, Module]  = cfg.funcs.map {
    case Sqrt       => {Sqrt       -> Module(new SqrtPreProcess      (spec, bitW, stage.preStage))}
    case InvSqrt    => {InvSqrt    -> Module(new InvSqrtPreProcess   (spec, bitW, stage.preStage))}
    case Reciprocal => {Reciprocal -> Module(new ReciprocalPreProcess(spec, bitW, stage.preStage))}
  }.toMap[FuncKind, Module]

  val tabs: Map[FuncKind, Module] = cfg.funcs.map {
    case Sqrt       => {Sqrt       -> Module(new SqrtTableCoeff      (spec, coeffs(Sqrt      ), stage.tableStage))}
    case InvSqrt    => {InvSqrt    -> Module(new InvSqrtTableCoeff   (spec, coeffs(InvSqrt   ), stage.tableStage))}
    case Reciprocal => {Reciprocal -> Module(new ReciprocalTableCoeff(spec, coeffs(Reciprocal), stage.tableStage))}
  }.toMap

  val nontabs: Map[FuncKind, Module] = cfg.funcs.map {
    case Sqrt       => {Sqrt       -> Module(new SqrtNonTable      (spec, bitW, stage.nonTableStage))}
    case InvSqrt    => {InvSqrt    -> Module(new InvSqrtNonTable   (spec, bitW, stage.nonTableStage))}
    case Reciprocal => {Reciprocal -> Module(new ReciprocalNonTable(spec, bitW, stage.nonTableStage))}
  }.toMap

  val posts: Map[FuncKind, Module] = cfg.funcs.map {
    case Sqrt       => {Sqrt       -> Module(new SqrtPostProcess      (spec, bitW, stage.postStage))}
    case InvSqrt    => {InvSqrt    -> Module(new InvSqrtPostProcess   (spec, bitW, stage.postStage))}
    case Reciprocal => {Reciprocal -> Module(new ReciprocalPostProcess(spec, bitW, stage.postStage))}
  }.toMap

  private def toModule[T <: Module](x: Module)(implicit ct: ClassTag[T]): T = {
    require(ct.runtimeClass.isInstance(x),
      s"expected ${ct.runtimeClass.getSimpleName} but got ${x.getClass.getSimpleName}")
    x.asInstanceOf[T]
  }

  // --------------------------------------------------------------------------

  if(cfg.has(Sqrt)) {
    val pre  = toModule[SqrtPreProcess ](pres   (Sqrt))
    val tab  = toModule[SqrtTableCoeff ](tabs   (Sqrt))
    val ntab = toModule[SqrtNonTable   ](nontabs(Sqrt))
    val post = toModule[SqrtPostProcess](posts  (Sqrt))
    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(Sqrt) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Sqrt) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    ntab.io.x := xTab

    post.io.in.valid     := ntab.io.z.valid && polynomialEval.io.out.valid
    post.io.in.bits.info := ntab.io.z.bits
    post.io.in.bits.zres := polynomialEval.io.out.bits

    zs(Sqrt) := post.io.z
  }

  // --------------------------------------------------------------------------

  if(cfg.has(InvSqrt)) {
    val pre  = toModule[InvSqrtPreProcess ](pres   (InvSqrt))
    val tab  = toModule[InvSqrtTableCoeff ](tabs   (InvSqrt))
    val ntab = toModule[InvSqrtNonTable   ](nontabs(InvSqrt))
    val post = toModule[InvSqrtPostProcess](posts  (InvSqrt))
    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(InvSqrt) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(InvSqrt) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    ntab.io.x := xTab

    post.io.in.valid     := ntab.io.z.valid && polynomialEval.io.out.valid
    post.io.in.bits.info := ntab.io.z.bits
    post.io.in.bits.zres := polynomialEval.io.out.bits

    zs(InvSqrt) := post.io.z
  }

  // --------------------------------------------------------------------------

  if(cfg.has(Reciprocal)) {
    val pre  = toModule[ReciprocalPreProcess ](pres   (Reciprocal))
    val tab  = toModule[ReciprocalTableCoeff ](tabs   (Reciprocal))
    val ntab = toModule[ReciprocalNonTable   ](nontabs(Reciprocal))
    val post = toModule[ReciprocalPostProcess](posts  (Reciprocal))
    pre.io.x   := xIn
    tab.io.adr := pre.io.adr

    polynomialCoefs(Reciprocal) := tab.io.cs
    if(bitW.nOrder != 0) {
      polynomialDxs.get(Reciprocal) := Pipe(pre.io.dx.get, stage.tableStage)
    }

    ntab.io.x := xTab

    post.io.in.valid     := ntab.io.z.valid && polynomialEval.io.out.valid
    post.io.in.bits.info := ntab.io.z.bits
    post.io.in.bits.zres := polynomialEval.io.out.bits

    zs(Reciprocal) := post.io.z
  }
}
