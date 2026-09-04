//% @file MultFPIntGeneric.scala
//
// multiply FP and int and return FP
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.arith

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._

/** Generic Floating-Point and Integer Multiplier.
 *
 * It takes `x` and `y` and returns `z` as `z = x * y`, where `y` is an Int.
 *
 * All the input/output floating-point spec are configurable through [[rial.arith.RealSpec]].
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"unpack"` after the operand unpacking
 *  - `"mul"`    after the mantissa product
 *  - `"norm"`   after the normalisation
 *  - `"round"`  after the rounding
 *  - `"out"`    the module output
 *
 * [[getStage]] gives the resulting latency, and `out.valid` is `in.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // three registers, placed by the module
 * new MultFPIntGeneric(spec, 32, true, spec, RoundSpec.roundToEven, 3)
 * // hand-picked positions
 * new MultFPIntGeneric(spec, 32, true, spec, RoundSpec.roundToEven,
 *   PipelineStageConfig.specify("mul" -> 1, "norm" -> 1, "round" -> 1))
 * // all the registers at the module output
 * new MultFPIntGeneric(spec, 32, true, spec, RoundSpec.roundToEven,
 *   PipelineStageConfig.atOut(3))
 * // combinational
 * new MultFPIntGeneric(spec, 32, true, spec, RoundSpec.roundToEven,
 *   PipelineStageConfig.none)
 * }}}
 *
 * {{{
 * class MultFPIntGeneric(...) extends Module {
 *   val io = IO(new Bundle{
 *     val in  = Flipped(Valid(new MultFPIntGeneric.Input(xSpec, yWidth)))
 *     val out = Valid(UInt(zSpec.W.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * @constructor create a chisel Module MultFPIntGeneric
 * @param xSpec          floating point spec of the input
 * @param yWidth         width of the integer input
 * @param ySigned        true if the input integer is signed.
 * @param zSpec          floating point spec of the output
 * @param roundSpec      rounding algorithm
 * @param stage          pipeline registers of this module.
 *
 * For parameters, see also:
 * @see [[rial.arith.RealSpec]]
 * @see [[rial.arith.RoundSpec]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 *
 */
class MultFPIntGeneric(
  xSpec : RealSpec, yWidth : Int, ySigned : Boolean, zSpec : RealSpec,
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = MultFPIntGeneric.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getParam = { (xSpec, yWidth, ySigned, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new MultFPIntGeneric.Input(xSpec, yWidth)))
    val out = Valid(UInt(zSpec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val xIn = cutIn(io.in.bits.x)
  val yIn = cutIn(io.in.bits.y)

  //----------------------------------------------------------------------
  // S1 : decompose / classify and normalisation of the integer operand
  val (xsgn,  xex,  xman) = FloatChiselUtil.decomposeWithLeading1(xSpec, xIn)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, xIn)

//   printf("x = %b|%b(%d)|%b\n", xsgn, xex, xex-xSpec.exBias.U, xman)

  val ysgn = Wire(UInt(1.W))
  val yval = Wire(UInt(yWidth.W))

  if (ySigned) {
    ysgn := yIn(yWidth-1)
    yval := Mux(ysgn === 0.U(1.W), yIn, ~yIn + 1.U) // abs
  } else {
    ysgn := 0.U(1.W)
    yval := yIn
  }

  val yLength = (yWidth-1).U - PriorityEncoder(Reverse(yval))
  val yman    = (yval << ((yWidth-1).U - yLength))(yWidth, 0)
  val yex     = yLength
//   printf("y = %b(%d)\n", io.in.bits.y, io.in.bits.y)
//   printf("yWidth  = %d\n", yWidth.U)
//   printf("yLength = %d\n", yLength)
//   printf("y = %b|%b|%b\n", ysgn, yex, yman)

  val zzero = xzero.asBool || yval === 0.U
  val zinf  = xinf.asBool
  val znan  = xnan.asBool
  val zsgn0 = (xsgn ^ ysgn)

//   printf("zzero, zinf, znan, zsgn0 = %b,%b,%b,%b\n",zzero, zinf, znan, zsgn0)

  //--------------------------- stage "unpack" ---------------------------
  private val cutUnpack = new SharedPipe("unpack", cutIn.valid, stages.at("unpack"))

  val xmanUnpack  = cutUnpack(xman)
  val ymanUnpack  = cutUnpack(yman)
  val xexUnpack   = cutUnpack(xex)
  val yexUnpack   = cutUnpack(yex)
  val zsgn0Unpack = cutUnpack(zsgn0)
  val zzeroUnpack = cutUnpack(zzero)
  val zinfUnpack  = cutUnpack(zinf)
  val znanUnpack  = cutUnpack(znan)

  //----------------------------------------------------------------------
  // S2 : Mantissa product

  val prod      = xmanUnpack * ymanUnpack

//   printf("xman * yman = %b\n", prod)

  //---------------------------- stage "mul" -----------------------------
  private val cutMul = new SharedPipe("mul", cutUnpack.valid, stages.at("mul"))

  val prodMul  = cutMul(prod)
  val xexMul   = cutMul(xexUnpack)
  val yexMul   = cutMul(yexUnpack)
  val zsgn0Mul = cutMul(zsgn0Unpack)
  val zzeroMul = cutMul(zzeroUnpack)
  val zinfMul  = cutMul(zinfUnpack)
  val znanMul  = cutMul(znanUnpack)

  //----------------------------------------------------------------------
  // S3 : normalisation shift and round / sticky selection; in parallel, the
  //      exponent sum (see the exponent section below)

  val bp        = xSpec.manW+yWidth-1
  val roundBits = bp-zSpec.manW
  val moreThan2 = prodMul(bp+1)
  val sticky    = prodMul(roundBits-2,0).orR | (moreThan2 & prodMul(roundBits-1))
  val round     = Mux(moreThan2, prodMul(roundBits), prodMul(roundBits-1))
  val prodShift = Mux(moreThan2, prodMul(bp, roundBits+1), prodMul(bp-1, roundBits))
  val lsb       = prodShift(0)
  // RoundSpec : truncate/ceil always towards 0/infinite
  val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky)

  //---------------------------- stage "norm" ----------------------------
  private val cutNorm = new SharedPipe("norm", cutMul.valid, stages.at("norm"))

  val prodShiftNorm = cutNorm(prodShift)
  val incNorm       = cutNorm(inc)
  val moreThan2Norm = cutNorm(moreThan2)
  val zsgn0Norm     = cutNorm(zsgn0Mul)
  val zzeroNorm     = cutNorm(zzeroMul)
  val zinfNorm      = cutNorm(zinfMul)
  val znanNorm      = cutNorm(znanMul)

  //----------------------------------------------------------------------
  // S4 : rounding increment and exponent adjustment
  val prodRound = prodShiftNorm +& incNorm
  val moreThan2AfterRound = prodRound(zSpec.manW)
  val resMan    = prodRound(zSpec.manW-1,0)

  //----------------------------------------------------------------------
  // Exponent

  val maxEx = ( maskI(xSpec.exW)-1-xSpec.exBias + yWidth - 1 + zSpec.exBias )
  val minEx = 1-xSpec.exBias+zSpec.exBias

  val exW0 = max( log2Up(maxEx+1), log2Up(abs(minEx)+1))
  val exW  = if (minEx<0) {exW0+1} else {exW0}

  val exAdd = (-xSpec.exBias + zSpec.exBias).S(exW.W).asUInt
  // the exponent sum runs in parallel with the normalisation / round selection
  val ex0   = xexMul.pad(exW) + yexMul.pad(exW) + exAdd
  val ex0Norm = cutNorm(ex0)

  val exInc = ex0Norm + (moreThan2Norm | moreThan2AfterRound)

  val exZero = !exInc.orR.asBool // needed?
  val exNeg  = (minEx<0).B && (exInc(exW-1)=/=0.U)
  val exZN   = exZero || exNeg || zzeroNorm
  val exInf  = zinfNorm || exInc(zSpec.exW-1,0).andR | exInc(exW-1, zSpec.exW).orR

  //--------------------------- stage "round" ----------------------------
  private val cutRound = new SharedPipe("round", cutNorm.valid, stages.at("round"))

  val resManRound  = cutRound(resMan)
  val zexBitsRound = cutRound(exInc(zSpec.exW-1,0))
  val exZNRound    = cutRound(exZN)
  val exInfRound   = cutRound(exInf)
  val znanRound    = cutRound(znanNorm)
  val zsgn0Round   = cutRound(zsgn0Norm)

  //----------------------------------------------------------------------
  // S5 : Final (special case muxes)
  val zex = Mux(exZNRound,                0.U(zSpec.exW.W),
            Mux(exInfRound | znanRound,   Fill(zSpec.exW,1.U(1.W)),
                                          zexBitsRound))

  // Final mantissa
  val zman = Mux(exZNRound || exInfRound || znanRound,
    Cat(znanRound, 0.U((zSpec.manW-1).W)), resManRound)
  val zsgn = ( !znanRound ) && zsgn0Round.asBool
  val z0 = if (zSpec.disableSign) {Cat(zex, zman)} else {Cat(zsgn, zex, zman)}

//   printf("z = %b|%b(%d)|%b\n", zsgn, zex, zex-zSpec.exBias.U, zman)

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutRound.valid, stages.at("out"))

  io.out.bits  := cutOut(z0)
  io.out.valid := cutOut.valid
  //printf("x=%x y=%x z=%x\n", io.in.bits.x, io.in.bits.y, io.out.bits)
}

/** Input bundle and pipeline stage points of [[rial.arith.MultFPIntGeneric]].
 *
 * @see [[rial.util.StagePoint]]
 */
object MultFPIntGeneric {

  /** Input of [[rial.arith.MultFPIntGeneric]]: the floating point factor `x`
   *  and the integer factor `y`, read as signed or unsigned according to the
   *  module's `ySigned`.
   *
   * @param xSpec  floating point spec of the input
   * @param yWidth width of the integer input
   */
  class Input(xSpec: RealSpec, yWidth: Int) extends Bundle {
    val x = UInt(xSpec.W.W)
    val y = UInt(yWidth.W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.arith.MultFPIntGeneric]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("unpack", 3), StagePoint("mul",   4), StagePoint("norm", 2),
    StagePoint("round",  3), StagePoint("out",   1))
}

/** Specialization of MultFPIntGeneric for FP64.
 *
 * @see [[rial.arith.MultFPIntGeneric]]
 */
class MultFPInt64( stage : PipelineStageConfig )
    extends MultFPIntGeneric( RealSpec.Float64Spec, 64, true, RealSpec.Float64Spec,
      RoundSpec.roundToEven, stage) {
}
