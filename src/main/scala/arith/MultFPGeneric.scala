//% @file BasicArithmetics.scala
//
// Basic Arithmetic Units
// Copyright (C) Makoto Taiji RIKEN BDR 2020
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

/** Generic Floating-Point Multiplier.
 *
 * It takes `x` and `y` and returns `z` as `z = x * y`.
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
 * new MultFPGeneric(spec, spec, spec, RoundSpec.roundToEven, 3)
 * // hand-picked positions
 * new MultFPGeneric(spec, spec, spec, RoundSpec.roundToEven,
 *   PipelineStageConfig.specify("mul" -> 1, "norm" -> 1, "round" -> 1))
 * // all the registers at the module output
 * new MultFPGeneric(spec, spec, spec, RoundSpec.roundToEven, PipelineStageConfig.atOut(3))
 * // combinational
 * new MultFPGeneric(spec, spec, spec, RoundSpec.roundToEven, PipelineStageConfig.none)
 * }}}
 *
 * {{{
 * class MultFPGeneric(...) extends Module {
 *   val io = IO(new Bundle{
 *     val in  = Flipped(Valid(new MultFPGeneric.Input(xSpec, ySpec)))
 *     val out = Valid(UInt(zSpec.W.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * @constructor create a chisel Module MultFPGeneric
 * @param xSpec          floating point spec of the left hand side input
 * @param ySpec          floating point spec of the right hand side input
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
class MultFPGeneric(
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = MultFPGeneric.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getParam = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new MultFPGeneric.Input(xSpec, ySpec)))
    val out = Valid(UInt(zSpec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val xIn = cutIn(io.in.bits.x)
  val yIn = cutIn(io.in.bits.y)

  //----------------------------------------------------------------------
  // S1 : decompose / classify
  val (xsgn, xex, xman) = FloatChiselUtil.decomposeWithLeading1(xSpec, xIn)
  val (ysgn, yex, yman) = FloatChiselUtil.decomposeWithLeading1(ySpec, yIn)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, xIn)
  val (yzero, yinf, ynan) = FloatChiselUtil.checkValue(ySpec, yIn)
  val zsgn0 = xsgn ^ ysgn

  val zzero = (xzero | yzero).asBool
  val zinf  = (xinf | yinf).asBool
  val znan  = (xnan | ynan).asBool & !zzero
  //printf("%x %x\n",  io.in.bits.x, io.in.bits.y)
  //printf("%b %b %b\n",  zzero, zinf, znan)

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
  val prod = xmanUnpack * ymanUnpack

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
  val bp = xSpec.manW+ySpec.manW
  val roundBits = bp-zSpec.manW
  val moreThan2 = prodMul(bp+1)

  //---------------------------- stage "norm" ----------------------------
  // The normalisation shift and the rounding increment live in the same Scala
  // block below, so the registers of this point are inserted inside it, right
  // after the shift and the round / sticky selection.
  private val cutNorm = new SharedPipe("norm", cutMul.valid, stages.at("norm"))

  val (moreThan2AfterRound, resMan) =
    if(roundBits >= 0) {
      val stickey = if(roundBits >= 2) {
        prodMul(roundBits-2,0).orR | (moreThan2 & prodMul(roundBits-1))
      } else if(roundBits == 0) {
        moreThan2 & prodMul(roundBits-1)
      } else {
        0.U(1.W).asBool
      }
      val round = if(roundBits >= 1) {
        Mux(moreThan2, prodMul(roundBits), prodMul(roundBits-1))
      } else {
        moreThan2 & prodMul(0)
      }
      val prodShift = Mux(moreThan2, prodMul(bp, roundBits+1), prodMul(bp-1, roundBits))
      val lsb   = prodShift(0)
      // RoundSpec : truncate/ceil always towards 0/infinite
      val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, stickey)
      val prodShiftNorm = cutNorm(prodShift)
      val incNorm       = cutNorm(inc)
      val prodRound = prodShiftNorm +& incNorm
      val moreThan2AfterRound = prodRound(zSpec.manW)
      (moreThan2AfterRound, prodRound(zSpec.manW-1, 0))
    } else {
      val prodShift = Mux(moreThan2, prodMul, Cat(prodMul(prodMul.getWidth-2, 0), 0.U(1.W)))
      assert(prodShift.getWidth == prodMul.getWidth)
      assert(!cutMul.valid || prodShift(prodShift.getWidth-1) === 1.U(1.W))
      assert(-1 >= roundBits)
      val prodRound = Cat(prodShift, Fill(-roundBits-1, 0.U(1.W)))
      assert(!cutMul.valid || prodRound(zSpec.manW) === 1.U(1.W))
      // no rounding increment on this path, so the cut is a plain rewiring
      (0.U(1.W), cutNorm(prodRound(zSpec.manW-1, 0)))
    }

  val moreThan2Norm = cutNorm(moreThan2)
  val zsgn0Norm     = cutNorm(zsgn0Mul)
  val zzeroNorm     = cutNorm(zzeroMul)
  val zinfNorm      = cutNorm(zinfMul)
  val znanNorm      = cutNorm(znanMul)

  //----------------------------------------------------------------------
  // Exponent
  val maxEx = ( xSpec.exMax + ySpec.exMax + 1 + zSpec.exBias )
  //println(f"${xSpec.exBias}%x ${ySpec.exBias}%x ${zSpec.exBias}%x")
  //println(f"${xSpec.exW}%d ${ySpec.exW}%d")
  //println(f"${maskI(xSpec.exW)}%x ${maskI(ySpec.exW)}%x")
  val minEx = xSpec.exMin + ySpec.exMin + zSpec.exBias
  val exW0 = max( log2Up(maxEx+1), log2Up(abs(minEx)+1))
  val exW  = if (minEx<0) {exW0+1} else {exW0} // width for calculation
  //println(f"maxEX=$maxEx%x minEx=-${-minEx}%x exW=$exW%d")
  val exAdd = (-xSpec.exBias-ySpec.exBias+zSpec.exBias).S(exW.W).asUInt
  // the exponent sum runs in parallel with the normalisation / round selection
  val ex0 = xexMul.pad(exW) + yexMul.pad(exW) + exAdd
  val ex0Norm = cutNorm(ex0)

  //----------------------------------------------------------------------
  // S4 : rounding increment (above) and exponent adjustment
  val exInc = ex0Norm + (moreThan2Norm | moreThan2AfterRound)

  val exZero = !exInc.orR.asBool
  val exNeg  = (minEx<0).B && (exInc(exW-1)=/=0.U)
  val exZN   = exZero || exNeg || zzeroNorm
  val exInf  = if(exW > zSpec.exW) {
    zinfNorm || exInc(zSpec.exW-1,0).andR | exInc(exW-1, zSpec.exW).orR
  } else if(exW == zSpec.exW) {
    zinfNorm || exInc(zSpec.exW-1,0).andR
  } else {
    zinfNorm
  }

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
  val zex =
    Mux(exZNRound, 0.U(zSpec.exW.W),
      Mux(exInfRound | znanRound, Fill(zSpec.exW,1.U(1.W)), zexBitsRound) )

  // Final mantissa
  val zman =
    Mux(exZNRound || exInfRound || znanRound, Cat(znanRound, 0.U((zSpec.manW-1).W)), resManRound)
  val zsgn = ( !znanRound ) && zsgn0Round.asBool
  //printf("%x\n", zman)

  assert(zman.getWidth == zSpec.manW)

  val z0 = if (zSpec.disableSign) {Cat(zex, zman)} else {Cat(zsgn, zex, zman)}

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutRound.valid, stages.at("out"))

  io.out.bits  := cutOut(z0)
  io.out.valid := cutOut.valid
  //printf("x=%x y=%x z=%x\n", io.in.bits.x, io.in.bits.y, io.out.bits)
}

/** Input bundle and pipeline stage points of [[rial.arith.MultFPGeneric]].
 *
 * @see [[rial.util.StagePoint]]
 */
object MultFPGeneric {

  /** Input of [[rial.arith.MultFPGeneric]]: the two factors.
   *
   * @param xSpec floating point spec of the left hand side input
   * @param ySpec floating point spec of the right hand side input
   */
  class Input(xSpec: RealSpec, ySpec: RealSpec) extends Bundle {
    val x = UInt(xSpec.W.W)
    val y = UInt(ySpec.W.W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.arith.MultFPGeneric]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("unpack", 2), StagePoint("mul",   4), StagePoint("norm", 2),
    StagePoint("round",  3), StagePoint("out",   1))
}

/** Specialization of MultFPGeneric for FP64.
 *
 * @see [[rial.arith.MultFPGeneric]]
 */
class MultFP64( stage : PipelineStageConfig )
    extends MultFPGeneric( RealSpec.Float64Spec,
      RealSpec.Float64Spec,  RealSpec.Float64Spec,
      RoundSpec.roundToEven, stage) {
}
