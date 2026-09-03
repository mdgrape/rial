//% @file ScaledFixedToFloatGeneric.scala
//
// Fixed point * floating point -> floating point
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.fixed

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import rial.arith._
import rial.fixed._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._


/** Calculate Fixed * Float and convert it into another Float.
 *
 * {{{
 * class ScaledFixedToFloatGeneric(...) extends Module {
 *   val io = IO(new Bundle{
 *     val x = Input (UInt(xSpec.W.W))
 *     val y = Input (UInt(w.W)) // defined only if `yintW` is defined
 *     val z = Output(UInt(zSpec.W.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * It takes a fixedpoint `x` and floating point `y` and returns a floating point `z = x * y`.
 *
 * The latency of this module can be configured via [[rial.util.PipelineStageConfig]].
 * The module declares the boundaries of its datapath, so the registers can be
 * placed inside it - either automatically or by name:
 *
 *  - `"in"`    the module input
 *  - `"norm"`  after the absolute value, the leading zero count and the
 *              normalisation shift, and after the decompose / classify of the
 *              floating point operand that runs in parallel with them
 *  - `"mul"`   after the mantissa product
 *  - `"round"` after the rounding increment and the exponent sum
 *  - `"out"`   the module output
 *
 * {{{
 * // three registers, placed by the module
 * new ScaledFixedToFloatGeneric(xSpec, ySpec, zSpec, RoundSpec.roundToEven, 3)
 * // hand-picked positions
 * new ScaledFixedToFloatGeneric(xSpec, ySpec, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.specify("norm" -> 1, "mul" -> 1, "round" -> 1))
 * // everything at the output, as before
 * new ScaledFixedToFloatGeneric(xSpec, ySpec, zSpec, RoundSpec.roundToEven,
 *   PipelineStageConfig.atOut(3))
 * }}}
 *
 * @param xSpec spec of the input fixed point number.
 * @param ySpec spec of the input floating point number.
 * @param zSpec spec of the output floating point number.
 * @param roundSpec rounding algorithm.
 * @param stage latency of this module.
 *
 * For parameters, see also:
 * @see [[rial.fixed.FixedSpec]]
 * @see [[rial.arith.RealSpec]]
 * @see [[rial.arith.RoundSpec]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class ScaledFixedToFloatGeneric(
  xSpec : FixedSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  /** Boundaries of the datapath at which pipeline registers can be placed.
   * The weights are rough, target-independent estimates of the delay of the
   * logic that ends at each boundary. */
  val stagePoints = ScaledFixedToFloatGeneric.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getParam = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x = Input (UInt(xSpec.W.W))
    val y = Input (UInt(ySpec.W.W))
    val z = Output(UInt(zSpec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  val xIn = ShiftRegister(io.x, stages.at("in"))
  val yIn = ShiftRegister(io.y, stages.at("in"))

  //----------------------------------------------------------------------
  // S1 : decompose / classify the floating point operand and, in parallel,
  //      the absolute value / leading zero count / normalisation shift of
  //      the fixed point operand
  val (ysgn,  yex,  ymanW1) = FloatChiselUtil.decomposeWithLeading1(ySpec, yIn)
  val (yzero, yinf, ynan)   = FloatChiselUtil.checkValue(ySpec, yIn)

  val xzero = ~xIn.orR
  val xsgn  = Wire(UInt(1.W))
  val xabs  = Wire(UInt(xSpec.W.W))
  if (xSpec.signed) {
    xsgn := xIn(xSpec.W-1)
    xabs := Mux(xIn(xSpec.W-1), ~xIn + 1.U, xIn)
  } else {
    xsgn := 0.U(1.W)
    xabs := xIn
  }

  val xclz = PriorityEncoder(Reverse(xabs))
  val xShifted = (xabs << xclz)(xSpec.W-1, 0)

  //---------------------------- stage "norm" ----------------------------
  private def cutNorm[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("norm"))
  }

  val xShiftedNorm = cutNorm(xShifted)
  val xclzNorm     = cutNorm(xclz)
  val xsgnNorm     = cutNorm(xsgn)
  val xzeroNorm    = cutNorm(xzero)
  val ymanW1Norm   = cutNorm(ymanW1)
  val yexNorm      = cutNorm(yex)
  val ysgnNorm     = cutNorm(ysgn)
  val yzeroNorm    = cutNorm(yzero)
  val yinfNorm     = cutNorm(yinf)
  val ynanNorm     = cutNorm(ynan)

  // -------------------------------------------------------------------------
  // S2 : mantissa product

  val prod = xShiftedNorm * ymanW1Norm
  val prodWidth = xSpec.W + ySpec.manW+1
  assert(prodWidth == prod.getWidth)

  //---------------------------- stage "mul" -----------------------------
  private def cutMul[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("mul"))
  }

  val prodMul  = cutMul(prod)
  val xclzMul  = cutMul(xclzNorm)
  val xsgnMul  = cutMul(xsgnNorm)
  val xzeroMul = cutMul(xzeroNorm)
  val yexMul   = cutMul(yexNorm)
  val ysgnMul  = cutMul(ysgnNorm)
  val yzeroMul = cutMul(yzeroNorm)
  val yinfMul  = cutMul(yinfNorm)
  val ynanMul  = cutMul(ynanNorm)

  // -------------------------------------------------------------------------
  // S3 : round / sticky selection and the rounding increment; in parallel,
  //      the exponent sum

  val prodMoreThan2 = prodMul(prodWidth - 1)

  //     decimal point
  //     |
  //     v  xSpec.W + y.manW + 1 = prodWidth
  // .---+------------------------------.
  // | | |          L|R|S.....          |
  //  ^'-------------' '----------------'
  //  1 :   z.manW    :     sticky
  //  : :             :
  //  : :             roundBit = prodWidth - 1 - manW
  //  : prodWidth - 2
  //  prodWidth - 1

  val roundBit = prodWidth - zSpec.manW - 1

  val lsb    = Mux(prodMoreThan2, prodMul(roundBit),   prodMul(roundBit-1))
  val round  = Mux(prodMoreThan2, prodMul(roundBit-1), prodMul(roundBit-2))
  val sticky = (prodMoreThan2 & prodMul(roundBit-2)) | prodMul(roundBit-3, 0).orR

  val man0 = Mux(prodMoreThan2, prodMul(prodWidth-2, roundBit),
                                prodMul(prodWidth-3, roundBit-1))
  val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb.asBool, round.asBool, sticky)
  val manRounded = man0 +& inc.asUInt

  val moreThan2AfterRound = manRounded(zSpec.manW)

  // -------------------------------------------------------------------------
  // exponent

  val maxEx = xSpec.intW-1 + ySpec.exMax
  val minEx = -xSpec.fracW + ySpec.exMin

  val zExInc = (prodMoreThan2 | moreThan2AfterRound).asUInt
  val yExNobias = yexMul.zext - ySpec.exBias.S
  val xExNobias = (xSpec.intW - 1).S - xclzMul.zext
  val zExNobias = yExNobias +& xExNobias + zExInc.zext
  val zEx0      = zExNobias + zSpec.exBias.S

  val znan  = ynanMul || (yinfMul && xzeroMul)
  val zinf  = (zSpec.exMax.S < zExNobias) || yinfMul
  val zzero = yzeroMul || xzeroMul || zEx0 <= 0.S

  //--------------------------- stage "round" ----------------------------
  private def cutRound[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("round"))
  }

  val zManBitsRound = cutRound(manRounded(zSpec.manW-1, 0))
  val zExBitsRound  = cutRound(zEx0(zSpec.exW-1, 0))
  val znanRound     = cutRound(znan)
  val zinfRound     = cutRound(zinf)
  val zzeroRound    = cutRound(zzero)
  val xsgnRound     = cutRound(xsgnMul)
  val ysgnRound     = cutRound(ysgnMul)

  // -------------------------------------------------------------------------
  // S4 : merge (special case muxes)

  val zEx = Mux(zzeroRound,            0.U(zSpec.exW.W),
            Mux(zinfRound | znanRound, Fill(zSpec.exW, 1.U(1.W)),
                                       zExBitsRound))

  val zMan = Mux(zzeroRound,            0.U(zSpec.manW.W),
             Mux(zinfRound | znanRound, Cat(znanRound.asUInt, Fill(zSpec.manW-1, 0.U(1.W))),
                                        zManBitsRound))

  val z0 = if (zSpec.disableSign) {
    Cat(zEx, zMan)
  } else {
    val zSgn = xsgnRound ^ ysgnRound
    Cat(zSgn, zEx, zMan)
  }
//   printf("z0    = %b(w=%d)\n", z0 , z0.getWidth.U)

  io.z := ShiftRegister(z0, stages.at("out"))
}

/** Pipeline stage points of [[rial.fixed.ScaledFixedToFloatGeneric]].
 *
 * @see [[rial.util.StagePoint]]
 */
object ScaledFixedToFloatGeneric {
  /** Ordered boundaries of the datapath, with a rough, target-independent
   *  estimate of the delay of the logic ending at each. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("norm",  3), StagePoint("mul", 4),
    StagePoint("round", 4), StagePoint("out", 1))
}
