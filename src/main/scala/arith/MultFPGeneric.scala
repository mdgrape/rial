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
 * The latency of this module can be configured via [[rial.util.PipelineStageConfig]].
 *
 * {{{
 * class MultFPGeneric(...) extends Module {
 *   val io = IO(new Bundle{
 *     val x = Input (UInt(xSpec.W.W))
 *     val y = Input (UInt(ySpec.W.W))
 *     val z = Output(UInt(zSpec.W.W))
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
 * @param stage          number of pipeline stages in this module.
 *
 * For parameters, see also:
 * @see [[rial.arith.RealSpec]]
 * @see [[rial.arith.RoundSpec]]
 * @see [[rial.util.PipelineStageConfig]]
 *
 */
class MultFPGeneric(
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  def getParam = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle {
    val x   = Input(UInt(xSpec.W.W))
    val y   = Input(UInt(ySpec.W.W))
    val z   = Output(UInt(zSpec.W.W))
  })

  val (xsgn, xex, xman) = FloatChiselUtil.decomposeWithLeading1(xSpec, io.x)
  val (ysgn, yex, yman) = FloatChiselUtil.decomposeWithLeading1(ySpec, io.y)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, io.x)
  val (yzero, yinf, ynan) = FloatChiselUtil.checkValue(ySpec, io.y)
  val zsgn0 = xsgn ^ ysgn

  val zzero = (xzero | yzero).asBool
  val zinf  = (xinf | yinf).asBool
  val znan  = (xnan | ynan).asBool & !zzero
  //printf("%x %x\n",  io.x, io.y)
  //printf("%b %b %b\n",  zzero, zinf, znan)

  //----------------------------------------------------------------------
  // Mantissa
  val prod = xman * yman
  val bp = xSpec.manW+ySpec.manW
  val roundBits = bp-zSpec.manW
  val moreThan2 = prod(bp+1)

  val (moreThan2AfterRound, resMan) =
    if(roundBits >= 0) {
      val stickey = if(roundBits >= 2) {
        prod(roundBits-2,0).orR | (moreThan2 & prod(roundBits-1))
      } else if(roundBits == 0) {
        moreThan2 & prod(roundBits-1)
      } else {
        0.U(1.W).asBool
      }
      val round = if(roundBits >= 1) {
        Mux(moreThan2, prod(roundBits), prod(roundBits-1))
      } else {
        moreThan2 & prod(0)
      }
      val prodShift = Mux(moreThan2, prod(bp, roundBits+1), prod(bp-1, roundBits))
      val lsb   = prodShift(0)
      // RoundSpec : truncate/ceil always towards 0/infinite
      val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, stickey)
      val prodRound = prodShift +& inc
      val moreThan2AfterRound = prodRound(zSpec.manW)
      (moreThan2AfterRound, prodRound(zSpec.manW-1, 0))
    } else {
      val prodShift = Mux(moreThan2, prod, Cat(prod(prod.getWidth-2, 0), 0.U(1.W)))
      assert(prodShift.getWidth == prod.getWidth)
      assert(prodShift(prodShift.getWidth-1) === 1.U(1.W))
      assert(-1 >= roundBits)
      val prodRound = Cat(prodShift, Fill(-roundBits-1, 0.U(1.W)))
      assert(prodRound(zSpec.manW) === 1.U(1.W))
      (0.U(1.W), prodRound(zSpec.manW-1, 0))
    }

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
  val ex0 = xex.pad(exW) + yex.pad(exW) + exAdd

  val exInc = ex0 + (moreThan2 | moreThan2AfterRound)

  val exZero = !exInc.orR.asBool
  val exNeg  = (minEx<0).B && (exInc(exW-1)=/=0.U)
  val exZN   = exZero || exNeg || zzero
  val exInf  = if(exW > zSpec.exW) {
    zinf || exInc(zSpec.exW-1,0).andR | exInc(exW-1, zSpec.exW).orR
  } else if(exW == zSpec.exW) {
    zinf || exInc(zSpec.exW-1,0).andR
  } else {
    zinf
  }

  val zex =
    Mux(exZN, 0.U(zSpec.exW.W),
      Mux(exInf | znan, Fill(zSpec.exW,1.U(1.W)), exInc(zSpec.exW-1,0)) )

  // Final mantissa
  val zman =
    Mux(exZN || exInf || znan, znan ## 0.U((zSpec.manW-1).W), resMan)
  val zsgn = ( !(exZN || znan) ) && zsgn0.asBool
  //printf("%x\n", zman)

  assert(zman.getWidth == zSpec.manW)

  val z0 = if (zSpec.disableSign) {(zex ## zman)} else {(zsgn ## zex ## zman)}

  io.z   := ShiftRegister(z0, nStage)
  //printf("x=%x y=%x z=%x\n", io.x, io.y, io.z)
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

