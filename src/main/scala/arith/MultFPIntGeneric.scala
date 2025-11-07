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
 * It takes `x` and `y` and returns `z` as `z = x + y` where y is an Int.
 *
 * All the input/output floating-point spec are configurable through [[rial.arith.RealSpec]].
 *
 * {{{
 * class MultFPIntGeneric(...) extends Module {
 *   val io = IO(new Bundle{
 *     val x = Input (UInt(xSpec.W.W))
 *     val y = Input (UInt(yWidth.W))
 *     val z = Output(UInt(zSpec.W.W))
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
 * @param stage          number of pipeline stages in this module.
 *
 * For parameters, see also:
 * @see [[rial.arith.RealSpec]]
 * @see [[rial.arith.RoundSpec]]
 * @see [[rial.util.PipelineStageConfig]]
 *
 */
class MultFPIntGeneric(
  xSpec : RealSpec, yWidth : Int, ySigned : Boolean, zSpec : RealSpec,
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  def getParam = { (xSpec, yWidth, ySigned, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle {
    val x   = Input (UInt(xSpec.W.W))
    val y   = Input (UInt(yWidth.W))
    val z   = Output(UInt(zSpec.W.W))
  })

  val (xsgn,  xex,  xman) = FloatChiselUtil.decomposeWithLeading1(xSpec, io.x)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, io.x)

//   printf("x = %b|%b(%d)|%b\n", xsgn, xex, xex-xSpec.exBias.U, xman)

  val ysgn = Wire(UInt(1.W))
  val yval = Wire(UInt(yWidth.W))

  if (ySigned) {
    ysgn := io.y(yWidth-1)
    yval := Mux(ysgn === 0.U(1.W), io.y, ~io.y + 1.U) // abs
  } else {
    ysgn := 0.U(1.W)
    yval := io.y
  }

  val yLength = (yWidth-1).U - PriorityEncoder(Reverse(yval))
  val yman    = (yval << ((yWidth-1).U - yLength))(yWidth, 0)
  val yex     = yLength
//   printf("y = %b(%d)\n", io.y, io.y)
//   printf("yWidth  = %d\n", yWidth.U)
//   printf("yLength = %d\n", yLength)
//   printf("y = %b|%b|%b\n", ysgn, yex, yman)

  val zzero = xzero.asBool || yval === 0.U
  val zinf  = xinf.asBool
  val znan  = xnan.asBool
  val zsgn0 = (xsgn ^ ysgn)

//   printf("zzero, zinf, znan, zsgn0 = %b,%b,%b,%b\n",zzero, zinf, znan, zsgn0)

  //----------------------------------------------------------------------
  // Mantissa

  val prod      = xman * yman

//   printf("xman * yman = %b\n", prod)

  val bp        = xSpec.manW+yWidth-1
  val roundBits = bp-zSpec.manW
  val moreThan2 = prod(bp+1)
  val sticky    = prod(roundBits-2,0).orR | (moreThan2 & prod(roundBits-1))
  val round     = Mux(moreThan2, prod(roundBits), prod(roundBits-1))
  val prodShift = Mux(moreThan2, prod(bp, roundBits+1), prod(bp-1, roundBits))
  val lsb       = prodShift(0)
  // RoundSpec : truncate/ceil always towards 0/infinite
  val inc = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky)
  val prodRound = prodShift +& inc
  val moreThan2AfterRound = prodRound(zSpec.manW)
  val resMan    = prodRound(zSpec.manW-1,0)

  //----------------------------------------------------------------------
  // Exponent

  val maxEx = ( maskI(xSpec.exW)-1-xSpec.exBias + yWidth - 1 + zSpec.exBias )
  val minEx = 1-xSpec.exBias+zSpec.exBias

  val exW0 = max( log2Up(maxEx+1), log2Up(abs(minEx)+1))
  val exW  = if (minEx<0) {exW0+1} else {exW0}

  val exAdd = (-xSpec.exBias + zSpec.exBias).S(exW.W).asUInt
  val ex0   = xex.pad(exW) + yex.pad(exW) + exAdd

  val exInc = ex0 + (moreThan2 | moreThan2AfterRound)

  val exZero = !exInc.orR.asBool // needed?
  val exNeg  = (minEx<0).B && (exInc(exW-1)=/=0.U)
  val exZN   = exZero || exNeg || zzero
  val exInf  = zinf || exInc(zSpec.exW-1,0).andR | exInc(exW-1, zSpec.exW).orR

  val zex = Mux(exZN,         0.U(zSpec.exW.W),
            Mux(exInf | znan, Fill(zSpec.exW,1.U(1.W)),
                              exInc(zSpec.exW-1,0)))

  // Final mantissa
  val zman = Mux(exZN || exInf || znan, znan ## 0.U((zSpec.manW-1).W), resMan)
  val zsgn = ( !(exZN || znan) ) && zsgn0.asBool
  val z0 = if (zSpec.disableSign) (zex ## zman) else (zsgn ## zex ## zman)

//   printf("z = %b|%b(%d)|%b\n", zsgn, zex, zex-zSpec.exBias.U, zman)

  io.z   := ShiftRegister(z0, nStage)
  //printf("x=%x y=%x z=%x\n", io.x, io.y, io.z)
}

/** Specialization of MultFPIntGeneric for FP64.
 *
 * @see [[rial.arith.MultFPIntGeneric]]
 */
class MultFPInt64( stage : PipelineStageConfig )
    extends MultFPIntGeneric( RealSpec.Float64Spec, 64, true, RealSpec.Float64Spec,
      RoundSpec.roundToEven, stage) {
}

