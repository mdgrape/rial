//% @file AddFPGeneric.scala
//
// Generic Floating-Point Adder Unit
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.arith

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import chisel3.util.experimental.BoringUtils
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.util.DebugControlSlave

class AddFPGeneric(
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig,
  val enableDebug : Boolean = false
) extends Module with DebugControlSlave {

  val nStage = stage.total

  def getParam() = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage() = nStage

  val io = IO(new Bundle{
    val x   = Input(UInt(xSpec.W.W))
    val y   = Input(UInt(ySpec.W.W))
    val z   = Output(UInt(zSpec.W.W))
  })

  val (xsgn, xex, xman) = FloatChiselUtil.decompose(xSpec, io.x)
  val (ysgn, yex, yman) = FloatChiselUtil.decompose(ySpec, io.y)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, io.x)
  val (yzero, yinf, ynan) = FloatChiselUtil.checkValue(ySpec, io.y)
  val diffSgn = (xsgn ^ ysgn).asBool

  val xyInf  = xinf || yinf
  val xyInfSgn = ! ((xinf && !xsgn) || (yinf && !ysgn))
  val xyNaN  = xnan || ynan
  val xyBothZero = xzero && yzero
  dbgPrintf("%x %x\n",  io.x, io.y)
  dbgPrintf("%b %b %b\n",  xyInf, xyNaN, xyBothZero)

  //----------------------------------------------------------------------
  // Exponent comparison
  val diffMaxXY = xSpec.exMax - ySpec.exMin // x:max, y:min
  val diffMaxYX = ySpec.exMax - xSpec.exMin // y:max, x:min
  val diffExW = log2Up(max(diffMaxXY,diffMaxYX)+1)+1 
  val diffExXY = ( xex.pad(diffExW).asSInt-yex.pad(diffExW).asSInt
    +(-xSpec.exBias+ySpec.exBias).S(diffExW.W) )
  // The following is just negative of diffExXY but in some case the separate calculation
  // decreases the latency
  val diffExYX = ( yex.pad(diffExW).asSInt-xex.pad(diffExW).asSInt
    +(xSpec.exBias-ySpec.exBias).S(diffExW.W) )
  val near = (!xzero) && (!yzero) && diffSgn && ((diffExXY(diffExW-1,1) === 0.U)|| (diffExXY=== -1.S)) // abs(xEx-yEx) <= 1
  val nearShift = diffExXY(0)
  val exXlargerThanY = diffExXY(diffExW-1) === 0.U // exponent of x >= y

  val zexBaseMax = max( xSpec.exMax + zSpec.exBias, ySpec.exMax + zSpec.exBias)
  val zexBaseMin = max( xSpec.exMin + zSpec.exBias, ySpec.exMin + zSpec.exBias)
  val zexBaseW = getWidthToRepresentNumbersAlwaysSigned( Seq(zexBaseMin,  zexBaseMax) )
  //println(f"diffExW=$diffExW%d zexBaseW = $zexBaseW%d")
  val zexNorm = Wire(SInt(zexBaseW.W))
  if ((zSpec.exBias == xSpec.exBias) && (zSpec.exBias == ySpec.exBias)) {
    zexNorm := (Mux (exXlargerThanY, xex, yex)).zext
  } else {
    val xexPad = xex.pad(zexBaseW).asSInt
    val yexPad = yex.pad(zexBaseW).asSInt
    val zexDiff = Mux (exXlargerThanY,
      (zSpec.exBias-xSpec.exBias).S(zexBaseW.W),
      (zSpec.exBias-ySpec.exBias).S(zexBaseW.W) )
    zexNorm := Mux (exXlargerThanY, xexPad, yexPad) + zexDiff
  }
  dbgPrintf("near=%b zexNorm=%x diffExXY/XY=%x / %x\n", near, zexNorm, diffExXY, diffExYX)

  //------------------------------------------------------------------------
  // Mantissa padding
  val minManW    = xSpec.manW.min(ySpec.manW)
  val maxManWxy  = xSpec.manW.max(ySpec.manW)
  val maxManWxyz = maxManWxy.max(zSpec.manW)
  //val manXlargerThanY = xman(xSpec.manW-1, xSpec.manW-minManW) > yman(ySpec.manW-1, ySpec.manW-minManW)

  val xmanPadXyz = Mux(xzero, 0.U((maxManWxyz+1).W), 1.U(1.W) ## padLSB(maxManWxyz, xman))
  val ymanPadXyz = Mux(yzero, 0.U((maxManWxyz+1).W), 1.U(1.W) ## padLSB(maxManWxyz, yman))
  

  //----------------------------------------------------------------------
  // Far path
  //   when exponent difference > 1 or sign is the same
  //   in this case, assume x is normalized as 1<=x<2 and abs(y)<x,
  //     -1/2 < y < 2, and then 1/2 < x+y < 4.
  //
  //   currently, separate barrel shifters/adders are used for far/near paths.
  val xFar = 0.U(1.W) ## Mux(exXlargerThanY, xmanPadXyz, ymanPadXyz)
  val yFar = Mux(exXlargerThanY, ymanPadXyz, xmanPadXyz)
  val yzeroFar = Mux(exXlargerThanY, xzero, yzero)
  val xFarSign = Mux(exXlargerThanY, xsgn, ysgn).asBool
  val yFarR = yFar ## 0.U(2.W)
  // Drop sign bit 
  val diffFar = Mux(exXlargerThanY, diffExXY(diffExW-2,0), diffExYX(diffExW-2,0))

  // Shift y
  //   Actually, shift amount here is always larger than or equal to 2.
  val shiftW0 = log2Up(maxManWxyz+1+2+1) // +2 for round bits, +1 for Leading 1
  // Note diffExW includes the width for sign bit
  val (shiftW, shiftOut) =
    if   (shiftW0<(diffExW-2)) (shiftW0, dropLSB(shiftW0,diffFar).orR)
    else (diffExW-1, false.B) // Usually this never happens
  val yFarShift = Mux(shiftOut, 0.U, yFarR >> diffFar(shiftW-1, 0))
  println(f"shiftW = $shiftW%d diffExW = $diffExW%d ${diffFar.getWidth}%d")
  when (!near) {
    dbgPrintf("xman=%x yman=%x\n",  xman, yman)
    dbgPrintf("xmanPadXyz=%x ymanPadXyz=%x\n", xmanPadXyz, ymanPadXyz)
    dbgPrintf("xFar=%x yFarShift=%x\n", xFar, yFarShift)
    //dbgPrintf("%b\n%b\n", yFar, yFarShift)
  }
  val diffFarMinus2    = diffFar(shiftW-1,0) - 2.U
  val diffFarLessThan2 = !diffFar(diffFar.getWidth-1,shiftW-1).orR && diffFarMinus2(shiftW-1)

  // Y shifted = (main value)-round0-round1-stickey
  //    these round0/round1 can be combined to stickey
  //    when zSpec.manW is smaller than x/ySpec.manW.
  //    However, we expect such situations are rare.
  //
  // Get stickey bit
  //   if stickey bit = Or ( all bits [shift-3:0] )
  //   PriorityEncoder in Chisel returns bit position of 1 from lsb
  //     0 returns all 1; unfortunately if bit width of input is 2^w,
  //     2^(w-1) will return the same result as 0.
  val onePosFromLSB = PriorityEncoder(yFar)
  val yFarW = yFar.getWidth
  // Since always a leading 1 bit exists, stickey == 1 when shiftOut, but round should be 0
  val stickeyFar = ( onePosFromLSB < diffFarMinus2) && !diffFarLessThan2
  // Increment for negative value occurs only rounded bits are all 0
  val guardFarNeg = (~(yFarShift(1,0) ## stickeyFar)) +& 1.U
  val guardFar    = Mux(diffSgn, guardFarNeg(2,1), yFarShift(1,0))
  val negIncFar   = diffSgn && guardFarNeg(3)
  val yAddFar = diffSgn ## Mux(diffSgn,~yFarShift(maxManWxyz+2,2),yFarShift(maxManWxyz+2,2))
  // Width = maxManWxyz+2
  val sumFar = xFar + yAddFar + negIncFar
  val sumFarWithGuard = sumFar ## guardFar
  when (!near) {
    dbgPrintf("yFar=%x onePosFromLSB=%d diffFar=%d(%x) shiftOut=%d\n", yFar, onePosFromLSB, diffFar, diffFar,shiftOut)
    dbgPrintf("guardFar=%x stickeyFar=%d\n", guardFar, stickeyFar)
    dbgPrintf("diffSgn=%d sumFar=%x yAddFar=%x negIncFar=%b\n", diffSgn, sumFar, yAddFar, negIncFar)
  }

  // Rounding
  val zmanFar     = Wire(UInt(zSpec.manW.W)) // This omit leading zero
  val roundFarSel = Wire(Bool())
  val additionalStickeyFar = Wire(Bool())
  val exSubFar    = Wire(SInt(2.W))
  val diffWz = maxManWxyz-zSpec.manW // >=0
  when (sumFar(maxManWxyz+1).asBool) { // >=2
    zmanFar      := sumFar(maxManWxyz, diffWz+1)
    roundFarSel  := sumFar(diffWz)
    additionalStickeyFar := orRLSB(diffWz+2, sumFarWithGuard)
    exSubFar     := -1.S
  }.elsewhen (sumFar(maxManWxyz)===1.U) { // >=1
    zmanFar      := sumFar(maxManWxyz-1, diffWz)
    roundFarSel  := sumFarWithGuard(diffWz+1)
    additionalStickeyFar := orRLSB(diffWz+1, sumFarWithGuard)
    exSubFar     := 0.S
  }.otherwise { // >=1/2
    zmanFar      := sumFarWithGuard(maxManWxyz, diffWz+1)
    roundFarSel  := sumFarWithGuard(diffWz)
    additionalStickeyFar := orRLSB(diffWz, sumFarWithGuard)
    exSubFar     := 1.S
  }
  when (!near) { dbgPrintf("zmanFar=%x\n", zmanFar) }

  val incFar = FloatChiselUtil.roundIncBySpec(roundSpec, zmanFar(0), roundFarSel,
    stickeyFar || additionalStickeyFar)

  //----------------------------------------------------------------------
  // Near path
  //   when exponent difference <= 1 and sign is different
  //   in this case, assume x is normalized as 1<=x<2 and abs(y)<x,
  //     -1/2 >= y > -2, and 0 < x+y < 1+1/2.
  //   if compare only exponent, the result can be negative.
  //     in this case, -1 < x+y < 1+1/2.
  //----------------------------------------------------------------------

  // Implementation without mantissa comparison - use XFar/YFar
  val xNear = getMSB(2+maxManWxy, xFar) ## 0.U(1.W)
  val yNear = getMSB(2+maxManWxy, yFar)
  
  val yNearShift = Mux( nearShift, 0.U(2.W) ## yNear, 0.U(1.W) ## yNear ## 0.U(1.W) )
  val sumNear = xNear - yNearShift
  // The below includes 1.xxxx+round bit
  //  width = maxManWxy+2
  val sumNearSign = getMSB1(0,sumNear)
  val sumNearAbs = Mux(sumNearSign, -sumNear(maxManWxy+1,0), sumNear(maxManWxy+1,0))

  val shiftSum    = PriorityEncoder(Reverse(sumNearAbs))
  // 0 .. >=1
  // 1 .. >=1/2 ...
  val sumNearNorm = (sumNearAbs << shiftSum)(maxManWxy,0) // include (possible) round bit
  val sumNearZero = !sumNear.orR
  val nearSign    = (!sumNearZero) && (sumNearSign ^ xFarSign)

  // Rounding
  val zmanNear = Wire(UInt((zSpec.manW).W))
  val incNear  = Wire(Bool())
  if (zSpec.manW < maxManWxy) { // Output precision is lower
    val roundNear   = getMSB1(zSpec.manW+1, sumNearNorm)
    val stickeyNear = orRLSB(maxManWxy-zSpec.manW, sumNearNorm) // This is usually zero
    incNear  := FloatChiselUtil.roundIncBySpec(roundSpec, getMSB1(zSpec.manW,sumNearNorm), roundNear,
      stickeyNear)
    zmanNear := sumNearNorm(maxManWxy, maxManWxy-zSpec.manW+1)
  } else if (zSpec.manW == maxManWxy) { // Output precision is the same as input
    incNear  := FloatChiselUtil.roundIncBySpec(roundSpec, sumNearNorm(1), sumNearNorm(0), false.B)
    zmanNear := sumNearNorm(maxManWxy, 1)
  } else { // Output precision is higher, no rounding
    incNear  := 0.U(1.W)
    zmanNear := padLSB( zSpec.manW, sumNearNorm(maxManWxy, 0) )
  }
  when (near) {
    dbgPrintf("xNear=%x yNear=%x nearShift=%d yNearShift=%x\n", xNear, yNear, nearShift, yNearShift) 
    dbgPrintf("sumNear=%x sumNearAbs=%x\n", sumNear, sumNearAbs)
    dbgPrintf("sumNearAbs=%b\n", sumNearAbs)
    dbgPrintf("shiftSum=%d(%x) sumNearNorm=%x incNear=%b\n", shiftSum, shiftSum, sumNearNorm, incNear)
  }

  //----------------------------------------------------------------------
  // Merge Far / Near 
  val zmanSel = Mux(near, zmanNear, zmanFar)
  val incSel  = Mux(near, incNear, incFar)
  val zeroSel = near && sumNearZero
  val exSub   = Mux(near, shiftSum.zext, exSubFar) // SInt

  val zman0   = zmanSel +& incSel

  // Calculate zex Range
  val zexMax = zexBaseMax+2
  // In cases when zSpec.manW is smaller, it will increase 2 (result become 4)
  val zexMin = zexBaseMin - sumNearAbs.getWidth
  val zexW   = getWidthToRepresentNumbersAlwaysSigned( Seq(zexMin,  zexMax) )
  val zex0   = zexNorm.pad(zexW) - exSub + zman0(zSpec.manW).zext
  val zeroAfterAdd = (zex0 < 0.S) || zeroSel
  val infAfterAdd = zex0 >= maskI(zSpec.exW).S

  val zsgn0 = Mux(near, nearSign, xFarSign)
  //----------------------------------------------------------------------
  // Final
  val infOrNaN = xyInf || xyNaN || infAfterAdd
  val zero     = zeroAfterAdd||xyBothZero
  val zex = Mux (infOrNaN, maskU(zSpec.exW),
    Mux (zero, 0.U(zSpec.exW), zex0(zSpec.exW-1,0)))
  val zman = Mux(infOrNaN || zero, xyNaN ## 0.U((zSpec.manW-1).W), zman0(zSpec.manW-1, 0))
  val zsgn = (!xyInf) && Mux(xyInf, xyInfSgn, zsgn0)

  val z0 = zsgn ## zex ## zman

  //val z0 = 0.U(zSpec.W.W)

  io.z   := ShiftRegister(z0, nStage)
  dbgPrintf("x=%x y=%x z=%x\n", io.x, io.y, io.z)
}

class AddFP64( stage : PipelineStageConfig )
    extends AddFPGeneric( RealSpec.Float64Spec,
      RealSpec.Float64Spec,  RealSpec.Float64Spec,
      RoundSpec.roundToEven, stage) {
}

object AddFP64_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new AddFP64(PipelineStageConfig.none)) ) )
}
