//% @file FusedMulAddFPGeneric.scala
//
// Generic Floating-Point FMA Unit
// Copyright (C) Toru Niina RIKEN BDR 2021
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

//
// w <= x * y + z
//
class FusedMulAddFPGeneric(
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, wSpec : RealSpec,
  roundSpec : RoundSpec,
  stage : PipelineStageConfig,
  val enableDebug : Boolean = false
) extends MultiIOModule with DebugControlSlave {

  val nStage = stage.total

  def getParam() = { (xSpec, ySpec, zSpec, wSpec, roundSpec, nStage) }

  def getStage() = nStage

  val io = IO(iodef = new Bundle {
    val x   = Input(UInt(xSpec.W.W))
    val y   = Input(UInt(ySpec.W.W))
    val z   = Input(UInt(zSpec.W.W))
    val w   = Output(UInt(wSpec.W.W))
  })

  val (xsgn, xex, xman) = FloatChiselUtil.decomposeWithLeading1(xSpec, io.x)
  val (ysgn, yex, yman) = FloatChiselUtil.decomposeWithLeading1(ySpec, io.y)
  val (zsgn, zex, zman) = FloatChiselUtil.decomposeWithLeading1(zSpec, io.z)

  dbgPrintf("x   = %d|%d|%d\n", xsgn, xex, xman)
  dbgPrintf("y   = %d|%d|%d\n", ysgn, yex, yman)
  dbgPrintf("z   = %d|%d|%d\n", zsgn, zex, zman)


  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, io.x)
  val (yzero, yinf, ynan) = FloatChiselUtil.checkValue(ySpec, io.y)
  val (zzero, zinf, znan) = FloatChiselUtil.checkValue(zSpec, io.z)

  val xysgn  = xsgn ^ ysgn
  val xyzop  = (xysgn ^ zsgn).asBool // 0->(+), 1->(-)
  val xyzero = xzero || yzero
  val xyinf  = xinf  || yinf
  val xynan  = xnan  || ynan || (xinf && yzero) || (xzero && yinf)
  val xyznan = xynan || znan || (xyinf && zinf && xyzop)
  val xyzinf = (xyinf && !znan && !zinf) || (xyinf && zinf && !xyzop)

  // -------------------------------------------------------------------------
  // exponent check

  val xyMaxEx  = xSpec.exMax + ySpec.exMax
  val xyMinEx  = xSpec.exMin + ySpec.exMin
  val xyExW0   = max(log2Up(xyMaxEx+1), log2Up(abs(xyMinEx) + 1))
  val xyExW    = if (xyMinEx<0) {xyExW0+1} else {xyExW0}
  val xyExBias = xSpec.exBias + ySpec.exBias

  dbgPrintf(f"xyMaxEx  = $xyMaxEx%d = xMaxEx(${xSpec.exMax}%d) + yMaxEx(${ySpec.exMax}%d)")
  dbgPrintf(f"xyMinEx  = $xyMinEx%d = xMaxEx(${xSpec.exMin}%d) + yMaxEx(${ySpec.exMin}%d)")
  dbgPrintf(f"xyExW0   = $xyExW0%d")
  dbgPrintf(f"xyExW    = $xyExW%d")
  dbgPrintf(f"xyExBias = $xyExBias%d")

  // nobias
  val xyEx0    = xex.pad(xyExW).asSInt + yex.pad(xyExW).asSInt - xyExBias.S(xyExW.W)

  dbgPrintf("xex   = %d\n", xex.pad(xyExW))
  dbgPrintf("yex   = %d\n", yex.pad(xyExW))
  dbgPrintf("xyEx0 = %d\n", xyEx0)
  dbgPrintf("zex   = %d\n", zex)

  val diffMaxXYMinZ = xyMaxEx - zSpec.exMin
  val diffMaxZMinXY = zSpec.exMax - xyMinEx
  val diffExW       = log2Up(max(diffMaxXYMinZ, diffMaxZMinXY) + 1) + 1

  dbgPrintf("zexnb = %d\n", zex.pad(diffExW).asSInt - zSpec.exBias.S(diffExW.W))

  dbgPrintf(f"diffMaxXYMinZ = $diffMaxXYMinZ%d")
  dbgPrintf(f"diffMaxZMinXY = $diffMaxZMinXY%d")
  dbgPrintf(f"diffExW       = $diffExW%d")

  val diffExXYZ = xyEx0.pad(diffExW).asSInt - zex.pad(diffExW).asSInt + zSpec.exBias.S(diffExW.W)
  val diffExZXY = zex.pad(diffExW).asSInt - xyEx0.pad(diffExW).asSInt - zSpec.exBias.S(diffExW.W)

  // later we may take abs of diffEx, so we calculate + and - beforehand.
  dbgPrintf("diffExXYZ = %d\n", diffExXYZ)
  dbgPrintf("diffExZXY = %d\n", diffExZXY)

  // zEx - xyEx == +2, +1, 0, -1 (for detail, see "naer path")
  val is_near = (!xyzero) && (!zzero) && (xysgn =/= zsgn) &&
                ((diffExZXY === 2.S) || (diffExZXY ===  1.S) ||
                 (diffExZXY === 0.S) || (diffExZXY === -1.S))

  val nearShift = diffExXYZ(0) // (+/-)1 or 0
  val exXYLargerThanZ = diffExXYZ(diffExW-1) === 0.U // exponent of xy >= z

  // -------------------------------------------------------------------------
  // multiply

  val xyprod = Mux(xyzero, 0.U((1 + xSpec.manW + 1 + ySpec.manW).W), xman * yman);

  // -------------------------------------------------------------------------
  // far path (product is larger)
  //
  // xyEx > zEx, so deltaEx = xyEx - zEx = diffExXYZ

  val zShiftW0 = log2Up(max(xSpec.manW+ySpec.manW, zSpec.manW) + 1 + 2 + 1)

  val (zShiftW, zShiftOut) = if (zShiftW0 < (diffExW-2)) {
      (zShiftW0, dropLSB(zShiftW0, diffExXYZ(diffExW-2, 0)).orR)
    } else {
      (diffExW-1, false.B)
    }

  // align addend z to (1+x.Spec.manW + 1+y.Spec.manW +1(carry))
  //
  //     .--x.W+y.W--.
  // 0**.************
  // 00z.zzzz00000000
  //     '--'
  //     z.W
  //
  // FIXME: if zSpec.manW > xSpec.manW + ySpec.manW ...?
  //
  val zManPad = 0.U(2.W) ## zman ## 0.U((xSpec.manW + ySpec.manW - zSpec.manW).W)
  val zManInverted = Mux(xyzop, ~zManPad+1.U, zManPad)
  val zManAligned  = Mux(zShiftOut, 0.U, zManInverted >> diffExXYZ(zShiftW-1, 0));
  val sumFarProd   = zManAligned + (0.U(1.W) ## xyprod)

  // -------------------------------------------------------------------------
  // far path (addend is larger)
  //
  // xyEx < zEx, so deltaEx = zEx - xyEx = diffExZXY


  // -------------------------------------------------------------------------
  // near path (signs are different)
  //
  // since x*y is a product, the integer part could be larger than 2.
  //
  //  dEx === 2  |  dEx === 1 |  dEx === 0 | dEx === -1  (dEx = zEx - xyEx)
  //   11.xxxxxx |  11.xxxxxx |  01.xxxxxx |  01.xxxxxx
  // -zzz.zzz    | -zz.zzzz   |  -z.zzzzz  |  -0.zzzzzz
  //
  //  2+xmanW+ymanW
  // 0011.xxxxxxxxx
  // 000z.zzz000000
  //  3+zmanW+Pad = 2+xmanW+ymanW
  //          Pad = xmanW+ymanW - zmanW - 2
  //
  // TODO: the case when wSpec.manW > xSpec.manW + ySpec.manW
  //            and when zSpec.manW + 1 > xSpec.manW + ySpec.manW
  //       the intermediate value width might not be correct because FMA
  //       requires intermediate values to have infinite precision.
  //       e.g. if zSpec.manW > xSpec.manW + ySpec.manW, then we need to add
  //            more bits to the least significant part to avoid truncation of
  //            LSBs of z while shifting

  val zNearWidth     = 2+max((1+xSpec.manW)+(1+ySpec.manW), 1+wSpec.manW)
  val zNearFracWidth = xSpec.manW + ySpec.manW
  val zNearPad       = 0.U(3.W) ## zman ## 0.U((zNearWidth - (1+zSpec.manW) - 3).W)
  val zNearShift     = Mux(diffExZXY(diffExW-1) === 0.U, (zNearPad << diffExZXY(1,0))(zNearWidth-1, 0), (0.U(1.W) ## zNearPad >> 1));
  val xyzNearSum     = (0.U(2.W) ## xyprod) - zNearShift
  val xyzNearSign    = Mux(xysgn.asBool, ~xyzNearSum(zNearWidth-1), xyzNearSum(zNearWidth-1))
  val xyzNearAbs     = Mux(xyzNearSum(zNearWidth-1), -xyzNearSum, xyzNearSum)
  val shiftNearSum   = PriorityEncoder(Reverse(xyzNearAbs))
  dbgPrintf(f"zNearWidth     = $zNearWidth%d")
  dbgPrintf(f"zNearFracWidth = $zNearFracWidth%d")
  dbgPrintf("zNearPad       = 0b%b(%d), width= %d\n",   zNearPad   , zNearPad   , zNearPad   .getWidth.U)
  dbgPrintf("zNearShift     = 0b%b(%d), width= %d\n",   zNearShift , zNearShift , zNearShift .getWidth.U)
  dbgPrintf("xyProd         = 0b  %b(%d), width= %d\n", xyprod     , xyprod     , xyprod     .getWidth.U)
  dbgPrintf("xyzNearSum     = 0b%b(%d), width= %d\n",   xyzNearSum , xyzNearSum , xyzNearSum .getWidth.U)
  dbgPrintf("xyzNearSign    = 0b%b(%d), width=%d\n",    xyzNearSign, xyzNearSign, xyzNearSign.getWidth.U)
  dbgPrintf("xyzNearAbs     = 0b%b(%d), width=%d\n",    xyzNearAbs , xyzNearAbs , xyzNearAbs .getWidth.U)
  dbgPrintf("shiftNearSum   = %b, width= %d\n", shiftNearSum , shiftNearSum.getWidth.U)

  val xyzNearNorm = (xyzNearAbs << shiftNearSum)(zNearWidth-1, 0)
  val xyzNearExDec = shiftNearSum.asSInt - 3.S
  dbgPrintf("xyzNearNorm  = %d, width=%d\n", xyzNearNorm , xyzNearNorm .getWidth.U)
  dbgPrintf("xyzNearExDec = %d, width=%d\n", xyzNearExDec, xyzNearExDec.getWidth.U)

  // 0xxx.xxxxx   : shiftNearSum 1, ex += 2 = 3 - shiftNearSum
  // 00xx.xxxxxx  : shiftNearSum 2, ex += 1
  // 000x.xxxxxxx : shiftNearSum 3, ex += 0
  // 0000.xxxxxxxx: shiftNearSum 4, ex += -1
  // |
  // v
  // 1xxx xxxxxxxx

  val wmanNear = Wire(UInt((wSpec.manW+1).W)) // has leading 1 to check rounding
  val wIncNear = Wire(Bool())
  if (wSpec.manW < zNearFracWidth) {
    val lsbNear    = xyzNearNorm(zNearWidth-1-wSpec.manW)
    val roundNear  = xyzNearNorm(zNearWidth-1-wSpec.manW-1)
    val stickyNear = xyzNearNorm(zNearWidth-1-wSpec.manW-2, 0).orR.asBool

    wIncNear := FloatChiselUtil.roundIncBySpec(roundSpec, lsbNear, roundNear, stickyNear)
    wmanNear := xyzNearNorm(zNearWidth-1, zNearWidth-1-wSpec.manW);
  } else {
    // w is wider, no rounding required
    wIncNear := false.asBool
    wmanNear := padLSB(wSpec.manW+1, xyzNearNorm)
  }
  dbgPrintf("wmanNear     = 0b%b(%d)\n", wmanNear, wmanNear)

  //----------------------------------------------------------------------
  // Merge Far(Prod/Addend) / Near

  val wman0 = wmanNear + wIncNear.asUInt;
  dbgPrintf("wman0     = 0b%b(%d)\n", wman0, wman0)
  // TODO shift and round once more, considering carry by this

  val wexMax = xyMaxEx + 2
  val wexMin = xyMinEx - zNearWidth
  val wexW   = getWidthToRepresentNumbersAlwaysSigned(Seq(wexMin, wexMax))
  val wex0   = xyEx0 - xyzNearExDec + wSpec.exBias.S

  dbgPrintf("xyEx0 =  %d\n", xyEx0)
  dbgPrintf("ExDec = -%d\n", xyzNearExDec)
  dbgPrintf("wex0  =  %d\n", wex0)

  val wZeroAfterAdd = (wex0.asSInt < 0.S)
  val wInfAfterAdd  = wex0 >= maskI(wSpec.exW).S

  val wsgn0 = xyzNearSign

  //----------------------------------------------------------------------
  // Final
  val infOrNaN = wInfAfterAdd  || xyzinf || xyznan
  val wzero    = wZeroAfterAdd || (xyzero && zzero)

  val wex  = Mux(infOrNaN, maskU(wSpec.exW), Mux(wzero, 0.U(wSpec.exW), wex0(wSpec.exW-1,0)))
  val wman = Mux(infOrNaN || wzero, xyznan ## 0.U((wSpec.manW-1).W), wman0(wSpec.manW-1, 0))
  val wsgn = Mux(xyzinf, xysgn, wsgn0)
  // +inf + inf = +inf, -inf - inf = inf, but inf - inf = nan.
  // so if xyzinf == true, then xysgn and zsgn should be the same.
  // if infAfterAdd == true, then xyzinf != true.

  dbgPrintf("w0    = %d|%d|%d\n", wsgn, wex, wman)

  val w0 = wsgn ## wex ## wman

  io.w   := ShiftRegister(w0, nStage)
  dbgPrintf("x=%x y=%x z=%x, w=%x\n", io.x, io.y, io.z, io.w)
}

class FusedMulAddFP64( stage : PipelineStageConfig )
    extends FusedMulAddFPGeneric( RealSpec.Float64Spec,
      RealSpec.Float64Spec,  RealSpec.Float64Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, stage) {
}

object FusedMulAddFP64_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new FusedMulAddFP64(PipelineStageConfig.none)) ) )
}
