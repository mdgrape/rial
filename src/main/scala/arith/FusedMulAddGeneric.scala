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

  val (xsgn, xex, xman0) = FloatChiselUtil.decomposeWithLeading1(xSpec, io.x)
  val (ysgn, yex, yman0) = FloatChiselUtil.decomposeWithLeading1(ySpec, io.y)
  val (zsgn, zex, zman0) = FloatChiselUtil.decomposeWithLeading1(zSpec, io.z)

  dbgPrintf("x   = %d|%d(%d)|%d\n", xsgn, xex, xex.asSInt-xSpec.exBias.S, xman0)
  dbgPrintf("y   = %d|%d(%d)|%d\n", ysgn, yex, yex.asSInt-ySpec.exBias.S, yman0)
  dbgPrintf("z   = %d|%d(%d)|%d\n", zsgn, zex, zex.asSInt-zSpec.exBias.S, zman0)

  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, io.x)
  val (yzero, yinf, ynan) = FloatChiselUtil.checkValue(ySpec, io.y)
  val (zzero, zinf, znan) = FloatChiselUtil.checkValue(zSpec, io.z)

  // if subnormal numbers are disabled, the mantissa of "zero(exponent == 0)"
  // may contain non-zero value and it can affect to the result. to avoid it,
  // we replace the whole mantissa of z with zeroes.
  val xman = Mux(xzero, 0.U, xman0)
  val yman = Mux(yzero, 0.U, yman0)
  val zman = Mux(zzero, 0.U, zman0)

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

  val diffMaxXYMinZ = xyMaxEx - zSpec.exMin
  val diffMaxZMinXY = zSpec.exMax - xyMinEx
  val diffExW       = log2Up(max(diffMaxXYMinZ, diffMaxZMinXY) + 1) + 1

  dbgPrintf(f"xyMaxEx  = $xyMaxEx%d = xMaxEx(${xSpec.exMax}%d) + yMaxEx(${ySpec.exMax}%d)\n")
  dbgPrintf(f"xyMinEx  = $xyMinEx%d = xMaxEx(${xSpec.exMin}%d) + yMaxEx(${ySpec.exMin}%d)\n")
  dbgPrintf(f"xyExW0   = $xyExW0%d\n")
  dbgPrintf(f"xyExW    = $xyExW%d\n")
  dbgPrintf(f"xyExBias = $xyExBias%d\n")

  val xyExNobias = xex.pad(xyExW).asSInt + yex.pad(xyExW).asSInt - xyExBias.S(xyExW.W)
  val zExNobias  = zex.pad(diffExW).asSInt - zSpec.exBias.S(diffExW.W)

  val zEx0  = Mux(zzero,  xyExNobias, zExNobias)
  val xyEx0 = Mux(xyzero, zExNobias, xyExNobias)

  dbgPrintf("xex   = %d\n", xex.pad(xyExW))
  dbgPrintf("yex   = %d\n", yex.pad(xyExW))
  dbgPrintf("zex   = %d\n", zex)
  dbgPrintf("xyEx0 = %d\n", xyEx0)
  dbgPrintf("zexnb = %d\n", zEx0)

  dbgPrintf(f"diffMaxXYMinZ = $diffMaxXYMinZ%d\n")
  dbgPrintf(f"diffMaxZMinXY = $diffMaxZMinXY%d\n")
  dbgPrintf(f"diffExW       = $diffExW%d\n")

  val diffExXYZ = xyEx0.pad(diffExW).asSInt - zEx0.pad(diffExW).asSInt
  val diffExZXY = zEx0.pad(diffExW).asSInt - xyEx0.pad(diffExW).asSInt

  // later we may take abs of diffEx, so we calculate + and - beforehand.
  dbgPrintf("diffExXYZ = %d\n", diffExXYZ)
  dbgPrintf("diffExZXY = %d\n", diffExZXY)

  // zEx - xyEx == +2, +1, 0, -1 and signs are different
  val isNear = (!xyzero) && (!zzero) && (xysgn =/= zsgn) &&
               ((diffExZXY === 2.S) || (diffExZXY ===  1.S) ||
                (diffExZXY === 0.S) || (diffExZXY === -1.S))
  val exXYLargerThanZ = diffExXYZ(diffExW-1) === 0.U // exponent of xy >= z

  // xyEx - zEx >= 0 (if 0 or 1, signs are the same)
  val isFarProd   = (!isNear) && exXYLargerThanZ
  val isFarAddend = (!isNear) && (!isFarProd)    // zEx - xyEx >= 1

  // -------------------------------------------------------------------------
  // multiply

  val xyprod = xman * yman;

  // -------------------------------------------------------------------------
  // far path (product is larger)
  //
  // xyEx > zEx, so |deltaEx| = xyEx - zEx = diffExXYZ >= 2.
  // if diffExXYZ == 1 || 0 and signs are different, use near path.
  // Note that, if diffExXYZ == 1 || 0 and signs are the same, use this path.
  //
  // +--carry
  // |      +---- 1+xSpec.manW + 1+ySpec.manW = 48
  // |      |
  // v.-----+-----.
  // 0**.**********00         |   001.0000000
  //     '---+---'^           |  -  0.00000zzzz
  //         |    |           |     0.11111wwww
  //         |    + round     |     ^
  //         + wSpec.manW+1   |     ex -= 1

  val diffExFarProd = diffExXYZ.asUInt
  val farProdWidth = max(2 + (1+xSpec.manW) + (1+ySpec.manW)+1,
                         3 + (1+wSpec.manW) + 1)
  val farProdFracWidth = farProdWidth - 4

  when(isFarProd) {
    dbgPrintf("diffExFarProd     = %d\n", diffExFarProd)
    dbgPrintf(f"farProdWidth      = $farProdWidth%d\n")
    dbgPrintf(f"farProdFracWidth  = $farProdFracWidth%d\n")
  }

  //      +-- 1+xSpec.manW(=3) + 1+ySpec.manW(=3)
  //  .---+---.|
  // 0**.******|
  // 00z.zzzz  |(no shift)
  //     zzzzz |(>>1)
  //           |zzzzz (>>7 = xyFracWidth + 1)
  //
  val zShiftOutLimit = farProdFracWidth + 1
  val zShiftOutW     = log2Up(zShiftOutLimit)
  val zShiftOut      = (diffExFarProd >= zShiftOutLimit.U)

  when(isFarProd) {
    dbgPrintf(f"zShiftOutLimit  = $zShiftOutLimit%d\n")
    dbgPrintf(f"zShiftOutW      = $zShiftOutW%d\n")
    dbgPrintf("zShiftOut       = %d\n", zShiftOut)
  }

  val zManPad       = 0.U(3.W) ## zman ## 0.U((farProdWidth - (1+zSpec.manW) - 3).W)
  val zManInverted  = Mux(xyzop, -(zManPad.asSInt), zManPad.asSInt)
  val zManAligned0  = zManInverted.asSInt >> diffExXYZ(zShiftOutW-1, 0)
  val zManAligned   = Mux(zShiftOut, 0.U,
    zManAligned0(zManAligned0.getWidth - 1, zManAligned0.getWidth - 1 - (farProdWidth-1)))
  val prodPad       = (0.U(2.W) ## xyprod ## 0.U(farProdWidth - 2 - xyprod.getWidth))
  val sumFarProd    = prodPad + zManAligned
  val shiftSticky   = PriorityEncoder(zManInverted) < diffExXYZ(zShiftOutW-1, 0) || zShiftOut

  when(isFarProd) {
    dbgPrintf("zManPad       = 0b%b(%d), width= %d\n", zManPad     , zManPad     , zManPad     .getWidth.U)
    dbgPrintf("zManInverted  = 0b%b(%d), width= %d\n", zManInverted, zManInverted, zManInverted.getWidth.U)
    dbgPrintf("zManAligned0  = 0b%b(%d), width= %d\n", zManAligned0, zManAligned0, zManAligned0.getWidth.U)
    dbgPrintf("xyProd        = 0b  %b(%d), width=  %d\n", xyprod     , xyprod     , xyprod     .getWidth.U)
    dbgPrintf("zManAligned   = 0b%b(%d), width= %d\n", zManAligned , zManAligned , zManAligned .getWidth.U)
    dbgPrintf("ProdPad       = 0b%b(%d), width= %d\n", prodPad     , prodPad     , prodPad     .getWidth.U)
    dbgPrintf("sumFarProd    = 0b%b(%d), width= %d\n", sumFarProd  , sumFarProd  , sumFarProd  .getWidth.U)
  }

  val farProdSgn    = Mux(xysgn.asBool, ~sumFarProd(sumFarProd.getWidth-1), sumFarProd(sumFarProd.getWidth-1))
  val farProdAbs    = Mux(sumFarProd(sumFarProd.getWidth-1), -sumFarProd, sumFarProd)
  val shiftFarProd  = PriorityEncoder(Reverse(farProdAbs))

  // xxxx.xxxx   : shiftNearSum 0, ex += 3
  // 0xxx.xxxxx  : shiftNearSum 1, ex += 2
  // 00xx.xxxxxx : shiftNearSum 2, ex += 1
  // 000x.xxxxxxx: shiftNearSum 3, ex += 0

  val farProdNorm   = (farProdAbs << shiftFarProd)(farProdWidth-1, 0)
  val farProdExInc  = 3.S - shiftFarProd.asSInt

  when(isFarProd) {
    dbgPrintf("farProdSgn    = 0b%b(%d), width= %d\n", farProdSgn  , farProdSgn  , farProdSgn  .getWidth.U)
    dbgPrintf("farProdAbs    = 0b%b(%d), width= %d\n", farProdAbs  , farProdAbs  , farProdAbs  .getWidth.U)
    dbgPrintf("shiftFarProd  = 0b%b(%d), width= %d\n", shiftFarProd, shiftFarProd, shiftFarProd.getWidth.U)
    dbgPrintf("farProdNorm   = 0b%b(%d), width= %d\n", farProdNorm , farProdNorm , farProdNorm .getWidth.U)
    dbgPrintf("farProdExInc  = 0b%b(%d), width= %d\n", farProdExInc, farProdExInc, farProdExInc.getWidth.U)
  }

  val wman0FarProd = Wire(UInt(wSpec.manW.W))
  val wex0FarProd  = Wire(SInt((xyExW+1).W))
  if (wSpec.manW < farProdFracWidth) {
    //                                             + leading 1
    //                                + MSB        |   + mantissa width
    //                              .-+----------. v .-+--------.
    val lsbFarProd    = farProdNorm(farProdWidth-1-1-wSpec.manW+1)
    val roundFarProd  = farProdNorm(farProdWidth-1-1-wSpec.manW)
    val stickyFarProd = farProdNorm(farProdWidth-1-1-wSpec.manW-1, 0).orR.asBool || shiftSticky
    when(isFarProd) {
      dbgPrintf("lsbFarProd    = 0b%b\n", lsbFarProd   )
      dbgPrintf("roundFarProd  = 0b%b\n", roundFarProd )
      dbgPrintf("stickyFarProd = 0b%b\n", stickyFarProd)
    }
    val wIncFarProd = FloatChiselUtil.roundIncBySpec(roundSpec, lsbFarProd, roundFarProd, stickyFarProd)
    val wmanFarProd = farProdNorm(farProdWidth-1, farProdWidth-1-wSpec.manW); // ?

    // 1.11111111 + Round = 10.00000000 -> ex+1 1.000000
    val carryRoundFarProd  = wIncFarProd && (wmanFarProd === maskL(wSpec.manW+1).U)

    wman0FarProd := Mux(carryRoundFarProd, 0.U, (wmanFarProd + wIncFarProd.asUInt)(wSpec.manW-1, 0))
    wex0FarProd  := Mux(carryRoundFarProd, xyEx0 + farProdExInc + (wSpec.exBias + 1).S((xyExW+1).W),
                                           xyEx0 + farProdExInc +  wSpec.exBias.S((xyExW+1).W))

    when(isFarProd) {
      dbgPrintf("wIncFarProd   = 0b%b(%d)\n", wIncFarProd, wIncFarProd )
      dbgPrintf("wmanFarProd   = 0b%b(%d)\n", wmanFarProd, wmanFarProd)
    }
  } else {
    // w is wider, no rounding required
    wman0FarProd := padLSB(wSpec.manW+1, farProdNorm)
    wex0FarProd  := xyEx0 + farProdExInc + wSpec.exBias.S((xyExW+1).W)
  }

  when(isFarProd) {
    dbgPrintf("wman0FarProd  = 0b%b(%d), width= %d\n", wman0FarProd  , wman0FarProd  , wman0FarProd.getWidth.U)
    dbgPrintf("wex0FarProd   = 0b%b(%d), width= %d\n", wex0FarProd   , wex0FarProd   , wex0FarProd .getWidth.U)
    dbgPrintf("xyEx0         = 0b%b(%d), width= %d\n", xyEx0         , xyEx0         , xyEx0       .getWidth.U)
    dbgPrintf("farProdExInc  = 0b%b(%d), width= %d\n", farProdExInc  , farProdExInc  , farProdExInc.getWidth.U)
  }

  // -------------------------------------------------------------------------
  // far path (addend is larger)
  //
  // xyEx < zEx, so deltaEx = zEx - xyEx = diffExZXY >= 3.
  // if deltaEx == 2 || 1 || 0 || -1, use near path.
  //
  // 00*.******
  // 000.000xxxxxxxxxxxxxxxx

  val diffExFarAdd = diffExZXY.asUInt // >= 3
  val farAddWidth = max(1 + (1+xSpec.manW) + (1+ySpec.manW), 2 + (1+wSpec.manW) + 1)
  val farAddFracWidth = farAddWidth - 2

  when(isFarAddend) {
    dbgPrintf("diffExFarAdd   = %d\n", diffExFarAdd)
  }

  val xyShiftOutLimit = farAddFracWidth + 2
  val xyShiftOutW     = log2Up(xyShiftOutLimit)
  val xyShiftOut      = (diffExFarAdd >= xyShiftOutLimit.U)

  val zManFarAddPad      = 0.U(2.W) ## zman ## 0.U((farAddWidth - (1+zSpec.manW) - 2).W)
  val zManFarAddInverted = Mux(xyzop, -(zManFarAddPad.asSInt), zManFarAddPad.asSInt)
  val xyManFarAddAligned = 0.U(1.W) ## (xyprod >> diffExFarAdd(xyShiftOutW-1, 0))
  val sumFarAdd          = Mux(xyShiftOut, zManFarAddInverted, xyManFarAddAligned.asSInt + zManFarAddInverted)
  val shiftStickyFarAdd  = PriorityEncoder(zManInverted) < diffExXYZ(xyShiftOutW-1, 0) || xyShiftOut

  when(isFarAddend) {
    dbgPrintf("xyProd              = 0b%b(%d), width= %d\n", xyprod     , xyprod     , xyprod     .getWidth.U)
    dbgPrintf("xyManFarAddAligned  = 0b%b(%d), width= %d\n", xyManFarAddAligned, xyManFarAddAligned, xyManFarAddAligned.getWidth.U)
    dbgPrintf("zManFarAddPad       = 0b%b(%d), width= %d\n", zManFarAddPad     , zManFarAddPad     , zManFarAddPad     .getWidth.U)
    dbgPrintf("zManFarAddInverted  = 0b%b(%d), width= %d\n", zManFarAddInverted, zManFarAddInverted, zManFarAddInverted.getWidth.U)
    dbgPrintf("sumFarAdd           = 0b%b(%d), width= %d\n", sumFarAdd         , sumFarAdd         , sumFarAdd         .getWidth.U)
    dbgPrintf("shiftStickyFarAdd   = 0b%b(%d), width= %d\n", shiftStickyFarAdd , shiftStickyFarAdd , shiftStickyFarAdd .getWidth.U)
  }

  val farAddSgn   = Mux(xysgn.asBool, ~sumFarAdd(sumFarAdd.getWidth-1), sumFarAdd(sumFarAdd.getWidth-1))
  val farAddAbs   = Mux(sumFarAdd(sumFarAdd.getWidth-1), -sumFarAdd, sumFarAdd)
  val shiftFarAdd = PriorityEncoder(Reverse(farAddAbs.asUInt))

  // 00*.****** (zex)
  // 000.000xxxxxxxxxxxxxxxx
  // 00w.wwwww
  // www.www  <<2 -> ex+=0
  //
  // 000.wwwwww
  // www.www  <<3 -> ex-=1
  //
  // 000.0wwwwww
  // www.www  <<4 -> ex-=2

  val farAddNorm   = (farAddAbs << shiftFarAdd)(farAddWidth-1, 0)
  val farAddExDec  = shiftFarAdd - 2.U

  when(isFarAddend) {
    dbgPrintf("shiftFarAdd   = 0b%b(%d), width= %d\n", shiftFarAdd  , shiftFarAdd , shiftFarAdd .getWidth.U)
    dbgPrintf("farAddNorm    = 0b%b(%d), width= %d\n", farAddNorm   , farAddNorm  , farAddNorm  .getWidth.U)
    dbgPrintf("farAddExDec   = 0b%b(%d), width= %d\n", farAddExDec  , farAddExDec , farAddExDec .getWidth.U)
  }

  val wman0FarAdd = Wire(UInt(wSpec.manW.W))
  val wex0FarAdd  = Wire(SInt((xyExW+1).W))
  if (wSpec.manW < farAddFracWidth) {
    //                                          + leading 1
    //                              + MSB       |   + mantissa width
    //                            .-+---------. v .-+--------.
    val lsbFarAdd    = farAddNorm(farAddWidth-1-1-wSpec.manW+1)
    val roundFarAdd  = farAddNorm(farAddWidth-1-1-wSpec.manW)
    val stickyFarAdd = farAddNorm(farAddWidth-1-1-wSpec.manW-1, 0).orR.asBool || shiftStickyFarAdd
    when(isFarAddend) {
      dbgPrintf("lsbFarAdd    = 0b%b\n", lsbFarAdd   )
      dbgPrintf("roundFarAdd  = 0b%b\n", roundFarAdd )
      dbgPrintf("stickyFarAdd = 0b%b\n", stickyFarAdd)
    }
    val wIncFarAdd = FloatChiselUtil.roundIncBySpec(roundSpec, lsbFarAdd, roundFarAdd, stickyFarAdd)
    val wmanFarAdd = farAddNorm(farAddWidth-1, farAddWidth-1-wSpec.manW);

    // 1.11111111 + Round = 10.00000000 -> ex+1 1.000000
    val carryRoundFarAdd = wIncFarAdd && (wmanFarAdd === maskL(wSpec.manW+1).U)

    wman0FarAdd := Mux(carryRoundFarAdd, 0.U, (wmanFarAdd + wIncFarAdd.asUInt)(wSpec.manW-1, 0))
    wex0FarAdd  := Mux(carryRoundFarAdd, zExNobias.asSInt - farAddExDec.asSInt + (wSpec.exBias+1).S,
                                         zExNobias.asSInt - farAddExDec.asSInt +  wSpec.exBias.S)
  } else {
    // w is wider, no rounding required
    wman0FarAdd := padLSB(wSpec.manW+1, farAddNorm)
    wex0FarAdd  := zExNobias.asSInt - farAddExDec.asSInt +  wSpec.exBias.S
  }

  when(isFarAddend) {
    dbgPrintf("zex         = 0b%b(%d), width = %d\n", zex        , zex        , zex        .getWidth.U)
    dbgPrintf("farAddExDec = 0b%b(%d), width = %d\n", farAddExDec, farAddExDec, farAddExDec.getWidth.U)
    dbgPrintf("wex0FarAdd  = 0b%b(%d), width = %d\n", wex0FarAdd , wex0FarAdd , wex0FarAdd .getWidth.U)
  }

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
  // if 2+x.manW+y.manW < w.manW
  //
  // 0011.xxxxxxxxx000
  // 000z.zzz000000000
  //      '--w.manW--'
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
  val nearSum        = (0.U(2.W) ## xyprod) - zNearShift
  val nearSgn        = Mux(xysgn.asBool, ~nearSum(zNearWidth-1), nearSum(zNearWidth-1))
  val nearAbs        = Mux(nearSum(zNearWidth-1), -nearSum, nearSum)
  val shiftNearSum   = PriorityEncoder(Reverse(nearAbs))
  when(isNear) {
    dbgPrintf(f"zNearWidth     = $zNearWidth%d\n")
    dbgPrintf(f"zNearFracWidth = $zNearFracWidth%d\n")
    dbgPrintf("zNearPad       = 0b%b(%d), width= %d\n",   zNearPad   , zNearPad   , zNearPad   .getWidth.U)
    dbgPrintf("zNearShift     = 0b%b(%d), width= %d\n",   zNearShift , zNearShift , zNearShift .getWidth.U)
    dbgPrintf("xyProd         = 0b  %b(%d), width= %d\n", xyprod     , xyprod     , xyprod     .getWidth.U)
    dbgPrintf("nearSum        = 0b%b(%d), width= %d\n",   nearSum    , nearSum    , nearSum    .getWidth.U)
    dbgPrintf("nearSgn        = 0b%b(%d), width=%d\n",    nearSgn    , nearSgn    , nearSgn    .getWidth.U)
    dbgPrintf("nearAbs        = 0b%b(%d), width=%d\n",    nearAbs    , nearAbs    , nearAbs    .getWidth.U)
    dbgPrintf("shiftNearSum   = %b, width= %d\n", shiftNearSum , shiftNearSum.getWidth.U)
  }

  val nearNorm = (nearAbs << shiftNearSum)(zNearWidth-1, 0)
  val nearExDec = shiftNearSum.asSInt - 3.S
  when(isNear) {
    dbgPrintf("nearNorm  = %d, width=%d\n", nearNorm , nearNorm .getWidth.U)
    dbgPrintf("nearExDec = %d, width=%d\n", nearExDec, nearExDec.getWidth.U)
  }

  // 01xx.xxxxx   : shiftNearSum 1, ex += 2 = 3 - shiftNearSum
  // 001x.xxxxxx  : shiftNearSum 2, ex += 1
  // 0001.xxxxxxx : shiftNearSum 3, ex += 0
  // 0000.1xxxxxxx: shiftNearSum 4, ex += -1
  // |
  // v
  // 1xxx xxxxxxxx

  val wman0Near = Wire(UInt(wSpec.manW.W))
  val wex0Near  = Wire(SInt((xyExW+1).W))
  if (wSpec.manW < zNearFracWidth) {
    val lsbNear    = nearNorm(zNearWidth-1-wSpec.manW)
    val roundNear  = nearNorm(zNearWidth-1-wSpec.manW-1)
    val stickyNear = nearNorm(zNearWidth-1-wSpec.manW-2, 0).orR.asBool

    val wIncNear = FloatChiselUtil.roundIncBySpec(roundSpec, lsbNear, roundNear, stickyNear)
    val wmanNear = nearNorm(zNearWidth-1, zNearWidth-1-wSpec.manW); // ?

    val carryRoundNear = wIncNear && (wmanNear === maskL(wSpec.manW+1).U)

    wman0Near := Mux(carryRoundNear, 0.U, (wmanNear + wIncNear.asUInt)(wSpec.manW-1, 0))
    wex0Near  := Mux(carryRoundNear, xyEx0 - nearExDec + (wSpec.exBias+1).S,
                                     xyEx0 - nearExDec + wSpec.exBias.S)
    when(isNear) {
      dbgPrintf("wmanNear  = 0b%b(%d)\n", wmanNear, wmanNear)
      dbgPrintf("wman0Near = 0b%b(%d)\n", wmanNear, wmanNear)
    }
  } else {
    // w is wider, no rounding required
    wman0Near := padLSB(wSpec.manW+1, nearNorm)
    wex0Near  := xyEx0 - nearExDec + wSpec.exBias.S
  }

  when(isNear) {
    dbgPrintf("wman0Near = 0b%b(%d)\n", wman0Near, wman0Near)
    dbgPrintf("xyEx0     =  %d\n", xyEx0)
    dbgPrintf("NearExDec = -%d\n", nearExDec)
    dbgPrintf("w.exBias  =  %d\n", wSpec.exBias.U)
    dbgPrintf("wex0Near  = 0b%b(%d)\n", wex0Near,  wex0Near)
  }

  //----------------------------------------------------------------------
  // Merge Far(Prod/Addend) / Near

  val wsgn0 = Mux(isNear, nearSgn,   Mux(isFarProd, farProdSgn  , farAddSgn  ))
  val wman0 = Mux(isNear, wman0Near, Mux(isFarProd, wman0FarProd, wman0FarAdd))
  val wex0  = Mux(isNear, wex0Near,  Mux(isFarProd, wex0FarProd , wex0FarAdd ))

  val wZeroAfterAdd = (wex0.asSInt < 0.S)
  val wInfAfterAdd  = wex0 >= maskI(wSpec.exW).S

  //----------------------------------------------------------------------
  // Final
  val infOrNaN = wInfAfterAdd  || xyzinf || xyznan
  val wzero    = wZeroAfterAdd || (xyzero && zzero)

  val wex  = Mux(infOrNaN, maskU(wSpec.exW), Mux(wzero, 0.U(wSpec.exW), wex0(wSpec.exW-1,0)))
  val wman = Mux(infOrNaN || wzero, xyznan ## 0.U((wSpec.manW-1).W), wman0(wSpec.manW-1, 0))
  val wsgn = Mux(xyznan, 0.U, Mux(xyzinf, xysgn, wsgn0)) // XXX sign of NaN?
  // +inf + inf = +inf, -inf - inf = inf, but inf - inf = nan.
  // so if xyzinf == true, then xysgn and zsgn should be the same.
  // if infAfterAdd == true, then xyzinf != true.

  dbgPrintf("wsgn  = %b(%d), width = %d\n", wsgn, wsgn, wsgn.getWidth.U)
  dbgPrintf("wex   = %b(%d), width = %d\n", wex,  wex,  wex .getWidth.U)
  dbgPrintf("wman  = %b(%d), width = %d\n", wman, wman, wman.getWidth.U)
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
