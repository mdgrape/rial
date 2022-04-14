//% @file 3opAddFPGeneric.scala
//
// Generic 3-op Floating-Point Adder Unit
// Copyright (C) Toru Niina RIKEN BDR 2022
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

class ThreeOpAddFPGeneric(
  xSpec: RealSpec, ySpec: RealSpec, zSpec: RealSpec, wSpec: RealSpec, // Input / Output floating spec
  roundSpec: RoundSpec, // Rounding spec
  stage: PipelineStageConfig,
  val enableDebug: Boolean = false
) extends Module with DebugControlSlave {

  val nStage = stage.total

  def getParam = { (xSpec, ySpec, zSpec, wSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x   = Input(UInt(xSpec.W.W))
    val y   = Input(UInt(ySpec.W.W))
    val z   = Input(UInt(zSpec.W.W))
    val w   = Output(UInt(wSpec.W.W))
  })

  val (xsgn, xex, xmanW1) = FloatChiselUtil.decomposeWithLeading1(xSpec, io.x)
  val (ysgn, yex, ymanW1) = FloatChiselUtil.decomposeWithLeading1(ySpec, io.y)
  val (zsgn, zex, zmanW1) = FloatChiselUtil.decomposeWithLeading1(zSpec, io.z)

  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, io.x)
  val (yzero, yinf, ynan) = FloatChiselUtil.checkValue(ySpec, io.y)
  val (zzero, zinf, znan) = FloatChiselUtil.checkValue(zSpec, io.z)

  val wnan = xnan || ynan || znan ||
            (xinf && yinf && (xsgn =/= ysgn)) ||
            (yinf && zinf && (ysgn =/= zsgn)) ||
            (zinf && xinf && (zsgn =/= xsgn))
  val winf0 = (xinf || yinf || zinf) && !wnan
  val winf0sgn = ! ((xinf && !xsgn) || (yinf && !ysgn) || (zinf && !zsgn))
  val wzero0 = xzero && yzero && zzero
  val wzero0sgn = if(roundSpec == RoundSpec.truncate) {
    xsgn.asBool || ysgn.asBool || zsgn.asBool
  } else {
    xsgn.asBool && ysgn.asBool && zsgn.asBool
  }

  //----------------------------------------------------------------------
  // Exponent comparison

  val xexBias = xSpec.exBias // just an alias
  val yexBias = ySpec.exBias
  val zexBias = zSpec.exBias

  def calcdEx(aexBias: Int, aex: UInt, bexBias: Int, bex: UInt) = {
    if(aexBias == bexBias) {
      (aex -& bex, bex -& aex)
    } else if(aexBias > bexBias) {
      (aex -& bex - (aexBias - bexBias).U, bex -& aex + (aexBias - bexBias).U)
    } else {
      (aex -& bex + (bexBias - aexBias).U, bex -& aex - (bexBias - aexBias).U)
    }
  }

  val (dexXY, dexYX) = calcdEx(xexBias, xex, yexBias, yex)
  val (dexYZ, dexZY) = calcdEx(yexBias, yex, zexBias, zex)
  val (dexZX, dexXZ) = calcdEx(zexBias, zex, xexBias, xex)

  val fracW = Seq(xSpec.manW, ySpec.manW, zSpec.manW, wSpec.manW).max

  val xleqy = dexYX.head(1).asBool || !dexYX.orR // x >= y <=> 0 >= y - x
  val yleqz = dexZY.head(1).asBool || !dexZY.orR
  val zleqx = dexXZ.head(1).asBool || !dexXZ.orR
  assert(xleqy || yleqz || zleqx)

  val ovfBits = 3 // bits to protect from overflow
  val intermedW = ovfBits + 1 + 2*fracW + 3

  def shiftAndCollectSticky(spec: RealSpec, manW1: UInt, shf: UInt, isZero: Bool) = {
    val padded   = Cat(Seq(0.U(ovfBits.W), manW1, 0.U((intermedW - (ovfBits+1+spec.manW) - 1).W),
                           0.U(manW1.getWidth.W)))
    val shiftOut = shf >= log2Up(intermedW-ovfBits).U
    val shifted  = padded >> shf
    val sticky   = shifted(manW1.getWidth-1, 0).orR || shiftOut
    val retval   = Mux(isZero, 0.U, shifted(shifted.getWidth-1, manW1.getWidth))
    Cat(retval, sticky)
  }

  val shiftXY = Mux(xleqy, shiftAndCollectSticky(ySpec, ymanW1, dexXY, yzero), shiftAndCollectSticky(xSpec, xmanW1, dexYX, xzero))
  val shiftYZ = Mux(yleqz, shiftAndCollectSticky(zSpec, zmanW1, dexYZ, zzero), shiftAndCollectSticky(ySpec, ymanW1, dexZY, yzero))
  val shiftZX = Mux(zleqx, shiftAndCollectSticky(xSpec, xmanW1, dexZX, xzero), shiftAndCollectSticky(zSpec, zmanW1, dexXZ, zzero))
  val xpadded = Cat(Seq(0.U(ovfBits.W), Mux(xzero, 0.U((xSpec.manW+1).W), xmanW1), 0.U((intermedW - (ovfBits+1+xSpec.manW)).W)))
  val ypadded = Cat(Seq(0.U(ovfBits.W), Mux(yzero, 0.U((ySpec.manW+1).W), ymanW1), 0.U((intermedW - (ovfBits+1+ySpec.manW)).W)))
  val zpadded = Cat(Seq(0.U(ovfBits.W), Mux(zzero, 0.U((zSpec.manW+1).W), zmanW1), 0.U((intermedW - (ovfBits+1+zSpec.manW)).W)))
  assert(shiftXY.getWidth == intermedW)
  assert(shiftYZ.getWidth == intermedW)
  assert(shiftZX.getWidth == intermedW)
  assert(xpadded.getWidth == intermedW)
  assert(ypadded.getWidth == intermedW)
  assert(zpadded.getWidth == intermedW)

  val exWidth   = Seq(xSpec.exW, ySpec.exW, zSpec.exW).max + 3
  val xexNobias = xex.zext.pad(exWidth) - xSpec.exBias.S(exWidth.W)
  val yexNobias = yex.zext.pad(exWidth) - ySpec.exBias.S(exWidth.W)
  val zexNobias = zex.zext.pad(exWidth) - zSpec.exBias.S(exWidth.W)

  val exCmp    = Cat(Seq(xleqy, yleqz, zleqx))
  val exMax = MuxLookup(exCmp, xexNobias, Array(
    "b000".U -> 0.S,
    "b001".U -> zexNobias,
    "b010".U -> yexNobias,
    "b011".U -> yexNobias,
    "b100".U -> xexNobias,
    "b101".U -> zexNobias,
    "b110".U -> xexNobias,
    "b111".U -> xexNobias,
  ))
  val xAligned = MuxLookup(exCmp, xpadded, Array(
    "b000".U -> 0.U,
    "b001".U -> shiftZX,
    "b010".U -> shiftXY,
    "b011".U -> shiftXY,
    "b100".U -> xpadded,
    "b101".U -> shiftZX,
    "b110".U -> xpadded,
    "b111".U -> xpadded
    ))
  val yAligned = MuxLookup(exCmp, ypadded, Array(
    "b000".U -> 0.U,
    "b001".U -> shiftYZ,
    "b010".U -> ypadded,
    "b011".U -> ypadded,
    "b100".U -> shiftXY,
    "b101".U -> shiftYZ,
    "b110".U -> shiftXY,
    "b111".U -> ypadded
    ))
  val zAligned = MuxLookup(exCmp, zpadded, Array(
    "b000".U -> 0.U,
    "b001".U -> zpadded,
    "b010".U -> shiftYZ,
    "b011".U -> shiftYZ,
    "b100".U -> shiftZX,
    "b101".U -> zpadded,
    "b110".U -> shiftZX,
    "b111".U -> zpadded
    ))

  val xEffSgn = 0.U
  val yEffSgn = xsgn ^ ysgn
  val zEffSgn = xsgn ^ zsgn
  val xNegSgn = ~xEffSgn
  val yNegSgn = ~yEffSgn
  val zNegSgn = ~zEffSgn
  dbgPrintf("xEffSgn = %b\n", xEffSgn)
  dbgPrintf("yEffSgn = %b\n", yEffSgn)
  dbgPrintf("zEffSgn = %b\n", zEffSgn)

  def genComplLSBs(sgns: UInt) = {
    MuxLookup(sgns, "b000000".U, Array(
      "b000".U -> "b00_00_00".U,
      "b001".U -> "b10_00_10".U,
      "b010".U -> "b00_10_10".U,
      "b011".U -> "b10_11_11".U,
      "b100".U -> "b10_10_00".U,
      "b101".U -> "b11_10_11".U,
      "b110".U -> "b11_11_10".U,
      "b111".U -> "b00_00_00".U
      ))
  }

  val posLSBs = genComplLSBs(Cat(xEffSgn, yEffSgn, zEffSgn))
  val xPosLSBs = posLSBs(5, 4)
  val yPosLSBs = posLSBs(3, 2)
  val zPosLSBs = posLSBs(1, 0)
  val negLSBs = genComplLSBs(Cat(xNegSgn, yNegSgn, zNegSgn))
  val xNegLSBs = negLSBs(5, 4)
  val yNegLSBs = negLSBs(3, 2)
  val zNegLSBs = negLSBs(1, 0)

  val xAlignedPos = Cat(Mux(xEffSgn.asBool, ~xAligned, xAligned), xPosLSBs)
  val yAlignedPos = Cat(Mux(yEffSgn.asBool, ~yAligned, yAligned), yPosLSBs)
  val zAlignedPos = Cat(Mux(zEffSgn.asBool, ~zAligned, zAligned), zPosLSBs)
  dbgPrintf("xAlignedPos = %b\n", xAlignedPos)
  dbgPrintf("yAlignedPos = %b\n", yAlignedPos)
  dbgPrintf("zAlignedPos = %b\n", zAlignedPos)

  val xAlignedNeg = Cat(Mux(xNegSgn.asBool, ~xAligned, xAligned), xNegLSBs)
  val yAlignedNeg = Cat(Mux(yNegSgn.asBool, ~yAligned, yAligned), yNegLSBs)
  val zAlignedNeg = Cat(Mux(zNegSgn.asBool, ~zAligned, zAligned), zNegLSBs)
  dbgPrintf("xAlignedNeg = %b\n", xAlignedNeg)
  dbgPrintf("yAlignedNeg = %b\n", yAlignedNeg)
  dbgPrintf("zAlignedNeg = %b\n", zAlignedNeg)

  // we assume that the aligned mantissas have overflow bit at their MSB
  def csa3to2(x: UInt, y: UInt, z: UInt) = {
    assert(x.getWidth == y.getWidth && x.getWidth == z.getWidth)
    val xyz0 = x ^ y ^ z
    val xyz1 = (x & y) | (y & z) | (z & x)
    (xyz0, Cat(xyz1, 0.U(1.W)).tail(1))
  }

  val (sum0Pos, sum1Pos) = csa3to2(xAlignedPos, yAlignedPos, zAlignedPos)
  val (sum0Neg, sum1Neg) = csa3to2(xAlignedNeg, yAlignedNeg, zAlignedNeg)

  val sumPos = sum0Pos + sum1Pos
  val sumNeg = sum0Neg + sum1Neg
  dbgPrintf("sum0Pos = %b\n", sum0Pos)
  dbgPrintf("sum1Pos = %b\n", sum1Pos)
  dbgPrintf("sumPos  = %b\n", sumPos)
  dbgPrintf("sum0Neg = %b\n", sum0Neg)
  dbgPrintf("sum1Neg = %b\n", sum1Neg)
  dbgPrintf("sumNeg  = %b\n", sumNeg)

  val useNeg = sumPos.head(1).asBool
  val sum = Mux(useNeg, sumNeg.asUInt, sumPos.asUInt)

  val sumShift = PriorityEncoder(Reverse(sum))
  val sumShifted  = (sum << sumShift)(sum.getWidth-1, 0)
  dbgPrintf("sumShift = %d\n", sumShift)
  assert(sumShifted(sumShifted.getWidth-1) === 1.U || wzero0)

  val sumRoundBit = sum.getWidth-1 - (1+wSpec.manW)
  val sumLSB      = sumShifted(sumRoundBit+1)
  val sumRound    = sumShifted(sumRoundBit)
  val sumSticky   = sumShifted(sumRoundBit-1, 0).orR
  dbgPrintf("sumShift = %d\n", sumShift - 2.U)

  val sumRoundInc = FloatChiselUtil.roundIncBySpec(roundSpec, sumLSB, sumRound, sumSticky)

  val wman0 = sumShifted(sum.getWidth-1, sum.getWidth - (1+wSpec.manW)) +& sumRoundInc
  val sumMoreThan2AfterRound = wman0(wSpec.manW+1)
  val wexNobias = exMax -& sumShift.zext + ovfBits.S + sumMoreThan2AfterRound.zext
  val wex0  = (wexNobias + wSpec.exBias.S)(wSpec.exW-1, 0)
  dbgPrintf("exMax       = %d\n", exMax)
  dbgPrintf("wexNobias   = %d\n", wexNobias)

  val wzero = (wexNobias < wSpec.exMin.S) || wzero0
  val winf  = (wexNobias > wSpec.exMax.S) || winf0

  val wsgn = Mux(wnan, 0.U,
             Mux(winf0, winf0sgn,
             Mux(wzero0, wzero0sgn,
             Mux(useNeg, !xsgn, xsgn))))

  val wex  = Mux(wzero, 0.U(wSpec.exW.W),
             Mux(winf || wnan,  maskI(wSpec.exW).U(wSpec.exW.W), wex0))
  val wman = Mux(wzero || winf || wnan, Cat(wnan, 0.U((wSpec.manW-1).W)),
                                        wman0(wSpec.manW-1, 0))

  val w0 = Cat(wsgn, wex, wman)
  assert(w0.getWidth == wSpec.W)

  io.w   := ShiftRegister(w0, nStage)
}
