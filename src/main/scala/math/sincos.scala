//% @file sincos.scala
//
// x -> sin(x) or cos(x)
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import spire.math.SafeLong
import spire.math.Numeric
import spire.math.Real
import spire.implicits._

import chisel3._
import chisel3.util._

import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.FloatChiselUtil
import rial.arith._

import rial.math._

//
// Evaluate sin(x) or cos(x) for the full range of x. Both function shares
// most of the part by converting cos(x) to sin(pi/2-x).
// So only SinCosPreProcess has conditional branch that depends on `isSin:Bool`.
//
// XXX Basically, the result depends on the fractional part of x / pi. If x is
// too large compared to the constant `Pi` used here, the fractional part will
// be trancated and the resulting value will lose the precision. Currently, we
// tested that the current implementation has enough precision (err < 2) with
// |x| <= 64pi. Until |x| < 128pi, the error increases a bit, but the result
// looks still makes sense (err < the last 3bit). In case of 128pi < |x|,
// the error increases (like, the last 8 bits) and is intolerable.
//

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

object SinCosPreMulArgs {
  def lhsW(spec: RealSpec): Int = {
    1 + spec.manW + Seq(spec.manW+1, 12).max
  }

  def lhs(spec: RealSpec, w: Int): UInt = {
    val diffW = w - lhsW(spec)
    assert(diffW >= 0, f"PreProcMultiplier lhsW(${w}) should be wider or " +
      f"equal to lhsW(${lhsW(spec)})")

    val coef = Real.one / Real.pi
    if(diffW != 0) {
      Cat(coef(lhsW(spec)+1).toBigInt.U(lhsW(spec).W), 0.U(diffW.W))
    } else {
      coef(lhsW(spec)+1).toBigInt.U(lhsW(spec).W)
    }
  }

  def prodW(spec: RealSpec): Int = {
    (1 + spec.manW) + (1 + spec.manW + Seq(spec.manW+1, 12).max)
  }
  def prod(spec: RealSpec, out: UInt): UInt = {
    val outW = out.getWidth
    out(outW-1, outW - prodW(spec))
  }
}

// this contains x/pi rounded into [0, 2) and zsgn.
// zsgn can be calculated from x/pi.
class SinCosPreProcessOutput(val spec: RealSpec) extends Bundle {
  val znan = Output(Bool())
  val zone = Output(Bool())
  val zsgn = Output(UInt(1.W))
  val yex  = Output(UInt(spec.exW.W))
  val yman = Output(UInt(spec.manW.W))
}

class SinCosPreProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val io = IO(new Bundle {
    val en    = Input (UInt(1.W))
    val isSin = Input(Bool())
    val x     = Flipped(new DecomposedRealOutput(spec))

    val adr   = Output(UInt(adrW.W))
    val dx    = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
    val out   = new SinCosPreProcessOutput(spec)
  })

  val xsgn = enable(io.en, io.x.sgn)
  val xex  = enable(io.en, io.x.ex )
  val xman = enable(io.en, io.x.man)

  // --------------------------------------------------------------------------
  // calc x/pi

  val oneOverPiPad = Seq(manW+1, 12).max
  val oneOverPi = (Real.one / Real.pi)(1+manW+1+oneOverPiPad).toBigInt
  val oneOverPiW = 1 + manW + oneOverPiPad

  val xOverPiProd          = Mux(xex === 0.U, 0.U, Cat(1.U(1.W), xman) * oneOverPi.U(oneOverPiW.W))
  val xOverPiProdMoreThan2 = xOverPiProd((1+manW) + oneOverPiW - 1)
  val xOverPiEx            = xex -& 2.U + xOverPiProdMoreThan2
  val xOverPiIsZero        = xOverPiEx(exW)
  // always x.ex > xOverPiEx because xOverPiProdMoreThan2 <= 1. Overflow does
  // not mean inf.

  val xOverPi              = Mux(xOverPiProdMoreThan2, xOverPiProd,
                                 Cat(xOverPiProd.tail(1), 0.U(1.W)))
  val xOverPiFracW         = (1+manW) + oneOverPiW - 1
  assert(xOverPi.getWidth == xOverPiFracW+1) // fraction bits + 1 integer bit
  assert(xOverPi(xOverPiFracW) === 1.U || xOverPiIsZero || io.en === 0.U) // normalized

//   printf("xOverPi      = %d\n", xOverPi)
//   printf("xOverPiFracW = %d\n", xOverPiFracW.U)

  // --------------------------------------------------------------------------
    // convert full range x/pi into (0, 2)

  val xOverPiExMoreThan1 = (xOverPiEx >= exBias.U)
  val xOverPiExNobiasPos = (xOverPiEx - exBias.U)(log2Up(xOverPiFracW)-1, 0)
  val xOverPiExNobiasNeg = (exBias.U - xOverPiEx)(log2Up(xOverPiFracW)-1, 0)

  val xOverPiAlignPos = (xOverPi << xOverPiExNobiasPos)(xOverPiFracW, 0)
  val xOverPiAlignNeg = (xOverPi >> xOverPiExNobiasNeg)(xOverPiFracW, 0)
  val xOverPiAligned = Mux(xOverPiExMoreThan1, xOverPiAlignPos, xOverPiAlignNeg)
//   printf("cir: xOverPiAlignPos = %d\n", xOverPiAlignPos)
//   printf("cir: xOverPiAlignNeg = %d\n", xOverPiAlignNeg)
//   printf("cir: xOverPiAligned  = %d\n", xOverPiAligned )

  // --------------------------------------------------------------------------
  // calculate zsgn from aligned x/pi

  // this determines the range of x.
  // 0 -> [0,0.5), 1 -> [0.5, 1), 2 -> [1, 1.5), 3 -> [1.5, 2)
  val xOverPiAligned2MSBs = xOverPiAligned.head(2)
//   printf("cir:xOverPiAligned2MSBs = %b\n", xOverPiAligned2MSBs)

  // we can already calculate the sign of return value from the value of x/pi.
  //
  // sin(x)                 cos(x)
  // 2MSBs 00 01 10 11      2MSBs 00 01 10 11
  //      |  :  :  :  :          |  :  :  :  :
  //      | .-. :  :  :        .-|-.:  :  :.-:
  // _____|'_:_': _:__:     __.__|__:__:__:__:
  //   _.'|  :  :'._.':          |  :'.:.':  :
  //      |  :  :  :  :          |  :  :  :  :
  //            '-----'             '-----'
  //            negative            negative
  //
  val zsgn = Mux(io.isSin, xOverPiAligned2MSBs(1) ^ xsgn,
                           xOverPiAligned2MSBs(1) ^ xOverPiAligned2MSBs(0)) & io.en
  io.out.zsgn := ShiftRegister(zsgn, nStage)

  // --------------------------------------------------------------------------
  // special value check: if x is inf or nan, then z is nan.

  val znan = xex.andR
  io.out.znan := ShiftRegister(znan, nStage)

  // --------------------------------------------------------------------------
  // calculate x/pi in sin/cos [0, 2) into sin [0, 1/2)

  assert(xOverPiAligned.getWidth == 1 + xOverPiFracW)
  // remove bits that represent larger than 2
  val ymanPos =  xOverPiAligned.tail(2)
  val ymanNeg = ~xOverPiAligned.tail(2)// + 1.U

  // sin(x)                 cos(x)
  // 2MSBs 00 01 10 11      2MSBs 00 01 10 11     MSB(0) | isSin  |   res
  //      |  :  :  :  :          |  :  :  :  :   --------+--------+--------
  //      | .-. :  :  :        .-|-.:  :  :.-:      0    | sin(1) | ymanPos
  // _____|'_:_': _:__:     __.__|__:__:__:__:      1    | sin(1) | ymanNeg
  //   _.'|  :  :'._.':          |  :'.:.':  :      0    | cos(0) | ymanNeg
  //      |  :  :  :  :          |  :  :  :  :      1    | cos(0) | ymanPos
  //         '--'  '--'          '--'  '--'
  //         pi-x                pi/2-x

  val ymanAligned = Mux(xOverPiAligned2MSBs(0) ^ io.isSin, ymanPos, ymanNeg)
//   printf("cir: ymanAligned = %d(%b)\n", ymanAligned , ymanAligned )

  // if 0 < x < pi/2, then we use x/pi itself to avoid
  // loss of precision due to alignment to pi.
  val ymanAsIs = xOverPi(xOverPi.getWidth-1, 2)
  val yexAsIs  = xOverPiEx
  val nonAlign = (!xOverPiExMoreThan1 && xOverPiAligned2MSBs === 0.U && io.isSin)

  assert(ymanAsIs.getWidth == ymanPos.getWidth)
  val yman0    = Mux(nonAlign, ymanAsIs, ymanAligned)
  val yex0     = Mux(nonAlign, yexAsIs, (exBias-2).U(exW.W))
  val ymanIsNonZero = yman0.orR // XXX

//   printf("cir: xOverPi  = %d(%b)\n", xOverPi , xOverPi )
//   printf("cir: ymanAsIs = %d(%b)\n", ymanAsIs, ymanAsIs)
//   printf("cir: ymanPos  = %d(%b)\n", ymanPos , ymanPos )
//   printf("cir: ymanNeg  = %d(%b)\n", ymanNeg , ymanNeg )
//   printf("cir: yman0    = %d\n", yman0   )

  // TODO here we can reduce the area by checking yman0Shift and roundbit

  // yman0 removes its MSB, so here the width becomes xOverPiFracW,
  // not xOverPiFracW+1. But yman0 should be aligned to xOverPiFracW+1.
  // So here we need to add 1.
  val yman0W        = yman0.getWidth
  val yman0Shift0   = PriorityEncoder(Reverse(yman0))
  val yman0ShiftW   = log2Up(yman0W)
  val yman0Shift    = yman0Shift0(yman0ShiftW-1, 0)
  val yman0Shifted  = (yman0 << yman0Shift)(yman0W-1, 0)
  // if x == 0, then postprocess (polynomial * x) makes the result 0.
  assert(yman0Shifted(yman0W-1) === 1.U || xOverPiIsZero || io.en === 0.U)
  val yman0RoundBit = yman0W - manW - 1
  val yman0Rounded  = yman0Shifted(yman0W-2, yman0RoundBit) +&
                      yman0Shifted(yman0RoundBit-1)
  val yman0MoreThan2 = yman0Rounded(manW)

//   printf("cir:yman0Shift0  = %d\n", yman0Shift0)
//   printf("cir:yman0.W      = %d\n", yman0.getWidth.U)
//   printf("cir:yman0ShiftW  = %d\n", yman0ShiftW.U)
//   printf("cir:yman0Shift   = %d\n", yman0Shift)
//   printf("cir:yman0Shifted = %b(%d)\n", yman0Shifted, yman0Shifted)
//   printf("cir:yman0Rounded = %b\n", yman0Rounded)

  val ymanNonZero = yman0Rounded(manW-1, 0)
  val yexNonZero = yex0 - yman0Shift + yman0MoreThan2
  // if it's zero, return zero.
  val yex  = Mux(!ymanIsNonZero, 0.U, yexNonZero )
  val yman = Mux(!ymanIsNonZero, 0.U, ymanNonZero)
//   printf("cir:yex  = %d\n", yex)
//   printf("cir:yman = %b\n", yman)

  io.out.yex  := ShiftRegister(enable(io.en, yex ), nStage)
  io.out.yman := ShiftRegister(enable(io.en, yman), nStage)

  val zone = (yex === (exBias-1).U) && yman === 0.U
  io.out.zone := ShiftRegister(io.en & zone, nStage)

  assert(yex =/= (exBias-1).U || yman === 0.U || io.en === 0.U) // if yex == -1, then yman should be 0

  // -------------------------------------------------------------------------
  // determine address and dx from x/pi aligned

  val yaligned = Cat(1.U(1.W), yman) >> ((exBias-1).U(exW.W) - yex)

  val adr  = enable(io.en, yaligned(manW-1, manW-adrW))
  io.adr   := ShiftRegister(adr,  nStage)

  if(order != 0) {
    val dx   = enable(io.en, Cat(~yaligned(manW-adrW-1), yaligned(manW-adrW-2,0)))
    io.dx.get := ShiftRegister(dx, nStage)
  }
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class SinCosTableCoeff(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
) extends Module {

  val manW   = spec.manW

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input  (UInt(adrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {

    val tableI = SinCosSim.sincosTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    // sign mode 1 = 2's complement and no sign bit
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/1)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    assert(maxCbit(0) == fracW)
    assert(coeff(0).getWidth == fracW)

    io.cs.cs(0) := enable(io.en, coeff(0))

  } else {

    val tableI = SinCosSim.sincosTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(cbit(i) <= maxCbit(i))

      if(0 < diffWidth) {
        val ci  = coeff(i)
        val msb = ci(cbit(i)-1)
        coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffs.cs(i) := coeff(i)
      }
    }
    io.cs := enable(io.en, coeffs)
  }
}

object SinCosTableCoeff {
  def getCBits(
    spec:     RealSpec,
    polySpec: PolynomialSpec
  ): Seq[Int] = {

    val order     = polySpec.order
    val adrW      = polySpec.adrW
    val extraBits = polySpec.extraBits
    val fracW     = polySpec.fracW

    if(order == 0) {
      return Seq(fracW)
    } else {
      return SinCosSim.
        sincosTableGeneration( order, adrW, spec.manW, fracW, None, None ).
          getCBitWidth(/*signmode=*/0)
    }
  }
  def getCalcW(
    spec:     RealSpec,
    polySpec: PolynomialSpec
  ): Seq[Int] = {

    val order     = polySpec.order
    val adrW      = polySpec.adrW
    val extraBits = polySpec.extraBits
    val fracW     = polySpec.fracW

    if(order == 0) {
      return Seq(fracW)
    } else {
      return SinCosSim.
        sincosTableGeneration( order, adrW, spec.manW, fracW, None, None ).
          calcWidth
    }
  }
}

// -------------------------------------------------------------------------
//                        _        _     _                    _   _
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_|
//                                               |_|
// -------------------------------------------------------------------------

// sincos does not need any non-table path.
// sincos table approximates sin(x)/(x/pi) for x in [0, pi/2).
// postprocess calculates sin(x)/(x/pi) * x/pi to obtain sin(x).
// this x/pi term gains precision where x is around 0.

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class SinCosPostProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  val nStage = stage.total
  def getStage = nStage

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order
  val extraBits = fracW - manW

  val io = IO(new Bundle {
    val en   = Input(UInt(1.W))
    val pre  = Flipped(new SinCosPreProcessOutput(spec))
    val zman0  = Input(UInt(manW.W))
    val zexInc = Input(UInt(1.W))
    val z    = Output(UInt(spec.W.W))
  })

  val zman0  = io.zman0
  val zex0   = io.pre.yex + 1.U + io.zexInc

  val znan = io.pre.znan
  val zone = io.pre.zone
  val zzero = io.pre.yex === 0.U

  val zsgn = io.pre.zsgn
  val zex  = Mux(zzero, 0.U, Mux(znan, maskI(exW).U(exW.W), Mux(zone, exBias.U(exW.W), zex0)))
  val zman = Mux(zzero, 0.U, Mux(znan || zone, Cat(znan, 0.U((manW-1).W)), zman0))

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := ShiftRegister(enable(io.en, z), nStage)
}
