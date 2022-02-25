//% @file sincos.scala
//
// x -> sin(x) or cos(x)
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

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

import rial.mathfunc._

//
// Evaluate sin(x) or cos(x) for the full range of x. Both function shares
// most of the part by converting cos(x) to sin(pi/2-x).
// So only SinCosPreProcess has conditional branch that depends on `isSin:Bool`.
//

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

// this contains x/pi rounded into [0, 2) and zsgn.
// zsgn can be calculated from x/pi.
class SinCosPreProcessOutput(val spec: RealSpec) extends Bundle {
  val znan          = Output(Bool())
  val zsgn          = Output(UInt(1.W))
  val xConvertedEx  = Output(UInt(spec.exW.W))
  val xConvertedMan = Output(UInt(spec.manW.W))
}

class SinCosPreProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias
  val exAdrW = MathFuncSinSim.calcExAdrW(spec)

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val io = IO(new Bundle {
    val en    = Input (UInt(1.W))
    val isSin = Input(Bool())
    val x     = Input (UInt(spec.W.W))
    val adr   = Output(UInt((exAdrW+adrW).W))
    val dx    = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
    val xConverted = new SinCosPreProcessOutput(spec)
  })

  val (xsgn, xex, xman) = FloatChiselUtil.decompose(spec, io.x & Fill(spec.W, io.en))

  // --------------------------------------------------------------------------
  // calc x/pi

  val oneOverPiPad = 23
  val oneOverPi  = BigInt(math.round(1.0 / math.Pi * pow(2.0, manW+2+oneOverPiPad)))
  val oneOverPiW = 1 + manW + oneOverPiPad

  val xOverPiProd          = Cat(1.U(1.W), xman) * oneOverPi.U(oneOverPiW.W)
  val xOverPiProdMoreThan2 = xOverPiProd((1+manW) + oneOverPiW - 1)
  val xOverPiEx            = xex -& 2.U + xOverPiProdMoreThan2
  val xOverPiIsZero        = xOverPiEx(exW)
  // always x.ex > xOverPiEx because xOverPiProdMoreThan2 <= 1. Overflow does
  // not mean inf.

  val xOverPi              = Mux(xOverPiProdMoreThan2, xOverPiProd,
                                 Cat(xOverPiProd.tail(1), 0.U(1.W)))
  val xOverPiFracW         = (1+manW) + oneOverPiW - 1
  assert(xOverPi.getWidth == xOverPiFracW+1) // fraction bits + 1 integer bit

//   printf("xOverPi      = %d\n", xOverPi)
//   printf("xOverPiFracW = %d\n", xOverPiFracW.U)

  // --------------------------------------------------------------------------
    // convert full range x/pi into (0, 2)

  val xOverPiExMoreThan1 = (xOverPiEx >= exBias.U)
  val xOverPiExNobiasPos = (xOverPiEx - exBias.U)(log2Up(xOverPiFracW)-1, 0)
  val xOverPiExNobiasNeg = (exBias.U - xOverPiEx)(log2Up(xOverPiFracW)-1, 0)

  val xOverPiAlignPos = xOverPi << xOverPiExNobiasPos
  val xOverPiAlignNeg = xOverPi >> xOverPiExNobiasNeg
  val xOverPiAligned = Mux(xOverPiExMoreThan1, xOverPiAlignPos(xOverPiFracW, 0),
                                               xOverPiAlignNeg(xOverPiFracW, 0))
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
  io.xConverted.zsgn := ShiftRegister(zsgn, nStage)

  // --------------------------------------------------------------------------
  // special value check: if x is inf or nan, then z is nan.

  val znan = ~(xex.orR) // if xex == 1,
  io.xConverted.znan := ShiftRegister(znan, nStage)

  // --------------------------------------------------------------------------
  // calculate x/pi in sin/cos [0, 2) into sin [0, 1/2)

  val ymanPos0 =  xOverPiAligned(1+xOverPiFracW-2, 0)
  val ymanNeg0 = ((BigInt(1)<<(1+xOverPiFracW-1)).U - xOverPiAligned(1+xOverPiFracW-2, 0))(1+xOverPiFracW-2, 0)
  val ymanPos  = Cat(ymanPos0.head(1) & io.isSin, ymanPos0.tail(1))
  val ymanNeg  = Cat(ymanNeg0.head(1) & io.isSin, ymanNeg0.tail(1))
  val yman0    = Mux(xOverPiAligned2MSBs(0) ^ io.isSin, ymanPos, ymanNeg)
  val ymanIsNonZero = yman0.orR

//   printf("cir: ymanPos0 = %d\n", ymanPos0)
//   printf("cir: ymanNeg0 = %d\n", ymanNeg0)
//   printf("cir: ymanPos  = %d\n", ymanPos )
//   printf("cir: ymanNeg  = %d\n", ymanNeg )
//   printf("cir: yman0    = %d\n", yman0   )

  // TODO here we can reduce the area by checking yman0Shift and roundbit

  // yman0 removes its MSB, so here the width becomes xOverPiFracW,
  // not xOverPiFracW+1. But yman0 should be aligned to xOverPiFracW+1.
  // So here we need to add 1.
  val yman0Shift0   = PriorityEncoder(Reverse(yman0)) + 1.U
  val yman0ShiftW   = log2Up(yman0.getWidth)
  val yman0Shift    = yman0Shift0(yman0ShiftW-1, 0)
  val yman0Shifted  = (yman0 << yman0Shift)(1+xOverPiFracW-1, 0)
  val yman0RoundBit = xOverPiFracW - manW
  val yman0Rounded  = yman0Shifted(1+xOverPiFracW-1, yman0RoundBit) +&
                      yman0Shifted(yman0RoundBit-1)
  val yman0MoreThan2 = yman0Rounded(manW+1)

//   printf("cir:yman0Shift0  = %d\n", yman0Shift0)
//   printf("cir:yman0.W      = %d\n", yman0.getWidth.U)
//   printf("cir:yman0ShiftW  = %d\n", yman0ShiftW.U)
//   printf("cir:yman0Shift   = %d\n", yman0Shift)
//   printf("cir:yman0Shifted = %d\n", yman0Shifted)
//   printf("cir:yman0Rounded = %b\n", yman0Rounded)

  val ymanNonZero = yman0Rounded(manW-1, 0)
  val yexNonZero = exBias.U(exW.W) - yman0Shift + yman0MoreThan2

  val yex  = Fill(exW,  ymanIsNonZero) & yexNonZero  // if it's zero, return zero.
  val yman = Fill(manW, ymanIsNonZero) & ymanNonZero

  val exOfs = (exBias-2).U - yex
  val exAdr = exOfs(exAdrW-1, 0)

  val adr0 = Cat(exAdr, yman(manW-1, manW-adrW))
  val adr  = adr0 & Fill(adr0.getWidth, io.en)
  io.adr   := ShiftRegister(adr,  nStage)

  if(order != 0) {
    val dx0  = Cat(~yman(manW-adrW-1), yman(manW-adrW-2,0))
    val dx   = dx0 & Fill(dx0.getWidth, io.en)
    io.dx.get := ShiftRegister(dx, nStage)
  }

  val xConvertedEx  = yex  & Fill(yex .getWidth, io.en)
  val xConvertedMan = yman & Fill(yman.getWidth, io.en)

//   printf("cir:yex  = %d\n", yex)
//   printf("cir:yman = %b\n", yman)

  io.xConverted.xConvertedEx  := ShiftRegister(xConvertedEx,  nStage)
  io.xConverted.xConvertedMan := ShiftRegister(xConvertedMan, nStage)
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
  val maxAdrW  : Int,      // max address width among all math funcs
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
  val stage    : PipelineStageConfig
) extends Module {

  val manW   = spec.manW
  val nStage = stage.total

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val taylorThreshold = MathFuncSinSim.calcTaylorThreshold(manW)

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input  (UInt(maxAdrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  val exAdrW = MathFuncSinSim.calcExAdrW(spec)
  val exAdr  = io.adr(adrW+exAdrW-1, adrW)
  val manAdr = io.adr(adrW-1, 0)
  assert(adrW+exAdrW <= maxAdrW)

  if(order == 0) {
    val tbl = VecInit( (-2 to taylorThreshold by -1).map( exponent => {
      VecInit((0L to (1L << adrW)).map(
        n => {
          val x = n.toDouble / (1L<<adrW)
          val y = round(scalb(math.sin(Pi * scalb(1.0 + x, exponent)), -exponent-3) * (1L<<fracW))
          assert(y < (1L<<fracW))
          y.U(fracW.W)
        }))
      })
    )

    assert(maxCbit(0) == fracW)

    val c0 = tbl(exAdr)(manAdr)
    val c  = c0 & Fill(c0.getWidth, io.en)
    io.cs.cs(0) := ShiftRegister(c, nStage) // width should be manW + extraBits

  } else {

    val cbit = MathFuncSinSim.sinTableGeneration( order, adrW, manW, fracW )
      .map( t => {t.getCBitWidth(/*sign mode = */0)} )
      .reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

    val tableIs = VecInit(
      MathFuncSinSim.sinTableGeneration( order, adrW, manW, fracW ).map(t => {
        t.getVectorWithWidth(cbit, /*sign mode = */ 0)
      })
    )
    val tableI = tableIs(exAdr)

    val coeff = getSlices(tableI(manAdr), cbit)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(cbit(i) <= maxCbit(i))
      val ci  = coeff(i)
      val msb = ci(cbit(i)-1)
      if(0 < diffWidth) {
        coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffs.cs(i) := ci
      }
    }
    val cs = coeffs.asUInt & Fill(coeffs.asUInt.getWidth, io.en)
    io.cs := ShiftRegister(cs.asTypeOf(new TableCoeffInput(maxCbit)), nStage)
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

class SinCosNonTableOutput(val spec: RealSpec) extends Bundle {
  val zIsNonTable = Output(Bool())
  val zsgn        = Output(UInt(1.W))
  val zex         = Output(UInt(spec.exW.W))
  val zman        = Output(UInt(spec.manW.W))
}

class SinCosOtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW      = polySpec.adrW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val io = IO(new Bundle {
    val x          = Flipped(new DecomposedRealOutput(spec))
    val xConverted = Flipped(new SinCosPreProcessOutput(spec))
    val zother     = new SinCosNonTableOutput(spec)
  })

  val zSgn = io.xConverted.zsgn
  val yex  = io.xConverted.xConvertedEx
  val yman = io.xConverted.xConvertedMan

  val znan  = io.xConverted.znan
  val zzero = yex === 0.U
  val zone  = yex === (exBias-1).U && yman === 0.U

  val coefPad    = 2 // for precision after rounding
  val fracW      = manW+coefPad

  // -------------------------------------------------------------------------
  // Here, we also check if the linear approx is available.
  // Since taylorExpansion uses 5th order term, the exponent might overflows
  // while calculation. But with FP32, fortunately, linear approx becomes
  // usable (x^2 < -(fracW+2) <=> x.ex < -14) before exponent overflow
  // (x^5 < -exBias <=> x < -127/5 ~ -25). This holds for most of the "normal"
  // floating point spec. Like: in case of FP64, (x^2 < -(manW+1) <=> x.ex < -27)
  // and (x^5 < -exBias <=> x < -1023 / 5 ~ -200).
  //     This quick workaround will be broken if we have an abnormal FP spec
  // something like exponent has extremely narrower bit width than that of
  // mantissa.
  //     Note that, pi * y does not underflow because pi > 2. Its exponent is
  // positive (without exBias), so it never overflows but underflows. And here,
  // y is rounded into [0, 1/2) range. In this range, pi * y never overflows.
  // So here we don't need to be afraid of over/under flow of exponent.
  //     The threshold exponent is the maximum y such that y^2 does not affect
  // to the rounding bit of 1 - pi^2/6 y^2. That means that
  //     pi^2/6 y^2 < 2^(-fracW-1)
  //     pi^2/6 y^2 < 2 y^2 < 2^(-fracW-1)
  //                    y^2 < 2^(-fracW-2)
  //                    y   < 2^floor((-fracW-2)/2)
  //                    yex < floor((-fracW-2)/2) - 1 (because: 1 <= man < 2)
  //
  //
  val linearThreshold = math.floor(-(fracW+2)/2).toInt - 1
  val isLinear        = yex < (linearThreshold + exBias).U(exW.W)
//   printf("cir: isLinear        = %b\n", isLinear)
//   printf("cir: linearThreshold = %d\n", (linearThreshold + exBias).U(exW.W))

  // --------------------------------------------------------------------------
  // taylor expansion

  val taylorThreshold = MathFuncSinSim.calcTaylorThreshold(manW)
  val isTaylor        = yex < (taylorThreshold + exBias).U(exW.W)

  // XXX if y does not use Taylor path, should we need to zero clear yex and
  // yman to avoid needless voltage change for the sake of energy efficiency?

  val coef1Ex    = 1 + exBias
  val coef1ManW1 = math.round(Pi * (1<<(fracW-(coef1Ex-exBias)))).toLong
  val coef3Ex    = 0 + exBias
  val coef3ManW1 = math.round(Pi * Pi / 6.0 * (1<<(fracW-(coef3Ex-exBias)))).toLong
  val coef5Ex    = -1 + exBias
  val coef5ManW1 = math.round(pow(Pi, 4) / 120.0 * (1<<(fracW-(coef5Ex-exBias)))).toLong

  // helper function because we will multiply a lot in this path
  val multiply = (xFracW: Int, xmanW1: UInt, yFracW: Int, ymanW1: UInt) => {
    assert(xmanW1.getWidth == xFracW + 1)
    assert(ymanW1.getWidth == yFracW + 1)
    assert(xmanW1(xFracW) === 1.U)
    assert(ymanW1(yFracW) === 1.U)

    val zProd      = xmanW1 * ymanW1
    val zMoreThan2 = zProd((xFracW+1) + (yFracW+1) - 1)

    // .---- zProd(1+xFracW+1+yFracW-1)
    // |.--- zProd(1+xFracW+1+yFracW-2)
    // || .- zProd(1+xFracW+1+yFracW-3)
    // vv v
    // xx.xxxxxxx
    //
    val zShifted   = Mux(zMoreThan2,
      zProd((xFracW+1)+(yFracW+1)-2, (xFracW+1)+(yFracW+1)-2-fracW),
      zProd((xFracW+1)+(yFracW+1)-3, (xFracW+1)+(yFracW+1)-3-fracW))
    assert(zShifted.getWidth == fracW+1)

    val zRounded   = zShifted(fracW, 1) +& zShifted(0)
    assert(zRounded.getWidth == fracW+1)

    val zLessThan2AfterRound = ~zRounded(fracW)
    val zmanW1     = Cat(1.U(1.W), Fill(fracW, zLessThan2AfterRound) & zRounded(fracW-1, 0))
    val zexInc     = zMoreThan2 + (~zLessThan2AfterRound) // does not overflow.
    (zexInc, zmanW1)
  }

  // .........................................................................
  // 1st step

  // y^2
  val ymanW1 = Cat(1.U(1.W), yman)
  val (ySqExInc, ySqManW1) = multiply(manW, ymanW1, manW, ymanW1)
  val ySqEx = yex +& yex - exBias.U + ySqExInc
  assert(ySqManW1(fracW) === 1.U)
//   printf("cir:ySqEx    = %d\n", ySqEx)
//   printf("cir:ySqManW1 = %b\n", ySqManW1)

  // pi*y
  val (piyExInc, piyManW1) = multiply(manW, ymanW1, fracW, coef1ManW1.U((1+fracW).W))
  val piyEx = yex +& (coef1Ex - exBias).U + piyExInc
  assert(coef1Ex > exBias)
  assert(piyManW1(fracW) === 1.U)
//   printf("cir:yex      = %d\n", yex)
//   printf("cir:coef1ex  = %d\n", coef1Ex.U)
//   printf("cir:piyExInc = %d\n", piyExInc)
//   printf("cir:piyEx    = %d\n", piyEx)
//   printf("cir:piyManW1 = %b\n", piyManW1)

  // .........................................................................
  // 2nd step

  // y^4
  val (yQdExInc, yQdManW1) = multiply(fracW, ySqManW1, fracW, ySqManW1)
  val yQdEx = ySqEx + ySqEx - exBias.U + yQdExInc
  assert(yQdManW1(fracW) === 1.U)
//   printf("cir:yQdEx    = %d\n", yQdEx)
//   printf("cir:yQdManW1 = %b\n", yQdManW1)

  // pi^2y^2/6
  val (c3ExInc, c3ManW1) = multiply(fracW, ySqManW1, fracW, coef3ManW1.U((1+fracW).W))
  val c3Ex = ySqEx + (coef3Ex - exBias).U + c3ExInc
  assert(c3ManW1(fracW) === 1.U)
//   printf("cir:c3Ex    = %d\n", c3Ex)
//   printf("cir:c3ManW1 = %b\n", c3ManW1)

  // .........................................................................
  // 3rd step

  // pi^4y^4/120
  val (c5ExInc, c5ManW1) = multiply(fracW, yQdManW1, fracW, coef5ManW1.U((1+fracW).W))
  val c5Ex = yQdEx - (exBias - coef5Ex).U + c5ExInc
  assert(coef5Ex < exBias)
  assert(c5ManW1(fracW) === 1.U)
//   printf("cir:c5Ex    = %d\n", c5Ex)
//   printf("cir:c5ManW1 = %b\n", c5ManW1)

  // 1 - pi^2y^2/6
  assert(c3Ex+1.U < exBias.U || isLinear) // detect ex overflow
  assert(exBias - (1+fracW) > 0)
  val c3Shift0   = ((exBias-1).U - c3Ex)
  val c3ShiftOut = c3Shift0(exW-1, log2Up(1+fracW)).orR
  val c3Shift    = c3Shift0(log2Up(1+fracW)-1, 0)
  val c3Aligned0 = Fill(fracW+1, ~c3ShiftOut) & (c3ManW1 >> c3Shift) // if shiftout, fill 0.
  val c3Aligned  = c3Aligned0(c3Aligned0.getWidth-1, 1) +& c3Aligned0(0)
  assert(c3Aligned(fracW) === 0.U)

  val oneMinusC3 = (1L << fracW).U((fracW+1).W) - c3Aligned
  assert(oneMinusC3(fracW) === 0.U || oneMinusC3(fracW-1,0).orR === 0.U)

//   printf("cir:1-c3 = %b\n", oneMinusC3)

  // .........................................................................
  // 4th step

  // 1 - pi^2y^2/6 + pi^4y^4/120
  // ~ 1 - 1.645y^2 + 0.8117y^4 <= 1
  assert(c5Ex+1.U < exBias.U || isLinear)
  val c5Shift0   = ((exBias-1).U - c5Ex)
  val c5ShiftOut = c5Shift0(exW-1, log2Up(1+fracW)).orR
  val c5Shift    = c5Shift0(log2Up(1+fracW)-1, 0)
  val c5Aligned0 = Fill(fracW+1, ~c5ShiftOut) & (c5ManW1 >> c5Shift)
  val c5Aligned  = c5Aligned0(c5Aligned0.getWidth-1, 1) +& c5Aligned0(0)
  assert(c3Aligned >= c5Aligned)

  val oneMinusC3PlusC5 = oneMinusC3 + c5Aligned
  assert(oneMinusC3PlusC5 <= (1<<fracW).U)

//   printf("cir:1-c3+c5 = %b\n", oneMinusC3PlusC5)

  // in case of x < 2^-manW, the result just 1
  val oneMinusC3PlusC5MoreThan1 = oneMinusC3PlusC5(fracW)
  val oneMinusC3PlusC5ManW1 = Mux(oneMinusC3PlusC5MoreThan1, oneMinusC3PlusC5,
    Cat(oneMinusC3PlusC5(fracW-1, 0), 0.U(1.W))) // normalize
  val oneMinusC3PlusC5Ex = (exBias - 1).U + oneMinusC3PlusC5MoreThan1
  assert(oneMinusC3PlusC5ManW1.getWidth == fracW+1)

//   printf("cir:1-c3+c5 Ex    = %d\n", oneMinusC3PlusC5Ex)
//   printf("cir:1-c3+c5 ManW1 = %b\n", oneMinusC3PlusC5ManW1)

  // .........................................................................
  // 5th step

  // piy * (1 - pi^2y^2/6 + pi^4y^4/120)
  val (taylorExInc, taylorManW1) = multiply(fracW, piyManW1, fracW, oneMinusC3PlusC5ManW1)
  val taylorManW1Rounded = taylorManW1(fracW, coefPad) +& taylorManW1(coefPad-1)
  val taylorManW1MoreThan2AfterRound = taylorManW1Rounded(2+manW-1)

  val zExTaylor  = piyEx + oneMinusC3PlusC5Ex - exBias.U + taylorExInc + taylorManW1MoreThan2AfterRound
  val zManTaylor = taylorManW1Rounded(manW-1, 0)
//   printf("cir:taylorEx  = %d\n", zExTaylor)
//   printf("cir:taylorMan = %b\n", zManTaylor)

  // --------------------------------------------------------------------------
  // linear approx

  val zManLinear0 = piyManW1(fracW-1, coefPad) +& piyManW1(coefPad-1)
  assert(zManLinear0.getWidth == 1+manW)
  val zManLinearMoreThan2AfterRound = zManLinear0(manW)
  val zManLinear = zManLinear0(manW-1, 0)
  val zExLinear  = piyEx + zManLinearMoreThan2AfterRound
//   printf("cir:linearEx  = %d\n", zExLinear)
//   printf("cir:linearMan = %b\n", zManLinear)

  // --------------------------------------------------------------------------
  // merge the results (special value or taylor expansion)

  val zIsNonTable = znan || zzero || zone || isTaylor
  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)

  val zEx = Mux(isLinear, zExLinear,
            Mux(isTaylor, zExTaylor,
            Mux(znan,     Fill(exW, 1.U(1.W)),
            Mux(zone,     exBias.U,
            Mux(zzero,    0.U,
                          yex)))))
  // polynomial uses yex as the starting point and add some correction later.
  // If `y` does not match any of special cases, return yex.

  val zMan = Mux(isLinear, zManLinear,
             Mux(isTaylor, zManTaylor,
                           Cat(znan, 0.U((manW-1).W))))

  io.zother.zsgn := ShiftRegister(zSgn, nStage)
  io.zother.zex  := ShiftRegister(zEx,  nStage)
  io.zother.zman := ShiftRegister(zMan, nStage)
}

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
  def getStage() = nStage

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val io = IO(new Bundle {
    val en = Input(UInt(1.W))
    val zother = Flipped(new SinCosNonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val z      = Output(UInt(spec.W.W))
  })

  // --------------------------------------------------------------------------
  // postprocess the polynomial result.

  val resLessThanHalf = io.zres(fracW-1) === 0.U
//   printf("cir:zres = %b\n", io.zres)

  val zEx0Table  = Mux(resLessThanHalf, io.zother.zex + 1.U, io.zother.zex + 2.U)
  val zExTable   = zEx0Table(exW-1, 0)

  val zMan0Table = Mux(resLessThanHalf, Cat(io.zres, 0.U(2.W))(fracW-1, 0),
                                        Cat(io.zres, 0.U(1.W))(fracW-1, 0))
  val zManTable  = zMan0Table(fracW-1, fracW-manW) + zMan0Table(fracW-manW-1)
//   printf("cir:zmanTable = %b\n", zManTable)

  // --------------------------------------------------------------------------
  // merge with the result from non-table path

  val zSgn = io.zother.zsgn
  val zEx  = Mux(io.zother.zIsNonTable, io.zother.zex,  zExTable)
  val zMan = Mux(io.zother.zIsNonTable, io.zother.zman, zManTable)

  val z0  = Cat(zSgn, zEx, zMan)
  val z   = z0 & Fill(z0.getWidth, io.en)
  io.z   := ShiftRegister(z, nStage)
}
