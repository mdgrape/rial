//% @file atan2.scala
//
// ATan2 function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import java.lang.Math.scalb
import scala.language.reflectiveCalls
import scala.math._

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
import rial.arith._

import rial.math._

// ATan2 Stage1 calculates min(x,y)/max(x,y).
//       Stage2 calculates atan(min(x,y)/max(x,y)) +/- constant.
// Some flags are needed to be saved.

object ATan2Status {
  val W = 2
  val xIsPosIsLarger  = 0.U(W.W) // ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|
  val xIsNegIsLarger  = 1.U(W.W) // ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y|
  val xIsPosIsSmaller = 2.U(W.W) // ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y|
  val xIsNegIsSmaller = 3.U(W.W) // ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y|
}
object ATan2SpecialValue {
  val W = 3
  val zNormal     = 0.U(W.W)
  val zNaN        = 1.U(W.W) // == nan
  val zZero       = 2.U(W.W) // == zero
  val zPi         = 3.U(W.W) // == pi
  val zHalfPi     = 4.U(W.W) // == pi/2
  val zQuarterPi  = 5.U(W.W) // == pi/4
  val z3QuarterPi = 6.U(W.W) // == 3pi/4
}
class ATan2Flags extends Bundle {
  val status  = UInt(ATan2Status.W.W)
  val special = UInt(ATan2SpecialValue.W.W)
  val ysgn    = UInt(1.W)
}

// =========================================================================
//      _                     _
//  ___| |_ __ _  __ _  ___  / |
// / __| __/ _` |/ _` |/ _ \ | |
// \__ \ || (_| | (_| |  __/ | |
// |___/\__\__,_|\__, |\___| |_|
//               |___/
// =========================================================================

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------
//
// Stage1 preprocess only checks the special cases.
// Stage1 calculates min(x,y)/max(x,y), so ReciprocalPreProcess is re-used.
//
class ATan2Stage1PreProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val dxW    = polySpec.dxW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val en      = Input (UInt(1.W))
    val x       = Flipped(new DecomposedRealOutput(spec))
    val y       = Flipped(new DecomposedRealOutput(spec))
    val special = Output(UInt(ATan2SpecialValue.W.W))
  })

  val samexy    = io.x.ex === io.y.ex && io.x.man === io.y.man

  val xpos      = io.x.sgn === 0.U
  val xneg      = io.x.sgn === 1.U
  val znan      =  (io.x.nan ||  io.y.nan) || ( io.x.zero &&  io.y.zero)
  val zzero     = ((io.x.inf && !io.y.inf) || (!io.x.zero &&  io.y.zero)) && xpos
  val zpi       = ((io.x.inf && !io.y.inf) || (!io.x.zero &&  io.y.zero)) && xneg
  val zhalfpi   = (!io.x.inf &&  io.y.inf) || ( io.x.zero && !io.y.zero)
  val z1piover4 = ((io.x.inf &&  io.y.inf) || (samexy && !io.x.nan && !io.y.nan)) && xpos
  val z3piover4 = ((io.x.inf &&  io.y.inf) || (samexy && !io.x.nan && !io.y.nan)) && xneg

  val special0 = MuxCase(ATan2SpecialValue.zNormal, Seq(
    znan      -> ATan2SpecialValue.zNaN,
    zzero     -> ATan2SpecialValue.zZero,
    zpi       -> ATan2SpecialValue.zPi,
    zhalfpi   -> ATan2SpecialValue.zHalfPi,
    z1piover4 -> ATan2SpecialValue.zQuarterPi,
    z3piover4 -> ATan2SpecialValue.z3QuarterPi
  ))
  val special = enable(io.en, special0)

  io.special := ShiftRegister(special, nStage)
}

// -------------------------------------------------------------------------
//                        _        _     _                    _   _
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_|
//                                               |_|
// -------------------------------------------------------------------------

class ATan2Stage1NonTableOutput(val spec: RealSpec) extends Bundle {
  val zex       = Output(UInt(spec.exW.W)) // sign of min(x,y)/max(x,y)
  val maxXYMan0 = Output(Bool())           // max(x,y).man === 0.U
  val xySameMan = Output(Bool())
}

class ATan2Stage1OtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val yIsLarger = Input(Bool())
    val x       = Flipped(new DecomposedRealOutput(spec))
    val y       = Flipped(new DecomposedRealOutput(spec))

    val zother  = new ATan2Stage1NonTableOutput(spec)
  })

  val maxXYMan0 = Mux(io.yIsLarger, (!io.y.man.orR.asBool), (!io.x.man.orR.asBool))
//   printf("x.man = %b\n", io.x.man)
//   printf("y.man = %b\n", io.y.man)
//   printf("x < y = %b\n", io.yIsLarger)

  val xySameMan = io.x.man === io.y.man

  io.zother.maxXYMan0 := ShiftRegister(maxXYMan0, nStage)
  io.zother.xySameMan := ShiftRegister(xySameMan, nStage)

  // --------------------------------------------------------------------------
  // Here we don't need to check if 1/max(x,y) is a special value because
  // Stage1PreProcess.io.special covers all the cases.

  // 1/x = 2^(-e-1) * 2/1.m
  // ex = -(x.ex - exBias) - 1 + exBias
  //    = -x.ex + exBias - 1 + exBias
  //    = exBias * 2 - 1 - x.ex
  //
  // y/x = 2^(y.e) * y.m * 2^(-x.e-1) * 2/x.m
  //     = 2^(y.e - x.e - 1) * y.m * (2/x.m)
  //
  // y/x.ex = (y.e - exBias - x.e + exBias - 1 + exBias
  //        = (y.e - x.e - 1 + exBias)

  val xexBiased = Mux(io.yIsLarger, io.y.ex, io.x.ex) // if y>x, swap x and y
  val yexBiased = Mux(io.yIsLarger, io.x.ex, io.y.ex)
  val zex0 = Wire(UInt(exW.W))

  // Since y < x, ex is always smaller than exBias. In normal cases, exBias is
  // around a half of 2^exW-1, so we have almost a half of the space of the
  // output port, UInt(exW.W). We re-interpret it as a signed integer to keep
  // information. If we round the negative value to zero, the postprocess
  // might consider the result is a small but non-zero value, though actually
  // that is less than the minimum.
  zex0 := yexBiased - xexBiased + (exBias-1).U(exW.W)

  // exponent of min(x,y)/max(x,y). we will later correct +/- 1 by checking
  // the mantissa of min(x,y) and max(x,y)
  io.zother.zex := ShiftRegister(zex0, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

// Stage1: takes 1/max(x, y) and min(x, y), returns min(x, y) / max(x, y)

class ATan2Stage1PostProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val nStage = stage.total
  def getStage = nStage

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order
  val extraBits = polySpec.extraBits

  val io = IO(new Bundle {
    val en = Input(UInt(1.W))
    val zother = Flipped(new ATan2Stage1NonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val minxy  = Flipped(new DecomposedRealOutput(spec))
    val z      = Output(UInt(spec.W.W))
  })

  val zsgn      = 0.U(1.W)
  val zex0      = io.zother.zex
  val maxXYMan0 = io.zother.maxXYMan0
  val xySameMan = io.zother.xySameMan

//   printf("cir: zex0 = %b\n", zex0)
//   printf("cir: zres = %b\n", io.zres)

  val denomW1 = Cat(1.U(1.W), Mux(maxXYMan0, 0.U, io.zres))
  val numerW1 = Cat(1.U(1.W), io.minxy.man)

//   printf("cir: denomW1 = %b\n", denomW1)
//   printf("cir: numerW1 = %b\n", numerW1)

  val zProd     = denomW1 * numerW1
  val bp        = fracW + manW
  val roundBits = fracW + manW - manW
  val zProdMoreThan2 = zProd((fracW+1)+(manW+1)-1)
  val zProdSticky    = zProd(roundBits-2, 0).orR | (zProdMoreThan2 & zProd(roundBits-1))
  val zProdRound     = Mux(zProdMoreThan2, zProd(roundBits),       zProd(roundBits-1))
  val zProdShifted   = Mux(zProdMoreThan2, zProd(bp, roundBits+1), zProd(bp-1, roundBits))
  assert(zProdShifted.getWidth == manW)
  val zProdLSB       = zProdShifted(0)
  val zProdInc       = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
    zProdLSB, zProdRound, zProdSticky)
  val zProdRounded   = zProdShifted +& zProdInc
  assert(zProdRounded.getWidth == manW+1)
  val zProdMoreThan2AfterRound = zProdRounded(manW)

  val zex = Wire(UInt(exW.W))

  // the result of OtherPath might cause underflow.
  val zex0Inc = Mux(xySameMan || maxXYMan0, zex0 + 1.U,
    zex0 + zProdMoreThan2 + zProdMoreThan2AfterRound)

  val canUnderflow = (spec.exMin - spec.exMax + exBias < 0)
  if (canUnderflow) {
    zex := Mux(zex0Inc(exW-1), 0.U, zex0Inc)
  } else {
    zex := zex0Inc
  }

  val zman = Mux(~zex.orR || xySameMan, 0.U(manW.W),
             Mux(maxXYMan0, io.minxy.man, zProdRounded(manW-1, 0)))

  val z0 = Cat(zsgn, zex, zman)
  val z = enable(io.en, z0)

  io.z := ShiftRegister(z, nStage)
}

// =========================================================================
//      _                     ____
//  ___| |_ __ _  __ _  ___  |___ \
// / __| __/ _` |/ _` |/ _ \   __) |
// \__ \ || (_| | (_| |  __/  / __/
// |___/\__\__,_|\__, |\___| |_____|
//               |___/
// =========================================================================

// -------------------------------------------------------------------------
//                                                     ____
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___  |___ \
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|   __) |
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \  / __/
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/ |_____|
// |_|            |_|
// -------------------------------------------------------------------------
//
// Stage2 calculates atan(x), so we need to extract the address value

class ATan2Stage2PreProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val dxW    = polySpec.dxW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val en  = Input (UInt(1.W))
    val x   = Flipped(new DecomposedRealOutput(spec))
    val adr = Output(UInt(adrW.W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  // stage1 calculates min(x,y)/max(x,y), so we can assume that x <= 1.
  // also, if x == 1, stage1 sets special value flag. so we can ignore the case.

  val xex    = io.x.ex
  val xmanW1 = Cat(1.U(1.W), io.x.man)

  val shiftW    = log2Up(1+manW)
  val xShift0   = exBias.U(exW.W) - xex
  val xShiftOut = xShift0(exW-1, shiftW).orR
  val xShift    = Mux(xShiftOut, maskL(shiftW).U(shiftW.W), xShift0(shiftW-1, 0))
  val xFixed    = xmanW1 >> xShift

//   printf("cir: x       = %b|%d|%b\n", io.x.sgn, xex, io.x.man)
//   printf("cir: xShift  = %b\n", xShift)
//   printf("cir: xFixed  = %b\n", xFixed)

  // xFixed(manW) === 1 if and only if x == 1.0.
  assert(xFixed(manW) =/= 1.U || xFixed(manW-1, 0) === 0.U || io.en === 0.U,
         "xFixed = %b", xFixed)

  val adr  = enable(io.en, xFixed(manW-1, dxW))
  io.adr := ShiftRegister(adr, nStage)

//   printf("cir: adr = %b\n", adr)

  if(order != 0) {
    val dx   = enable(io.en, Cat(~xFixed(manW-adrW-1), xFixed(manW-adrW-2,0)))
//     printf("cir: dx  = %b\n", dx)
    io.dx.get := ShiftRegister(dx, nStage)
  }
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __   ____
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _| |___ \
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_    __) |
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|  / __/
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|   |_____|
// -------------------------------------------------------------------------
//
// Stage1 re-use the reciprocal table. we don't need to implement it for atan2.
//
class ATan2Stage2TableCoeff(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
) extends Module {

  val manW   = spec.manW
  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input  (UInt(adrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  val adr = io.adr

  if(order == 0) {

    val tableI = ATan2Stage2Sim.atanTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    // sign mode 0: always include sign
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode = */1)
    val coeff = getSlices(coeffTable(adr), coeffWidth)

    io.cs.cs(0) := enable(io.en, coeff(0))

  } else {

    val tableI = ATan2Stage2Sim.atanTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    // sign mode 0: always include sign
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode = */0)
    val coeff = getSlices(coeffTable(adr), cbit)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(0 <= diffWidth)
      if(diffWidth != 0) {
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

object ATan2Stage2TableCoeff {
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
      return ATan2Stage2Sim.atanTableGeneration( order, adrW, spec.manW, fracW ).
        getCBitWidth(/*sign mode = */0)
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
      return ATan2Stage2Sim.atanTableGeneration( order, adrW, spec.manW, fracW ).calcWidth
    }
  }
}


// -----------------------------------------------------------------------------
//                        _        _     _                    _   _       ____
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__   |___ \
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \    __) |
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |  / __/
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_| |_____|
//                                               |_|
// -----------------------------------------------------------------------------

class ATan2Stage2NonTableOutput(val spec: RealSpec) extends Bundle {
  val zsgn        = Output(UInt(1.W))
  val zex         = Output(UInt(spec.exW.W))
  val zman        = Output(UInt(spec.manW.W))
  val zIsLinear   = Output(Bool())
  val zIsNonTable = Output(Bool())
  val correctionNeeded = Output(Bool())
}

class ATan2Stage2OtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val flags   = Input(new ATan2Flags())
    val x       = Flipped(new DecomposedRealOutput(spec))
    val zother  = new ATan2Stage2NonTableOutput(spec)
  })

  val pi         = new RealGeneric(spec, Pi)
  val halfPi     = new RealGeneric(spec, Pi * 0.5)
  val quarterPi  = new RealGeneric(spec, Pi * 0.25)
  val quarter3Pi = new RealGeneric(spec, Pi * 0.75)

  val linearThreshold = (ATan2Stage2Sim.calcLinearThreshold(manW) + exBias)
  val isLinear = io.x.ex <= linearThreshold.U(exW.W)

  val xzero = !io.x.ex.orR

  val defaultEx  = io.x.ex // isLinear includes xzero.
  val defaultMan = Mux(xzero, 0.U(exW), io.x.man) // we need to re-set man if x is zero
  // non-linear mantissa is calculated by table.

  io.zother.zIsLinear        := isLinear
  io.zother.zIsNonTable      := isLinear || (io.flags.special =/= ATan2SpecialValue.zNormal)
  io.zother.correctionNeeded := isLinear || (io.flags.special === ATan2SpecialValue.zNormal)

  io.zother.zsgn := io.flags.ysgn
  io.zother.zex  := MuxCase(defaultEx, Seq(
    (io.flags.special === ATan2SpecialValue.zNaN)        -> Fill(exW, 1.U(1.W)),
    (io.flags.special === ATan2SpecialValue.zZero)       -> 0.U(exW.W),
    (io.flags.special === ATan2SpecialValue.zPi)         -> pi.ex.U(exW.W),
    (io.flags.special === ATan2SpecialValue.zHalfPi)     -> halfPi.ex.U(exW.W),
    (io.flags.special === ATan2SpecialValue.zQuarterPi)  -> quarterPi.ex.U(exW.W),
    (io.flags.special === ATan2SpecialValue.z3QuarterPi) -> quarter3Pi.ex.U(exW.W)
    ))
  io.zother.zman := MuxCase(defaultMan, Seq(
    (io.flags.special === ATan2SpecialValue.zNaN)        -> Fill(manW, 1.U(1.W)),
    (io.flags.special === ATan2SpecialValue.zZero)       -> 0.U(manW.W),
    (io.flags.special === ATan2SpecialValue.zPi)         -> pi        .man.toBigInt.U(manW.W),
    (io.flags.special === ATan2SpecialValue.zHalfPi)     -> halfPi    .man.toBigInt.U(manW.W),
    (io.flags.special === ATan2SpecialValue.zQuarterPi)  -> quarterPi .man.toBigInt.U(manW.W),
    (io.flags.special === ATan2SpecialValue.z3QuarterPi) -> quarter3Pi.man.toBigInt.U(manW.W)
    ))
}

// -------------------------------------------------------------------------
//                  _                                       ____
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___  |___ \
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|   __) |
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \  / __/
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/ |_____|
// |_|                 |_|
// -------------------------------------------------------------------------
//
// Stage2: takes status flags, returns corrected result
//

class ATan2Stage2PostProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val nStage = stage.total
  def getStage = nStage

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order
  val extraBits = polySpec.extraBits

  val io = IO(new Bundle {
    val en = Input(UInt(1.W))
    val x      = Flipped(new DecomposedRealOutput(spec))
    val zother = Flipped(new ATan2Stage2NonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val flags  = Input(new ATan2Flags()) // need status (|x|<|y|, xsgn)
    val z      = Output(UInt(spec.W.W))
  })

  val zSgn = io.zother.zsgn

  val zres = Cat(io.zres, 0.U(1.W))
//   printf("cir: io.zres = %b\n", io.zres)

  // atan2Stage2Sim: atanEx  = 115(-12)

  // sim: zres  = 11111111111111111111110111
  // cir: zres  = 11111111111111111111111000
  // cir: atanEx  = 115
  // sim: atanMan =  1100100100011111111100(3295228)
  // cir: atanMan = 01100100100011111111101

//   printf("cir: zres   = %b\n", zres  )
//   printf("cir: xmanW1 = %b\n", Cat(1.U(1.W), io.x.man))

  val atanProd = zres * Cat(1.U(1.W), io.x.man)
  val atanProdMoreThan2 = atanProd((fracW+1)+(manW+1)-1)
  val atanShifted = Mux(atanProdMoreThan2, atanProd((fracW+1)+(manW+1)-2, fracW+1),
                                           atanProd((fracW+1)+(manW+1)-3, fracW  ))
  val atanRoundInc = Mux(atanProdMoreThan2, atanProd(fracW), atanProd(fracW-1))
  val atanRound = atanShifted +& atanRoundInc
  val atanRoundMoreThan2 = atanRound(manW)

  val atanMan0 = atanRound(manW-1, 0)
  val atanEx0  = io.x.ex +& (atanProdMoreThan2 + atanRoundMoreThan2) - 1.U

  val atanEx  = Mux(io.zother.zIsNonTable, io.zother.zex,  atanEx0)
  val atanMan = Mux(io.zother.zIsNonTable, io.zother.zman, atanMan0)
//   printf("cir: isNonTab= %d\n", io.zother.zIsNonTable)
//   printf("cir: atanEx  = %d\n", atanEx)
//   printf("cir: atanMan = %b\n", atanMan)
// 
//   printf("atan = %d|%b\n", atanEx, atanMan)

  // ==========================================================================
  // select correction by:
  //   ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y| : xpos, xlarger
  //   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y| : xneg, xlarger
  //   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y| : xpos, xsmaller
  //   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y| : xneg, xsmaller
  //          ^^^^^^^^^^^^^^^^^^^
  //          this part is always positive
  //

  val pi         = new RealGeneric(spec, Pi)
  val halfPi     = new RealGeneric(spec, Pi * 0.5)

  //        1 +  manW   + 3
  //         .----------.
  // pi   = 1x.xxxxxx...x000
  // pi/2 = 01.xxxxxx...xx00
  // atan = 01.xxxxxx...xx00 (before shift)
  //           '---------'
  //        2 +  manW   +  2

//  val piManW1        = Cat(1.U(1.W), pi.man.toLong.U(manW.W),     0.U(3.W))
//  val halfPiManW1    = Cat(1.U(2.W), halfPi.man.toLong.U(manW.W), 0.U(2.W))
//  val atanManW1      = Cat(1.U(2.W), atanMan,                     0.U(2.W))

  val piManW1        = Real.pi(manW+2).toBigInt.U((manW+4).W)
  val halfPiManW1    = (Real.pi / Real.two)(manW+2).toBigInt.U((manW+4).W)
  val atanManW1      = Cat(1.U(2.W), atanMan, 0.U(2.W))

  assert(piManW1    .getWidth == 2+manW+2)
  assert(halfPiManW1.getWidth == 2+manW+2)
  assert(atanManW1  .getWidth == 2+manW+2)

  // --------------------------------------------------
  // align atan mantissa first

  val log2ManW3   = log2Up(1+manW+3)
  val atanShift0  = exBias.U(exW.W) - atanEx
  val atanShift   = Mux(atanShift0(exW-1, log2ManW3).orR,
                        (1+manW+3).U(log2ManW3.W), atanShift0(log2ManW3-1, 0))
  val atanAligned = Mux(atanShift > (manW+1+3).U(log2ManW3.W), 0.U((manW+1+3).W),
                        atanManW1 >> atanShift)

  // -------------------------------------------------------------------
  //   ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y| : xpos, xlarger

  val zExPL  = atanEx
  val zManPL = atanMan(manW-1, 0)

  // -------------------------------------------------------------------
  //   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y| : xneg, xlarger
  //
  // atan(x) is in [0, pi/4)
  // pi - atan(x) is in (3pi/4, pi] ~ (2.35.., 3.14..], ex is always 1

  val piMinusATanMan0 = piManW1 - atanAligned
  val zExNL  = (exBias+1).U(exW.W)
  val zManNL = piMinusATanMan0((1+manW+3)-2, 3) + piMinusATanMan0(2)

  // -------------------------------------------------------------------
  //   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y| : xpos, xsmaller
  //
  // atan(x) is in [0, pi/4)
  // pi/2 - atan(x) is in (pi/4, pi/2] ~ (0.78.., 1.57..], ex is -1 or 0

  val halfPiMinusATanMan0 = halfPiManW1 - atanAligned
  val halfPiMinusATanMan0MoreThan1 = halfPiMinusATanMan0((1+manW+3)-2)

  val zExPS  = (exBias-1).U(exW.W) + halfPiMinusATanMan0MoreThan1
  val zManPS = Mux(halfPiMinusATanMan0MoreThan1,
      halfPiMinusATanMan0((1+manW+3)-3, 2) + halfPiMinusATanMan0(1),
      halfPiMinusATanMan0((1+manW+3)-4, 1) + halfPiMinusATanMan0(0)
    )

  // -------------------------------------------------------------------
  //   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y| : xneg, xsmaller
  //
  // atan(x) is in [0, pi/4)
  // pi/2 + atan(x) is in (pi/2, 3pi/4] ~ (1.57.., 2.35..], ex is 0 or 1
  //
  // 01.00000000
  // ^  '----'^^
  // |  manW  2
  // morethan2

  val halfPiPlusATanMan0         = halfPiManW1 + atanAligned
  val halfPiPlusATanManRounded   = dropLSB(2, halfPiPlusATanMan0) + halfPiPlusATanMan0(1)
  val halfPiPlusATanManMoreThan2 = halfPiPlusATanManRounded(manW+1)

  val zExNS  = exBias.U(exW.W) + halfPiPlusATanManMoreThan2
  val zManNS = Mux(halfPiPlusATanManMoreThan2,
    halfPiPlusATanManRounded(manW,   1),
    halfPiPlusATanManRounded(manW-1, 0)
    )

//   printf("cir: pi/2+atan = %b\n", halfPiPlusATanMan0)
//   printf("cir: rounded   = %b\n", halfPiPlusATanManRounded)
//   printf("cir: moreThan2 = %b\n", halfPiPlusATanManMoreThan2)
//   printf("cir: zExNS     = %b\n", zExNS)
//   printf("cir: zManNS    = %b\n", zManNS)

  // -------------------------------------------------------------------
//   printf("PosLarger  = %d|%b\n", zExPL, zManPL)
//   printf("NegLarger  = %d|%b\n", zExNL, zManNL)
//   printf("PosSmaller = %d|%b\n", zExPS, zManPS)
//   printf("NegSmaller = %d|%b\n", zExNS, zManNS)
//   printf("status = %d\n", io.flags.status)

//   assert(!io.enable || (zres.getWidth == fracW).B)
//   assert(!io.enable || zres(fracW-1) === 1.U)
//   assert(!io.enable || zresRounded(manW) === 1.U)

  val zEx = MuxCase(0.U(exW.W), Seq(
      (io.flags.status === ATan2Status.xIsPosIsLarger ) -> zExPL,
      (io.flags.status === ATan2Status.xIsNegIsLarger ) -> zExNL,
      (io.flags.status === ATan2Status.xIsPosIsSmaller) -> zExPS,
      (io.flags.status === ATan2Status.xIsNegIsSmaller) -> zExNS
    ))
  val zMan = MuxCase(0.U(manW.W), Seq(
      (io.flags.status === ATan2Status.xIsPosIsLarger ) -> zManPL(manW-1, 0),
      (io.flags.status === ATan2Status.xIsNegIsLarger ) -> zManNL(manW-1, 0),
      (io.flags.status === ATan2Status.xIsPosIsSmaller) -> zManPS(manW-1, 0),
      (io.flags.status === ATan2Status.xIsNegIsSmaller) -> zManNS(manW-1, 0)
    ))

//   printf("cir: zex  = %b\n", zEx)
//   printf("cir: zman = %b\n", zMan)

  val zNormal = Cat(zSgn, zEx, zMan)

  // -------------------------------------------
  // select special cases

  val nan        = RealGeneric.nan(spec)
  val zero       = new RealGeneric(spec, 0)
  val quarterPi  = new RealGeneric(spec, Pi * 0.25)
  val quarter3Pi = new RealGeneric(spec, Pi * 0.75)

  val special = io.flags.special
  val z0 = MuxCase(zNormal, Seq(
      (special === ATan2SpecialValue.zNormal    ) -> zNormal,
      (special === ATan2SpecialValue.zNaN       ) -> Cat(0.U(1.W), nan   .ex.U(exW.W), nan       .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.zZero      ) -> Cat(0.U(1.W), zero  .ex.U(exW.W), zero      .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.zPi        ) -> Cat(zSgn, pi        .ex.U(exW.W), pi        .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.zHalfPi    ) -> Cat(zSgn, halfPi    .ex.U(exW.W), halfPi    .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.zQuarterPi ) -> Cat(zSgn, quarterPi .ex.U(exW.W), quarterPi .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.z3QuarterPi) -> Cat(zSgn, quarter3Pi.ex.U(exW.W), quarter3Pi.man.toBigInt.U(manW.W))
    ))
  val z = enable(io.en, z0)

  io.z := ShiftRegister(z, nStage)
}


// -------------------------------------------------------------------------
//                      _     _                _
//   ___ ___  _ __ ___ | |__ (_)_ __   ___  __| |
//  / __/ _ \| '_ ` _ \| '_ \| | '_ \ / _ \/ _` |
// | (_| (_) | | | | | | |_) | | | | |  __/ (_| |
//  \___\___/|_| |_| |_|_.__/|_|_| |_|\___|\__,_|
// -------------------------------------------------------------------------

class ATan2Generic(
  val spec: RealSpec,
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage1: MathFuncPipelineConfig, // x/y stage config
  val stage2: MathFuncPipelineConfig, // atan stagel config
  val stageGap: Boolean,              // register between stage1 result to stage2
  val dxW0 : Option[Int] = None,
  val enableRangeCheck: Boolean = true,
  val enablePolynomialRounding: Boolean = false
) extends Module {

  val pcGap1 = if(stage1.preCalcGap ) {1} else {0}
  val cpGap1 = if(stage1.calcPostGap) {1} else {0}

  val nPreStage1  = stage1.preStage.total
  val nCalcStage1 = stage1.calcStage.total
  val nPostStage1 = stage1.postStage.total

  val pcGap2 = if(stage2.preCalcGap ) {1} else {0}
  val cpGap2 = if(stage2.calcPostGap) {1} else {0}

  val nPreStage2  = stage2.preStage.total
  val nCalcStage2 = stage2.calcStage.total
  val nPostStage2 = stage2.postStage.total

  val sGap = if(stageGap) {1} else {0}

  val nStage   = stage1.total + sGap + stage2.total
  def getStage = nStage

  val polySpec = new PolynomialSpec(spec, nOrder, adrW, extraBits, dxW0,
    enableRangeCheck, enablePolynomialRounding)
  val order = polySpec.order

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val recCbits  = ReciprocalTableCoeff.getCBits(spec, polySpec)
  val recCalcW  = ReciprocalTableCoeff.getCalcW(spec, polySpec)
  val atanCbits = ATan2Stage2TableCoeff.getCBits(spec, polySpec)
  val atanCalcW = ATan2Stage2TableCoeff.getCalcW(spec, polySpec)

  def getCbit  = Seq(recCbits, atanCbits).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )
  def getCalcW = Seq(recCalcW, atanCalcW).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  val maxCbits = getCbit
  val maxCalcW = getCalcW

  val io = IO(new Bundle {
    val en = Input(Bool())
    val x = Input (UInt(spec.W.W))
    val y = Input (UInt(spec.W.W))
    val z = Output(UInt(spec.W.W))
  })

  // --------------------------------------------------------------------------

  val xdecomp = Module(new DecomposeReal(spec))
  val ydecomp = Module(new DecomposeReal(spec))
  xdecomp.io.real := io.x
  ydecomp.io.real := io.y

  val yIsLarger = io.x(spec.W-2, 0) < io.y(spec.W-2, 0) // cmp without sign bit

  val yIsLargerPCReg    = ShiftRegister(yIsLarger, nPreStage1)
  val yIsLargerPCGapReg = ShiftRegister(yIsLargerPCReg, pcGap1)
  val yIsLargerCPReg    = ShiftRegister(yIsLargerPCGapReg, nCalcStage1)
  val yIsLargerCPGapReg = ShiftRegister(yIsLargerCPReg, cpGap1)

  val enPC1Reg    = ShiftRegister(io.en,       nPreStage1)
  val enPC1GapReg = ShiftRegister(enPC1Reg,    pcGap1)
  val enCP1Reg    = ShiftRegister(enPC1GapReg, nCalcStage1)
  val enCP1GapReg = ShiftRegister(enCP1Reg,    cpGap1)

  val xdecPCReg    = ShiftRegister(xdecomp.io.decomp, nPreStage1)
  val xdecPCGapReg = ShiftRegister(xdecPCReg,         pcGap1)
  val xdecCPReg    = ShiftRegister(xdecPCGapReg,      nCalcStage1)
  val xdecCPGapReg = ShiftRegister(xdecCPReg,         cpGap1)

  val ydecPCReg    = ShiftRegister(ydecomp.io.decomp, nPreStage1)
  val ydecPCGapReg = ShiftRegister(ydecPCReg,         pcGap1)
  val ydecCPReg    = ShiftRegister(ydecPCGapReg,      nCalcStage1)
  val ydecCPGapReg = ShiftRegister(ydecCPReg,         cpGap1)

  // --------------------------------------------------------------------------

  val atan2Stage1Pre   = Module(new ATan2Stage1PreProcess (spec, polySpec, stage1.preStage))
  val atan2Stage1Other = Module(new ATan2Stage1OtherPath  (spec, polySpec, stage1.calcStage))
  val atan2Stage1Post  = Module(new ATan2Stage1PostProcess(spec, polySpec, stage1.postStage))
  val recPre           = Module(new ReciprocalPreProcess(spec, polySpec, stage1.preStage))
  val recTab           = Module(new ReciprocalTableCoeff(spec, polySpec, maxCbits))

  // atan2Stage1Pre checks if x and y are special values.
  // for calculation, reciprocal is re-used.
  atan2Stage1Pre.io.en := io.en
  atan2Stage1Pre.io.x  := xdecomp.io.decomp
  atan2Stage1Pre.io.y  := ydecomp.io.decomp

  recPre.io.en  := io.en
  recPre.io.x   := Mux(yIsLarger, ydecomp.io.decomp, xdecomp.io.decomp)

  // ------ Preprocess-Calculate ------
  atan2Stage1Other.io.x  := xdecPCGapReg
  atan2Stage1Other.io.y  := ydecPCGapReg
  atan2Stage1Other.io.yIsLarger := yIsLargerPCGapReg
  recTab.io.en  := enPC1GapReg
  recTab.io.adr := ShiftRegister(recPre.io.adr, pcGap1)

  val atan2FlagReg = Wire(new ATan2Flags())
  atan2FlagReg.status  := Cat(yIsLargerPCReg, xdecPCReg.sgn)
  atan2FlagReg.special := atan2Stage1Pre.io.special
  atan2FlagReg.ysgn    := ydecPCReg.sgn

//   printf("cir: atan2FlagReg.status  = %b\n", atan2FlagReg.status )
//   printf("cir: atan2FlagReg.special = %b\n", atan2FlagReg.special)
//   printf("cir: atan2FlagReg.ysgn    = %b\n", atan2FlagReg.ysgn   )

  // --------------------------------------------------------------------------

  val polynomialEval1 = Module(new PolynomialEval(spec, polySpec, maxCbits, stage1.calcStage))

  if(order != 0) {
    polynomialEval1.io.dx.get := ShiftRegister(recPre.io.dx.get, pcGap1)
  }
  polynomialEval1.io.coeffs.cs := recTab.io.cs.cs

  val polynomialResult1CPGapReg = ShiftRegister(polynomialEval1.io.result, cpGap1)

  atan2Stage1Post.io.en     := enCP1GapReg
  atan2Stage1Post.io.zother := ShiftRegister(atan2Stage1Other.io.zother, cpGap1)
  atan2Stage1Post.io.zres   := polynomialResult1CPGapReg
  atan2Stage1Post.io.minxy  := Mux(yIsLargerCPGapReg, xdecCPGapReg, ydecCPGapReg)

  val w = atan2Stage1Post.io.z

  // ---------------------------------------------------------------------------

  val wStage2Reg  = ShiftRegister(w, sGap)
  val enStage2Reg = ShiftRegister(enCP1GapReg, nPostStage1 + sGap)

  val enPC2Reg    = ShiftRegister(enStage2Reg, nPreStage2)
  val enPC2GapReg = ShiftRegister(enPC2Reg,    pcGap2)
  val enCP2Reg    = ShiftRegister(enPC2GapReg, nCalcStage2)
  val enCP2GapReg = ShiftRegister(enCP2Reg,    cpGap2)

  val wdecomp = Module(new DecomposeReal(spec))
  wdecomp.io.real := wStage2Reg

  val wdecPCReg    = ShiftRegister(wdecomp.io.decomp, nPreStage2)
  val wdecPCGapReg = ShiftRegister(wdecPCReg,         pcGap2)
  val wdecCPReg    = ShiftRegister(wdecPCGapReg,      nCalcStage2)
  val wdecCPGapReg = ShiftRegister(wdecCPReg,         cpGap2)

  val atan2Stage2Pre   = Module(new ATan2Stage2PreProcess (spec, polySpec, stage2.preStage))
  val atan2Stage2Tab   = Module(new ATan2Stage2TableCoeff (spec, polySpec, maxCbits))
  val atan2Stage2Other = Module(new ATan2Stage2OtherPath  (spec, polySpec, stage2.calcStage))
  val atan2Stage2Post  = Module(new ATan2Stage2PostProcess(spec, polySpec, stage2.postStage))

  atan2Stage2Pre.io.en  := enStage2Reg
  atan2Stage2Pre.io.x   := wdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  atan2Stage2Tab.io.en  := enPC2GapReg
  atan2Stage2Tab.io.adr := ShiftRegister(atan2Stage2Pre.io.adr, pcGap2)
  atan2Stage2Other.io.x := wdecPCGapReg

  atan2Stage2Other.io.flags := ShiftRegister(atan2FlagReg, nCalcStage1 + cpGap1 + nPostStage1 + sGap + nPreStage2)
  atan2Stage2Post.io.flags  := ShiftRegister(atan2FlagReg, nCalcStage1 + cpGap1 + nPostStage1 + sGap + nPreStage2)

  assert(atan2Stage2Pre.io.adr === 0.U               || enPC2Reg)
  assert(atan2Stage2Pre.io.dx.getOrElse(0.U) === 0.U || enPC2Reg)
  assert(atan2Stage2Tab.io.cs.asUInt === 0.U         || enPC2GapReg)

  val polynomialEval2 = Module(new PolynomialEval(spec, polySpec, maxCbits, stage2.calcStage))

  if(order != 0) {
    polynomialEval2.io.dx.get := ShiftRegister(atan2Stage2Pre.io.dx.get, pcGap2)
  }
  polynomialEval2.io.coeffs.cs := atan2Stage2Tab.io.cs.cs

  val polynomialResult2CPGapReg = ShiftRegister(polynomialEval2.io.result, cpGap2)

  atan2Stage2Post.io.en     := enCP2GapReg
  atan2Stage2Post.io.zother := ShiftRegister(atan2Stage2Other.io.zother, cpGap2)
  atan2Stage2Post.io.zres   := polynomialResult2CPGapReg

  io.z := atan2Stage2Post.io.z
}
