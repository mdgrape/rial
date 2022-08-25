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

// ATan2 Phase1 calculates min(x,y)/max(x,y).
//       Phase2 calculates atan(min(x,y)/max(x,y)) +/- constant.
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
// Phase1 preprocess only checks the special cases.
// Phase1 calculates min(x,y)/max(x,y), so ReciprocalPreProcess is re-used.
//
class ATan2Phase1PreProcess(
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
    val en        = Input (UInt(1.W))
    val x         = Flipped(new DecomposedRealOutput(spec))
    val y         = Flipped(new DecomposedRealOutput(spec))
    val yIsLarger = Input(Bool())
    val special   = Output(UInt(ATan2SpecialValue.W.W))
  })

  val minex     = Mux(io.yIsLarger, io.x.ex, io.y.ex)
  val maxex     = Mux(io.yIsLarger, io.y.ex, io.x.ex)
  val diffexDec = Mux(io.yIsLarger, io.y.man > io.x.man, io.x.man > io.y.man)
  val zeroed    = minex +& exBias.U(exW.W) <= maxex +& diffexDec.asUInt

  // |min(x,y)| / |max(x,y)| = 0 means atan2(y,x) = (n/2)pi, n=0,1,2,3
  // case |y| << |x| && 0 < x : z = 0
  // case |y| << |x| && x < 0 : z = pi
  // case |x| << |y| && 0 < y : z = pi/2
  // case |x| << |y| && y < 0 : z = 3pi/2

  val tooLargeX = (zeroed && !io.yIsLarger) || ( io.x.inf && !io.y.inf) || (!io.x.zero &&  io.y.zero)
  val tooLargeY = (zeroed &&  io.yIsLarger) || (!io.x.inf &&  io.y.inf) || ( io.x.zero && !io.y.zero)

  val samexy    = io.x.ex === io.y.ex && io.x.man === io.y.man && !io.x.nan && !io.y.nan

  val xpos      = io.x.sgn === 0.U
  val xneg      = io.x.sgn === 1.U
  val zzero     = tooLargeX && xpos
  val zpi       = tooLargeX && xneg
  val zhalfpi   = tooLargeY
  val z1piover4 = samexy && xpos
  val z3piover4 = samexy && xneg
  val znan      = (io.x.nan || io.y.nan) || (io.x.zero && io.y.zero)

  val special0 = MuxCase(ATan2SpecialValue.zNormal, Seq(
    znan      -> ATan2SpecialValue.zNaN,
    zzero     -> ATan2SpecialValue.zZero,
    zpi       -> ATan2SpecialValue.zPi,
    zhalfpi   -> ATan2SpecialValue.zHalfPi,
    z1piover4 -> ATan2SpecialValue.zQuarterPi,
    z3piover4 -> ATan2SpecialValue.z3QuarterPi
  ))
  val special = enable(io.en, special0)

  assert(!io.en || !zeroed || special0 =/= ATan2SpecialValue.zNormal)

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

class ATan2Phase1NonTableOutput(val spec: RealSpec) extends Bundle {
  val zex       = Output(UInt(spec.exW.W)) // sign of min(x,y)/max(x,y)
  val maxXYMan0 = Output(Bool())           // max(x,y).man === 0.U
  val xySameMan = Output(Bool())
}

class ATan2Phase1OtherPath(
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

    val zother  = new ATan2Phase1NonTableOutput(spec)
  })

  val maxXYMan0 = Mux(io.yIsLarger, (!io.y.man.orR.asBool), (!io.x.man.orR.asBool))
//   printf("x.man = %b\n", io.x.man)
//   printf("y.man = %b\n", io.y.man)
//   printf("x < y = %b\n", io.yIsLarger)

  val xySameMan = io.x.man === io.y.man

  io.zother.maxXYMan0 := ShiftRegister(maxXYMan0, nStage)
  io.zother.xySameMan := ShiftRegister(xySameMan, nStage)

  val exDec  = Mux(Mux(io.yIsLarger, io.x.man < io.y.man, io.y.man < io.x.man),
                   1.U(1.W), 0.U(1.W))
  val maxEx  = Mux(io.yIsLarger, io.y.ex, io.x.ex)
  val minEx  = Mux(io.yIsLarger, io.x.ex, io.y.ex)
  val zeroed = minEx +& exBias.U(exW.W) <= maxEx + exDec.asUInt
  val zex0  = Mux(io.x.inf && io.y.inf, exBias.U(exW),
              Mux((io.x.inf && !io.y.inf) || (!io.x.inf && io.y.inf), 0.U(exW.W),
                  ((minEx +& exBias.U) - maxEx) - exDec.asUInt))
  val zex   = Mux(zeroed, 0.U, zex0)

  // exponent of min(x,y)/max(x,y). we will later correct +/- 1 by checking
  // the mantissa of min(x,y) and max(x,y)
  io.zother.zex := ShiftRegister(zex, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

// Phase1: takes 1/max(x, y) and min(x, y), returns min(x, y) / max(x, y)

class ATan2Phase1PostProcess(
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
    val en     = Input(Bool())
    val zother = Flipped(new ATan2Phase1NonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val minxy  = Flipped(new DecomposedRealOutput(spec))
    val z      = Output(UInt(spec.W.W))
  })

  val zsgn      = 0.U(1.W)
  val zex       = io.zother.zex
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
// Phase2 calculates atan(x), so we need to extract the address value

class ATan2Phase2PreProcess(
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
// Phase1 re-use the reciprocal table. we don't need to implement it for atan2.
//
class ATan2Phase2TableCoeff(
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

    val tableI = ATan2Phase2Sim.atanTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    // sign mode 0: always include sign
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode = */1)
    val coeff = getSlices(coeffTable(adr), coeffWidth)

    io.cs.cs(0) := enable(io.en, coeff(0))

  } else {

    val tableI = ATan2Phase2Sim.atanTableGeneration( order, adrW, manW, fracW )
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

object ATan2Phase2TableCoeff {
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
      return ATan2Phase2Sim.atanTableGeneration( order, adrW, spec.manW, fracW ).
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
      return ATan2Phase2Sim.atanTableGeneration( order, adrW, spec.manW, fracW ).calcWidth
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

class ATan2Phase2NonTableOutput(val spec: RealSpec) extends Bundle {
  val zsgn        = Output(UInt(1.W))
  val zex         = Output(UInt(spec.exW.W))
  val zman        = Output(UInt(spec.manW.W))
  val zIsNonTable = Output(Bool())
}

class ATan2Phase2OtherPath(
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
    val zother  = new ATan2Phase2NonTableOutput(spec)
  })

  val pi         = new RealGeneric(spec, Pi)
  val halfPi     = new RealGeneric(spec, Pi * 0.5)
  val quarterPi  = new RealGeneric(spec, Pi * 0.25)
  val quarter3Pi = new RealGeneric(spec, Pi * 0.75)

  val linearThreshold = (ATan2Phase2Sim.calcLinearThreshold(manW) + exBias)
  val isLinear = io.x.ex <= linearThreshold.U(exW.W)
  // non-linear mantissa is calculated by table.
//   printf("cir: Phase2OtherPath: linearThr = %b\n", linearThreshold.U)
//   printf("cir: Phase2OtherPath: isLinear  = %b\n", isLinear)
//   printf("cir: Phase2OtherPath: isspecial = %b\n", io.flags.special)

  io.zother.zIsNonTable := ShiftRegister(
    isLinear || (io.flags.special =/= ATan2SpecialValue.zNormal), nStage)

  val defaultEx  = io.x.ex // isLinear includes xzero.
  val defaultMan = io.x.man
//   printf("cir: Phase2OtherPath: defaultEx  = %b\n", defaultEx )
//   printf("cir: Phase2OtherPath: defaultMan = %b\n", defaultMan)

  io.zother.zsgn := ShiftRegister(io.flags.ysgn, nStage)
  io.zother.zex  := ShiftRegister(MuxCase(defaultEx, Seq(
    (io.flags.special === ATan2SpecialValue.zNaN)        -> maskL(exW).U(exW.W),
    (io.flags.special === ATan2SpecialValue.zZero)       -> 0.U(exW.W),
    (io.flags.special === ATan2SpecialValue.zPi)         -> pi.ex.U(exW.W),
    (io.flags.special === ATan2SpecialValue.zHalfPi)     -> halfPi.ex.U(exW.W),
    (io.flags.special === ATan2SpecialValue.zQuarterPi)  -> quarterPi.ex.U(exW.W),
    (io.flags.special === ATan2SpecialValue.z3QuarterPi) -> quarter3Pi.ex.U(exW.W)
    )), nStage)
  io.zother.zman := ShiftRegister(MuxCase(defaultMan, Seq(
    (io.flags.special === ATan2SpecialValue.zNaN)        -> Cat(1.U(1.W), 0.U((manW-1).W)),
    (io.flags.special === ATan2SpecialValue.zZero)       -> 0.U(manW.W),
    (io.flags.special === ATan2SpecialValue.zPi)         -> pi        .man.toBigInt.U(manW.W),
    (io.flags.special === ATan2SpecialValue.zHalfPi)     -> halfPi    .man.toBigInt.U(manW.W),
    (io.flags.special === ATan2SpecialValue.zQuarterPi)  -> quarterPi .man.toBigInt.U(manW.W),
    (io.flags.special === ATan2SpecialValue.z3QuarterPi) -> quarter3Pi.man.toBigInt.U(manW.W)
    )), nStage)
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
// Phase2: takes status flags, returns corrected result
//

class ATan2Phase2PostProcess(
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
    val zother = Flipped(new ATan2Phase2NonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val flags  = Input(new ATan2Flags()) // need status (|x|<|y|, xsgn)
    val z      = Output(UInt(spec.W.W))
  })

  val zSgn = io.zother.zsgn

  val zres = Cat(io.zres, 0.U(1.W))
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

  // XXX there is no overflow.
  // 1. io.x in atan2stage2 is the result of min(x,y)/max(x,y), so |io.x| <= 1.
  //   - it means io.x.ex <= 127. Since exW=8, atanEx0 never become larger than 2^exW
  // 2. we don't use subnormal here.
  //   - it means that io.x.ex >= 1 or io.x == 0.
  //   - if io.x is zero, it means the result is one of the special values.
  //     overflow can happen, but is ignorable.
  //   - if io.x.ex >= 1, atanEx0 does not overflow.
  assert(!io.en || io.x.ex =/= 0.U || io.flags.special =/= 0.U)
  assert(!io.en || io.flags.special =/= 0.U || atanEx0(exW) === 0.U)

  val atanEx  = Mux(io.zother.zIsNonTable, io.zother.zex,  atanEx0(exW-1, 0))
  val atanMan = Mux(io.zother.zIsNonTable, io.zother.zman, atanMan0)
//   printf("cir: isNonTab= %d\n", io.zother.zIsNonTable)
//   printf("cir: atanEx  = %d\n", atanEx)
//   printf("cir: atanMan = %b\n", atanMan)

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

  assert(zExPL.getWidth == spec.exW)
  assert(zExNL.getWidth == spec.exW)
  assert(zExPS.getWidth == spec.exW)
  assert(zExNS.getWidth == spec.exW)
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
  assert(zSgn.getWidth == 1)
  assert(zEx .getWidth == spec.exW)
  assert(zMan.getWidth == spec.manW)
  assert(zNormal.getWidth == spec.W)

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
