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

// ATan2 Phase1 calculates min(x,y)/max(x,y).
//       Phase2 calculates atan(min(x,y)/max(x,y)) +/- constant.

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

    val xdecoded = new DecomposedRealOutput(spec)
    val flags = Output(new ATan2Flags)
  })

  // stage1 calculates min(x,y)/max(x,y), so we can assume that x <= 1.
  // also, if x == 1, stage1 sets special value flag. so we can ignore the case.

  // --------------------------------------------------------------------------
  // decode flags and x

  val flags = WireDefault(0.U.asTypeOf(new ATan2Flags))

  flags.ysgn      := io.x.sgn
  flags.status    := Cat(io.x.ex(exW-1), io.x.man(0))
  flags.isSpecial := io.x.ex === ((1 << exW) - 1).U
  flags.special   := io.x.man

  val xdecoded = WireDefault(0.U.asTypeOf(new Bundle {
    val sgn  = UInt(1.W)
    val ex   = UInt(spec.exW.W)
    val man  = UInt(spec.manW.W)
    val zero = Bool()
    val inf  = Bool()
    val nan  = Bool()
  }))

  xdecoded.sgn  := 0.U
  xdecoded.ex   := Cat(0.U(1.W), io.x.ex(exW-2, 0))
  xdecoded.man  := Cat(io.x.man(manW-1, 1), 0.U(1.W))
  xdecoded.zero := 0.U
  xdecoded.inf  := 0.U
  xdecoded.nan  := 0.U

  io.flags    := ShiftRegister(enableIf(io.en, flags), nStage)
  io.xdecoded := ShiftRegister(enableIf(io.en, xdecoded), nStage)

  // --------------------------------------------------------------------------
  // remove flags from exponent and mantissa

  val xex    = Cat(0.U(1.W), io.x.ex(exW-2, 0))
  val xmanW1 = Cat(1.U(1.W), io.x.man(manW-1, 1), 0.U(1.W))

  val shiftW    = log2Up(1+manW)
  val xShift0   = exBias.U(exW.W) - xex
  val xShiftOut = xShift0(exW-1, shiftW).orR
  val xShift    = Mux(xShiftOut, maskL(shiftW).U(shiftW.W), xShift0(shiftW-1, 0))
  val xFixed    = xmanW1 >> xShift

  // xFixed(manW) === 1 if and only if x == 1.0.
  assert(xFixed(manW) =/= 1.U || xFixed(manW-1, 0) === 0.U || io.en === 0.U,
         "xFixed = %b", xFixed)

  val adr = xFixed(manW-1, dxW)
  io.adr := ShiftRegister(enableIf(io.en, adr), nStage)

  if(order != 0) {
    val dx = Cat(~xFixed(manW-adrW-1), xFixed(manW-adrW-2,0))
    io.dx.get := ShiftRegister(enableIf(io.en, dx), nStage)
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

    io.cs.cs(0) := enableIf(io.en, coeff(0))

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
    io.cs := enableIf(io.en, coeffs)
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

  io.zother.zIsNonTable := ShiftRegister(isLinear || io.flags.isSpecial, nStage)

  val defaultEx  = io.x.ex // isLinear includes xzero.
  val defaultMan = io.x.man
//   printf("cir: Phase2OtherPath: defaultEx  = %b\n", defaultEx )
//   printf("cir: Phase2OtherPath: defaultMan = %b\n", defaultMan)

  io.zother.zsgn := ShiftRegister(io.flags.ysgn, nStage)
  io.zother.zex  := ShiftRegister(MuxCase(defaultEx, Seq(
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zNaN)        -> maskL(exW).U(exW.W),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zZero)       -> 0.U(exW.W),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zPi)         -> pi.ex.U(exW.W),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zHalfPi)     -> halfPi.ex.U(exW.W),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zQuarterPi)  -> quarterPi.ex.U(exW.W),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.z3QuarterPi) -> quarter3Pi.ex.U(exW.W)
    )), nStage)
  io.zother.zman := ShiftRegister(MuxCase(defaultMan, Seq(
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zNaN)        -> Cat(1.U(1.W), 0.U((manW-1).W)),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zZero)       -> 0.U(manW.W),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zPi)         -> pi        .man.toBigInt.U(manW.W),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zHalfPi)     -> halfPi    .man.toBigInt.U(manW.W),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.zQuarterPi)  -> quarterPi .man.toBigInt.U(manW.W),
    (io.flags.isSpecial && io.flags.special === ATan2SpecialValue.z3QuarterPi) -> quarter3Pi.man.toBigInt.U(manW.W)
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
    val zman0  = Input(UInt(manW.W))
    val zexInc = Input(UInt(1.W))
    val flags  = Input(new ATan2Flags()) // need status (|x|<|y|, xsgn)
    val z      = Output(UInt(spec.W.W))
  })

  val zSgn = io.zother.zsgn

  val atanMan0 = io.zman0
  val atanEx0  = io.x.ex +& io.zexInc - 1.U

  // XXX there is no overflow.
  // 1. io.x in atan2stage2 is the result of min(x,y)/max(x,y), so |io.x| <= 1.
  //   - it means io.x.ex <= 127. Since exW=8, atanEx0 never become larger than 2^exW
  // 2. we don't use subnormal here.
  //   - it means that io.x.ex >= 1 or io.x == 0.
  //   - if io.x is zero, it means the result is one of the special values.
  //     overflow can happen, but is ignorable.
  //   - if io.x.ex >= 1, atanEx0 does not overflow.
  assert(!io.en || io.x.ex =/= 0.U || io.flags.isSpecial)
  assert(!io.en || io.flags.isSpecial || atanEx0(exW) === 0.U)

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

  val special = Mux(io.flags.isSpecial, io.flags.special, ATan2SpecialValue.zNormal)
  val z0 = MuxCase(zNormal, Seq(
      (special === ATan2SpecialValue.zNormal    ) -> zNormal,
      (special === ATan2SpecialValue.zNaN       ) -> Cat(0.U(1.W), nan   .ex.U(exW.W), nan       .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.zZero      ) -> Cat(zSgn, zero      .ex.U(exW.W), zero      .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.zPi        ) -> Cat(zSgn, pi        .ex.U(exW.W), pi        .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.zHalfPi    ) -> Cat(zSgn, halfPi    .ex.U(exW.W), halfPi    .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.zQuarterPi ) -> Cat(zSgn, quarterPi .ex.U(exW.W), quarterPi .man.toBigInt.U(manW.W)),
      (special === ATan2SpecialValue.z3QuarterPi) -> Cat(zSgn, quarter3Pi.ex.U(exW.W), quarter3Pi.man.toBigInt.U(manW.W))
    ))
  val z = enableIf(io.en, z0)

  io.z := ShiftRegister(z, nStage)
}
