//% @file exp.scala
//
// exp function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._

import spire.math.SafeLong
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
import rial.arith.RoundSpec
import rial.arith.FloatChiselUtil

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

private[rial] object ExpPreMulArgs {
  def lhsW(spec: RealSpec): Int = {
    1 + spec.manW + 14 // XXX the width is determined empirically
  }

  def lhs(spec: RealSpec, w: Int): UInt = {
    val diffW = w - lhsW(spec)
    assert(diffW >= 0, f"PreProcMultiplier lhsW(${w}) should be wider or " +
      f"equal to lhsW(${lhsW(spec)})")

    val coef = Real.one / Real.log2
    if(diffW != 0) {
      Cat(coef(lhsW(spec)-1).toBigInt.U(lhsW(spec).W), 0.U(diffW.W))
    } else {
      coef(lhsW(spec)-1).toBigInt.U(lhsW(spec).W)
    }
  }

  def prodW(spec: RealSpec): Int = {
    (1 + spec.manW) + (1 + spec.manW + 14)
  }
  def prod(spec: RealSpec, out: UInt): UInt = {
    val outW = out.getWidth
    out(outW-1, outW - prodW(spec))
  }
}

private[rial] class ExpPreProcess(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW  = spec.exW
  val manW = spec.manW

  val exBias = spec.exBias

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order
  val extraBits = polySpec.extraBits

  val padding   = extraBits

  val io = IO(new Bundle {
    val en        = Input(UInt(1.W))
    val x         = Flipped(new DecomposedRealOutput(spec))
    val prod      = Input(UInt(ExpPreMulArgs.prodW(spec).W))
    val adr       = Output(UInt(adrW.W))
    val dx        = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }

    // integer part of x, to calculate zex
    val xint      = Output(UInt(exW.W))

    // LSB of fractional part of x, for zman correction
    val xfracLSBs = if(padding != 0) {Some(Output(UInt(padding.W)))} else {None}

    // 1 if (x * log2e).ex is larger than x.ex.
    val xexd = Output(UInt(1.W))
  })

  val xsgn  = enableIf(io.en, io.x.sgn)
  val xex0  = enableIf(io.en, io.x.ex )
  val xman0 = enableIf(io.en, io.x.man)

  val log2 = (a:Double) => {log(a) / log(2.0)}
  val xExOvfLimit = math.ceil(log2((maskL(exW)-exBias).toDouble)).toLong // log2(255-127 = 128) = 7
  val xExUdfLimit = math.ceil(log2((abs(0 - exBias)  ).toDouble)).toLong // log2(|0-127| = 127) > 6

  val xIntW  = max(xExOvfLimit, xExUdfLimit).toInt
  val xFracW = manW + padding
  val xValW  = xIntW + xFracW + 1

  val extraMan = padding + xIntW

  // ------------------------------------------------------------------------
  // for exponential, multiply x by log2e.

  val log2eManW = ExpPreMulArgs.lhsW(spec) - 1
  val xprod  = io.prod
  val xprodW = ExpPreMulArgs.prodW(spec)

  val xprodMoreThan2 = xprod(xprodW - 1) === 1.U

  assert(log2eManW - extraMan > 0)

  val xprodbp        = spec.manW + log2eManW
  val xprodRoundBits = xprodbp - (spec.manW + extraMan)

  val xprodSticky  = (xprodMoreThan2 & xprod(xprodRoundBits-1)) |
                     xprod(xprodRoundBits-2, 0).orR
  val xprodRound   = Mux(xprodMoreThan2, xprod(xprodRoundBits), xprod(xprodRoundBits-1))
  val xprodShifted = Mux(xprodMoreThan2, xprod(xprodbp, xprodRoundBits+1),
                                         xprod(xprodbp-1, xprodRoundBits))
  val xprodLSB     = xprodShifted(0)

  val xprodInc = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
    xprodLSB, xprodRound, xprodSticky)

  val xprodRounded = xprodShifted +& xprodInc
  assert(xprodRounded.getWidth == manW+1 + extraMan)
  val xprodMoreThan2AfterRounded = xprodRounded(manW+extraMan)

  val xman = xprodRounded(manW+extraMan-1, 0)
  val xex  = xex0 + xprodMoreThan2AfterRounded + xprodMoreThan2

  assert(xprodMoreThan2AfterRounded +& xprodMoreThan2 =/= 2.U)

  val xexd0 = enableIf(io.en, xprodMoreThan2AfterRounded + xprodMoreThan2)
  io.xexd := ShiftRegister(xexd0, nStage)

  val xVal = Cat(1.U(1.W), xman)
  assert(xVal.getWidth == xValW)

  // ------------------------------------------------------------------------
  // do preprocess

  val xshift0 = (xIntW + exBias).U(exW.W) - xex
  val xshift  = xshift0(xIntW-1, 0)
  val xValShifted = Mux(xshift0 >= (xValW+1).U, 0.U, xVal >> (xshift-1.U(1.W)))
  val xValRounded = (xValShifted >> 1) + xValShifted(0)

  val xint0  = xValRounded(xIntW+xFracW-1, xFracW)
  val xfrac0 = xValRounded(xFracW-1, 0)

  // if xsgn == 1, we negate xint to calculate exbias. but we later do that in
  // ExpOtherPath.
  val xint = enableIf(io.en, (xint0 +& (xsgn.asBool && (xfrac0 =/= 0.U)).asUInt))
  io.xint := ShiftRegister(xint, nStage)

  val xfracNeg = (BigInt(1)<<xFracW).U((xFracW+1).W) - xfrac0
  val xfrac = Mux(xsgn === 0.U, xfrac0, xfracNeg(xFracW-1, 0))

  val adr  = enableIf(io.en, xfrac(xFracW-1, (xFracW-1)-adrW+1))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx   = enableIf(io.en, Cat(~xfrac(xFracW-1-adrW), xfrac(xFracW-1-adrW-1, padding)))
    io.dx.get := ShiftRegister(dx, nStage)
  }

  if(padding != 0) {
    val xfracLSBs  = enableIf(io.en, xfrac(padding-1, 0))
    io.xfracLSBs.get := ShiftRegister(xfracLSBs, nStage)
  }
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

private[rial] class ExpTableCoeff(
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
    val adr = Input  (UInt((1+adrW).W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {

    val tableI = ExpSim.expTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/1)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    assert(maxCbit(0) == fracW)
    assert(coeff(0).getWidth == fracW)

    io.cs.cs(0) := enableIf(io.en, coeff(0))

  } else {
    // x -> xInt + xFrac
    // pow2(x) -> 2^xInt * 2^xFrac
    // Since xFrac is in [0, 1), 2^xFrac is in [1, 2)

    val tableI = ExpSim.expTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    // both 1st and 2nd derivative of 2^x is larger than 0
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

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
private[rial] object ExpTableCoeff {
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
      return ExpSim.expTableGeneration( order, adrW, spec.manW, fracW ).cbit
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
      return ExpSim.expTableGeneration( order, adrW, spec.manW, fracW ).calcWidth
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

private[rial] class ExpNonTableOutput(val spec: RealSpec) extends Bundle {
  //  zsgn  = 0, always
  val zex   = Output(UInt(spec.exW.W))
  val znan  = Output(Bool())
  val zinf  = Output(Bool())
  val zzero = Output(Bool())
  val zIsNonTable = Output(Bool())
}

// No pathway (like taylor series) other than table interpolation.
// Here we calculate z.ex and correction to zman.
private[rial] class ExpOtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val extraBits = polySpec.extraBits
  val padding = extraBits

  val zCorrTmpW = min(10, 1+manW+padding) // XXX see expSim

  val io = IO(new Bundle {
    val x         = Flipped(new DecomposedRealOutput(spec))
    val xint      = Input(UInt(exW.W))
    val xexd      = Input(UInt(1.W)) // if exp & (x * log2e).ex is larger than x.ex, 1
    // for z correction
    val xfracLSBs = if(padding != 0) {Some(Input (UInt(padding.W  )))} else {None}
    val zCorrCoef = if(padding != 0) {Some(Output(UInt(zCorrTmpW.W)))} else {None}
    // normal output
    val zother    = new ExpNonTableOutput(spec)
  })

  // --------------------------------------------------------------------------
  // calc first half of zCorrection = z * xfracLSBs * ln2
  //                                      ^^^^^^^^^^^^^^^ this part
  //
  // note: z isn't calculated yet, so here we can only calculate xfracLSBs * ln2
  if(padding != 0) {
    val coefficient = (Real.log(Real.two))(spec.manW+1).toBigInt.U((1+manW).W)

    val xfracLSBs  = io.xfracLSBs.get
    val zCorrCoef0 = coefficient * xfracLSBs
    val zCorrCoefW = 1 + manW + padding
    val zCorrCoef  = zCorrCoef0(zCorrCoefW-1, zCorrCoefW-zCorrTmpW)

    assert(zCorrCoef.getWidth  == zCorrTmpW)
    assert(xfracLSBs.getWidth  == padding)
    assert(zCorrCoef0.getWidth == zCorrCoefW)

    io.zCorrCoef.get := ShiftRegister(zCorrCoef, nStage)
  }

  // --------------------------------------------------------------------------
  // calc zex

  val log2 = (a:Double) => {log(a) / log(2.0)}
  val xExOvfLimit = math.ceil(log2( spec.exMax)).toInt // ceil(log2(128)) = 7
  val xExUdfLimit = math.ceil(log2(-spec.exMin)).toInt // ceil(log2(127)) = 7

  val xExOvfLimBiased = (xExOvfLimit + exBias).U(exW.W)
  val xExUdfLimBiased = (xExUdfLimit + exBias).U(exW.W)

  val zexPos = io.xint +& exBias.U(exW.W)
  val zexNeg = exBias.U(exW.W) -& io.xint
  val zex0   = Mux(io.x.sgn === 0.U, zexPos(exW-1, 0), zexNeg(exW-1, 0))

  val xex = io.x.ex + io.xexd

  val znan  = io.x.nan
  val zinf  = (io.x.sgn === 0.U) &&
    ((xExOvfLimBiased <= xex) || (  spec.exMax.U  < io.xint) || io.x.inf)
  val zzero = (io.x.sgn === 1.U) &&
    ((xExUdfLimBiased <= xex) || ((-spec.exMin).U < io.xint) || io.x.inf)

  val zex = Mux(znan || zinf, Fill(exW, 1.U(1.W)),
            Mux(zzero, 0.U(exW.W),
                zex0))

  val zIsNonTable = znan || zinf || zzero

  io.zother.zex         := ShiftRegister(zex,   nStage)
  io.zother.znan        := ShiftRegister(znan,  nStage)
  io.zother.zinf        := ShiftRegister(zinf,  nStage)
  io.zother.zzero       := ShiftRegister(zzero, nStage)
  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

private[rial] class ExpPostProcess(
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

  val padding = extraBits

  val zCorrTmpW = min(10, 1+manW+padding) // XXX see expSim

  val io = IO(new Bundle {
    val en = Input(UInt(1.W))
    // ex and some flags
    val zother = Flipped(new ExpNonTableOutput(spec))
    // table interpolation results
    val zres   = Input(UInt(fracW.W))
    // for z correction
    val zCorrCoef = if(padding != 0) {Some(Input(UInt(zCorrTmpW.W)))} else {None}
    // output
    val z      = Output(UInt(spec.W.W))
  })

  val zex0  = io.zother.zex
  val zman0 = if(extraBits > 0) {
    dropLSB(extraBits, io.zres) + io.zres(extraBits-1)
  } else {
    io.zres
  }

  // --------------------------------------------------------------------------
  // calc zman correction coeff

  val zCorrection = Wire(UInt(1.W))

  if(padding != 0) {
    val zCorrCoef  = io.zCorrCoef.get
    val zCorrTerm0 = Cat(1.U(1.W), zman0) * zCorrCoef
    val zCorrTermW = 1 + manW + zCorrTmpW
    val zCorrTerm  = zCorrTerm0(zCorrTermW-1, zCorrTermW-zCorrTmpW)
    assert(zCorrTerm.getWidth == zCorrTmpW)

    zCorrection := (zCorrTerm(zCorrTmpW-1).asBool ||
                    zCorrTerm(zCorrTmpW-2).asBool).asUInt
  } else {
    // no correction needed.
    zCorrection := 0.U(1.W)
  }

  // --------------------------------------------------------------------------
  // calc z

  val zmanCorrected = zman0 +& zCorrection
  assert(zmanCorrected.getWidth == manW+1)

  val zexCorrection = zmanCorrected(manW) & (!io.zother.zIsNonTable)

  val zEx  = zex0 + zexCorrection
  val zMan = Mux(io.zother.zIsNonTable, Cat(io.zother.znan, 0.U((manW-1).W)),
                 zmanCorrected(manW-1, 0))
  val z0   = Cat(0.U(1.W), zEx, zMan)

  assert(zEx .getWidth == exW)
  assert(zMan.getWidth == manW)
  assert(z0  .getWidth == spec.W)

  val z = enableIf(io.en, z0)

  io.z   := ShiftRegister(z, nStage)
}
