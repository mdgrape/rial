//% @file Exp.scala
//
// exp(x) hardware modules. Pre-process consumes the shared PreProcMultiplier
// to compute x * log2(e); the result is split into integer/fractional parts.
// The fractional part feeds the polynomial table (which approximates 2^t-1
// in [0, 1)), and the integer part becomes the result's exponent.
//
// A small ln(2)*xfracLSBs correction term is computed in ExpNonTable and
// applied in ExpPostProcess to recover the precision lost when rounding the
// pre-product to the table address width.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._

import spire.math.Real
import spire.implicits._

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{RealSpec, RoundSpec, FloatChiselUtil}
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._

// -------------------------------------------------------------------------
//
// constants for the pre-process multiplier (`log2(e)`)
//

object ExpPreMulArgs {

  // Minimum lhs width the pre-mult must have for exp.
  // The trailing pad (14 bits) was determined empirically
  def lhsW(spec: RealSpec): Int = { 1 + spec.manW + 14 }

  def lhs(spec: RealSpec, w: Int): UInt = {
    val base  = lhsW(spec)
    val diffW = w - base
    assert(diffW >= 0, f"PreMul lhsW($w) must be >= exp lhsW(${base})")

    val coef = Real.one / Real.log2
    val v    = coef(base - 1).toBigInt

    if (diffW == 0) {
      v.U(base.W)
    } else {
      Cat(v.U(base.W), 0.U(diffW.W))
    }
  }

  def prodW(spec: RealSpec): Int = {
    (1 + spec.manW) + lhsW(spec)
  }
  def prod(spec: RealSpec, out: UInt): UInt = {
    val w = out.getWidth
    out(w - 1, w - prodW(spec))
  }
}

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

/** Information that flows from [[ExpPreProcess]] to [[ExpNonTable]]. */
class ExpPreProcInfo(val spec: RealSpec, val padding: Int) extends Bundle {
  val xint      = UInt(spec.exW.W)
  val xexd      = UInt(1.W)
  val xfracLSBs = if (padding != 0) Some(UInt(padding.W)) else None
}

class ExpPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW    = bitW.addrW
  val dxW     = bitW.inBitW
  val order   = bitW.nOrder
  val padding = bitW.outBitW - manW

  // width budget
  private val log2 = (a: Double) => log(a) / log(2.0)

  val xExOvfLimit = scala.math.ceil(log2((maskL(exW) - exBias).toDouble)).toLong
  val xExUdfLimit = scala.math.ceil(log2((abs(0 - exBias)).toDouble)).toLong
  val xIntW       = max(xExOvfLimit, xExUdfLimit).toInt
  val xFracW      = manW + padding
  val xValW       = xIntW + xFracW + 1
  val extraMan    = padding + xIntW

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val prod = Flipped(Valid(UInt(ExpPreMulArgs.prodW(spec).W)))
    val adr  = Valid(UInt(adrW.W))
    val dx   = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    val info = Valid(new ExpPreProcInfo(spec, padding))
  })

  val xsgn = io.x.bits.sgn
  val xex0 = io.x.bits.ex
  val xman = io.x.bits.man

  // ------------------------------------------------------------------------
  // round the pre-product to (manW + extraMan) bits

  val log2eManW      = ExpPreMulArgs.lhsW(spec) - 1
  val xprod          = io.prod.bits
  val xprodW         = ExpPreMulArgs.prodW(spec)
  val xprodMoreThan2 = xprod(xprodW - 1) === 1.U

  assert(log2eManW - extraMan > 0)

  val xprodbp        = manW + log2eManW
  val xprodRoundBits = xprodbp - (manW + extraMan)

  val xprodSticky  = (xprodMoreThan2 & xprod(xprodRoundBits - 1)) |
                     xprod(xprodRoundBits - 2, 0).orR
  val xprodRound   = Mux(xprodMoreThan2, xprod(xprodRoundBits), xprod(xprodRoundBits - 1))
  val xprodShifted = Mux(xprodMoreThan2, xprod(xprodbp,     xprodRoundBits + 1),
                                         xprod(xprodbp - 1, xprodRoundBits))
  val xprodLSB     = xprodShifted(0)
  val xprodInc     = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
                                                    xprodLSB, xprodRound, xprodSticky)
  val xprodRounded = xprodShifted +& xprodInc
  assert(xprodRounded.getWidth == manW + 1 + extraMan)
  val xprodMoreThan2AfterRounded = xprodRounded(manW + extraMan)

  val xman2 = xprodRounded(manW + extraMan - 1, 0)
  val xex   = xex0 + xprodMoreThan2AfterRounded + xprodMoreThan2

  assert(!io.prod.valid || xprodMoreThan2AfterRounded +& xprodMoreThan2 =/= 2.U)

  val xexd = xprodMoreThan2AfterRounded + xprodMoreThan2
  val xVal = Cat(1.U(1.W), xman2)
  assert(xVal.getWidth == xValW)

  // ------------------------------------------------------------------------
  // split into integer/fractional parts

  val xshift0     = (xIntW + exBias).U(exW.W) - xex
  val xshift      = xshift0(xIntW - 1, 0)
  val xValShifted = Mux(xshift0 >= (xValW + 1).U, 0.U, xVal >> (xshift - 1.U(1.W)))
  val xValRounded = (xValShifted >> 1) + xValShifted(0)

  val xint0  = xValRounded(xIntW + xFracW - 1, xFracW)
  val xfrac0 = xValRounded(xFracW - 1, 0)

  // (sign handling: ExpNonTable converts xint to a biased exponent.)
  val xint = xint0 +& (xsgn.asBool && (xfrac0 =/= 0.U)).asUInt

  val xfracNeg = (BigInt(1) << xFracW).U((xFracW + 1).W) - xfrac0
  val xfrac    = Mux(xsgn === 0.U, xfrac0, xfracNeg(xFracW - 1, 0))

  // ------------------------------------------------------------------------
  // outputs

  val adr = xfrac(xFracW - 1, (xFracW - 1) - adrW + 1)

  val infoOut = Wire(new ExpPreProcInfo(spec, padding))
  infoOut.xint := xint
  infoOut.xexd := xexd
  if (padding != 0) {
    infoOut.xfracLSBs.get := xfrac(padding - 1, 0)
  }

  val valid = io.x.valid && io.prod.valid
  io.info := Pipe(valid, infoOut, stage)
  io.adr  := Pipe(valid, adr,     stage)
  if (order != 0) {
    val dx = Cat(~xfrac(xFracW - 1 - adrW), xfrac(xFracW - 1 - adrW - 1, padding))
    io.dx.get := Pipe(valid, dx, stage)
  }
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class ExpTableCoeff(
  val spec:  RealSpec,
  val coeff: PolynomialCoefficients,
  val stage: Int,
) extends Module {

  val adrW = coeff.bitW.addrW

  val io = IO(new Bundle {
    val adr = Flipped(Valid(UInt(adrW.W)))
    val cs  = Valid(UInt(coeff.bitW.coefBits.reduce(_+_).W))
  })

  val packed = VecInit(PolynomialGenerator.packTable(coeff).map(c => {
    c.U(coeff.bitW.coefBits.reduce(_+_).W)
  }))
  io.cs := Pipe(io.adr.valid, packed(io.adr.bits), stage)
}

// -------------------------------------------------------------------------
//                        _        _     _                    _   _
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_|
//                                               |_|
// -------------------------------------------------------------------------

class ExpPostProcInfo(val spec: RealSpec, val padding: Int) extends Bundle {
  val zex         = UInt(spec.exW.W)
  val znan        = Bool()
  val zinf        = Bool()
  val zzero       = Bool()
  val zIsNonTable = Bool()
  val zCorrCoef   = if (padding != 0) Some(UInt(scala.math.min(6, 1 + spec.manW + padding).W)) else None
}

class ExpNonTable(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val padding = bitW.outBitW - manW

  val zCorrTmpW = scala.math.min(6, 1 + manW + padding)
  val coeffW    = scala.math.min(8, manW + 1)

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val info = Flipped(Valid(new ExpPreProcInfo(spec, padding)))
    val z    = Valid(new ExpPostProcInfo(spec, padding))
  })

  // ------------------------------------------------------------------------
  // small inline multiplier for the correction coefficient ln(2) * xfracLSBs

  val zout = Wire(new ExpPostProcInfo(spec, padding))

  if (padding != 0) {
    val coefficient = (Real.log(Real.two))(coeffW).toBigInt.U(coeffW.W)
    val xfracLSBs   = io.info.bits.xfracLSBs.get
    val zCorrCoef0  = coefficient * xfracLSBs
    val zCorrCoefW  = coeffW + padding
    val zCorrCoef   = zCorrCoef0(zCorrCoefW - 1, zCorrCoefW - zCorrTmpW)
    assert(zCorrCoef.getWidth == zCorrTmpW)
    zout.zCorrCoef.get := zCorrCoef
  }

  // ------------------------------------------------------------------------
  // exponent + special-flag computation

  private val log2 = (a: Double) => log(a) / log(2.0)
  val xExOvfLimit = scala.math.ceil(log2(spec.exMax.toDouble)).toInt
  val xExUdfLimit = scala.math.ceil(log2((-spec.exMin).toDouble)).toInt

  val xExOvfLimBiased = (xExOvfLimit + exBias).U(exW.W)
  val xExUdfLimBiased = (xExUdfLimit + exBias).U(exW.W)

  val xint = io.info.bits.xint
  val xexd = io.info.bits.xexd

  val zexPos = xint +& exBias.U(exW.W)
  val zexNeg = exBias.U(exW.W) -& xint
  val zex0   = Mux(io.x.bits.sgn === 0.U, zexPos(exW - 1, 0), zexNeg(exW - 1, 0))

  val xex   = io.x.bits.ex + xexd
  val znan  = io.x.bits.nan
  val zinf  = (io.x.bits.sgn === 0.U) &&
              ((xExOvfLimBiased <= xex) || (  spec.exMax.U  < xint) || io.x.bits.inf)
  val zzero = (io.x.bits.sgn === 1.U) &&
              ((xExUdfLimBiased <= xex) || ((-spec.exMin).U < xint) || io.x.bits.inf)

  zout.zex         := Mux(znan || zinf, Fill(exW, 1.U(1.W)),
                      Mux(zzero,        0.U(exW.W),
                                        zex0))
  zout.znan        := znan
  zout.zinf        := zinf
  zout.zzero       := zzero
  zout.zIsNonTable := znan || zinf || zzero

  io.z := Pipe(io.x.valid && io.info.valid, zout, stage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class ExpPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW       = spec.exW
  val manW      = spec.manW
  val padding   = bitW.outBitW - manW
  val extraBits = padding

  val zCorrTmpW = scala.math.min(6, 1 + manW + padding)

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info = new ExpPostProcInfo(spec, padding)
      val zres = UInt(bitW.outBitW.W)
    }))
    val z = Valid(UInt(spec.W.W))
  })

  val info = io.in.bits.info
  val zres = io.in.bits.zres

  val zman0 = if (extraBits > 0) (zres >> extraBits) +& zres(extraBits - 1) else zres

  val zCorrection = Wire(UInt(1.W))
  if (padding != 0) {
    val zCorrTerm0 = Cat(1.U(1.W), zman0(manW - 1, 0)) * info.zCorrCoef.get
    val zCorrTermW = 1 + manW + zCorrTmpW
    val zCorrTerm  = zCorrTerm0(zCorrTermW - 1, zCorrTermW - zCorrTmpW)
    assert(zCorrTerm.getWidth == zCorrTmpW)
    zCorrection := (zCorrTerm(zCorrTmpW - 1).asBool ||
                    zCorrTerm(zCorrTmpW - 2).asBool).asUInt
  } else {
    zCorrection := 0.U(1.W)
  }

  val zmanCorrected  = zman0 +& zCorrection
  val zexCorrection  = zmanCorrected(manW) & !info.zIsNonTable

  val zEx = info.zex + zexCorrection
  val zMan = Mux(info.zIsNonTable, Cat(info.znan, 0.U((manW - 1).W)),
                 zmanCorrected(manW - 1, 0))
  val z = Cat(0.U(1.W), zEx, zMan)
  assert(z.getWidth == spec.W)

  io.z := Pipe(io.in.valid, z, stage)
}
