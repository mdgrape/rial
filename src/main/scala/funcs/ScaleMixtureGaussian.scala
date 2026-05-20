//% @file ScaleMixtureGaussian.scala
//
// d/dx log p(x) for a 2-Gaussian mixture (with stddev sgmA / sgmB).
//
// Pipeline:
//   pre  : xIn   .. xTab    -- compute table addr, useTable flag, special
//                              flags; consume preMul to derive xsgmA2 (= x/sgmA^2)
//   tab  : xTab  .. xCalc
//   poly : xCalc .. xPost
//   mulArgs (comb) at xPost -- build  z0 = polynomial*2^maxBitDigit + 1
//                              and drive PostProcMultiplier with lhs = z0,
//                              rhs = xsgmA2.manW1
//   postMul : xPost   .. xPost+postMulStage
//   post    : xPost+postMulStage .. xOut
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.funcs

import scala.math._

import spire.math.Real
import spire.implicits._

import chisel3._
import chisel3.util._

import rial.approx._
import rial.arith.{RealSpec, RealGeneric, FloatChiselUtil, RoundSpec}
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._

// -------------------------------------------------------------------------
//   shared definitions
// -------------------------------------------------------------------------

object ScaleMixtureGaussianPreMulArgs {
  /** lhs width: just (1/sgmA^2).manW1, no extra precision pad. */
  def lhsW(spec: RealSpec): Int = { 1 + spec.manW }

  /** The (1/sgmA^2).manW1 constant as a UInt at width `w` (>= lhsW(spec)).
   *  When the shared multiplier is sized for a wider lhs (e.g. sincos), we
   *  pad the LSB side with zeros so the leading mantissa bits stay correct. */
  def lhs(sgmA: Double, spec: RealSpec, w: Int): UInt = {
    val base  = lhsW(spec)
    val diffW = w - base
    assert(diffW >= 0, f"PreMul lhsW($w) must be >= smg lhsW(${base})")
    val coef = new RealGeneric(spec, 1.0 / (sgmA * sgmA))
    val v    = coef.manW1.toBigInt
    if (diffW == 0) { v.U(base.W) } else { Cat(v.U(base.W), 0.U(diffW.W)) }
  }

  /** Width of the slice of the multiplier output that SMG consumes (the
   *  full product since we pre-shift the lhs into the high bits when the
   *  multiplier was sized wider for other functions). */
  def prodW(spec: RealSpec): Int = { (1 + spec.manW) + lhsW(spec) }

  def prod(spec: RealSpec, out: UInt): UInt = {
    val w = out.getWidth
    out(w - 1, w - prodW(spec))
  }
}

class ScaleMixtureGaussianInfo(val spec: RealSpec) extends Bundle {
  val xsgn       = UInt(1.W)
  val xnan       = Bool()
  val xinf       = Bool()
  val xzero      = Bool()
  val useTable   = Bool()
  val xsgmA2Ex   = UInt(spec.exW.W)
  val xsgmA2Man  = UInt(spec.manW.W)
}

// -------------------------------------------------------------------------
//  pre-process
// -------------------------------------------------------------------------

class ScaleMixtureGaussianPreProcess(
  val sgmA:  Double,
  val sgmB:  Double,
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int, // = preStage - preMulStage; pre-mult product arrives at xIn+preMulStage
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = bitW.addrW
  val fracW  = bitW.outBitW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  val xTableMaxDigit = ScaleMixtureGaussianTable.tableDomainDigit(fracW, sgmA, sgmB)
  assert(xTableMaxDigit <= spec.exMax,
    s"smg table domain digit ($xTableMaxDigit) must fit in the FP exponent (${spec.exMax})")

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val prod = Flipped(Valid(UInt(ScaleMixtureGaussianPreMulArgs.prodW(spec).W)))
    val adr  = Valid(UInt(adrW.W))
    val dx   = if (order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
    val info = Valid(new ScaleMixtureGaussianInfo(spec))
  })

  val x = io.x.bits

  // ---- table addr (uses x only) ------------------------------------------
  val xTableShift = (xTableMaxDigit + exBias).U(exW.W) - x.ex
  val useTable    = (xTableMaxDigit + exBias).U(exW.W) > x.ex

  val xTableMan = Cat(1.U(1.W), x.man) >> xTableShift
  val adr       = xTableMan(manW - 1, dxW)
  io.adr := Pipe(io.x.valid, adr, stage)

  if (order != 0) {
    val dx = Cat(~xTableMan(manW - adrW - 1), xTableMan(manW - adrW - 2, 0))
    io.dx.get := Pipe(io.x.valid, dx, stage)
  }

  // ---- xsgmA2 = x / sgmA^2 from preMul product ----------------------------
  // Mirror the rounding done by the old ScaleMixtureGaussianOtherPath.
  val rsgmA2 = new RealGeneric(spec, 1.0 / (sgmA * sgmA))
  val prodW  = ScaleMixtureGaussianPreMulArgs.prodW(spec)
  val prod   = io.prod.bits

  val moreThan2 = prod(prodW - 1)
  val roundBits = manW
  val sticky    = prod(roundBits - 2, 0).orR | (moreThan2 & prod(roundBits - 1))
  val round     = Mux(moreThan2, prod(roundBits),     prod(roundBits - 1))
  val shifted   = Mux(moreThan2, prod(prodW - 2, manW + 1),
                                 prod(prodW - 3, manW))
  val lsb       = shifted(0)
  val roundInc  = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven, lsb, round, sticky)
  val rounded   = shifted +& roundInc
  val moreThan2AfterRound = rounded(manW)
  val exInc     = moreThan2 | moreThan2AfterRound

  val xsgmA2Man = rounded(manW - 1, 0)

  assert(rsgmA2.toDouble > 1.0, "sgmA < 1.0, so 1 / sgmA^2 > 1.0")

  val xex0 = x.ex +& (rsgmA2.ex - spec.exBias).U +& exInc
  val xsgmA2Ex = Mux(xex0 > (spec.exMax + exBias).U,
    (spec.exMax + exBias + 1).U, xex0)

  // ---- info -- pipeline x.ex by preMulStage to align with the prod input -
  val info = Wire(new ScaleMixtureGaussianInfo(spec))
  info.xsgn      := x.sgn
  info.xnan      := x.nan
  info.xinf      := x.inf
  info.xzero     := x.zero
  info.useTable  := useTable
  info.xsgmA2Ex  := xsgmA2Ex
  info.xsgmA2Man := xsgmA2Man
  io.info := Pipe(io.x.valid && io.prod.valid, info, stage)
}

// -------------------------------------------------------------------------
//   table coefficient ROM (single table)
// -------------------------------------------------------------------------

class ScaleMixtureGaussianTableCoeff(
  val spec:  RealSpec,
  val coeff: PolynomialCoefficients,
  val stage: Int,
) extends Module {

  val adrW   = coeff.bitW.addrW
  val coefBW = coeff.bitW.coefBits.reduce(_+_)

  val io = IO(new Bundle {
    val adr = Flipped(Valid(UInt(adrW.W)))
    val cs  = Valid(UInt(coefBW.W))
  })

  val packed = VecInit(PolynomialGenerator.packTable(coeff).map(c => c.U(coefBW.W)))
  io.cs := Pipe(io.adr.valid, packed(io.adr.bits), stage)
}

// -------------------------------------------------------------------------
//   post-mult args (combinational)
// -------------------------------------------------------------------------
//
// Computes z0 = polynomial * 2^maxBitDigit + 1, normalized into a (1+fracW)
// mantissa.  Drives the shared PostProcMultiplier:
//   lhs = z0 mantissa with hidden 1
//   rhs = xsgmA2 mantissa with hidden 1
// also returns z0Ex (no-bias) for the post-process.

class ScaleMixtureGaussianPostMulArgsOutput(val spec: RealSpec, val fracW: Int) extends Bundle {
  val lhs    = UInt((1 + fracW).W)
  val rhs    = UInt((1 + spec.manW).W)
  val z0ex   = UInt(spec.exW.W)
}

class ScaleMixtureGaussianPostMulArgs(
  val sgmA: Double,
  val sgmB: Double,
  val spec: RealSpec,
  val bitW: PolynomialBitWidth,
) extends Module {
  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias
  val fracW  = bitW.outBitW

  val zTableScaleDigit = ScaleMixtureGaussianTable.tableMaxBitDigit(sgmA, sgmB)

  val io = IO(new Bundle {
    val useTable = Input(Bool())
    val xsgmA2Man = Input(UInt(manW.W))
    val zres   = Input(UInt(fracW.W))
    val valid  = Input(Bool())
    val mul    = Valid(new ScaleMixtureGaussianPostMulArgsOutput(spec, fracW))
  })

  // table value, gated to 0 outside the table domain
  val zTable       = Mux(io.useTable, io.zres, 0.U(fracW.W))
  val zTableScaled = Cat(zTable, 0.U(zTableScaleDigit.W)) // (fracW + zTableScaleDigit) bits

  // z0 = zTableScaled + 1.0 (where 1.0 = 1 << fracW)
  val z0 = zTableScaled +& Cat(1.U(1.W), 0.U(fracW.W))
  val z0W = (fracW + zTableScaleDigit + 1).U(log2Up(fracW + zTableScaleDigit + 2).W) -
            PriorityEncoder(Reverse(z0))
  val z0Shift = z0W - (1 + fracW).U

  // round to (1 + fracW)
  val z0Rounded0 = Cat(z0, 0.U(1.W)) >> z0Shift // pad LSB so the half-bit is captured
  val z0Rounded  = z0Rounded0(z0Rounded0.getWidth - 1, 1) +& z0Rounded0(0)
  val z0MoreThan2AfterRound = z0Rounded.head(1)

  val z0man = Mux(z0MoreThan2AfterRound === 1.U, z0Rounded >> 1, z0Rounded)

  val z0ex0 = z0W +& (exBias - 1 - fracW).U(exW.W) +& z0MoreThan2AfterRound
  val z0ex  = Mux(z0ex0.head(2).orR,
                  (spec.exMax + 1 + exBias).U(exW.W),
                  z0ex0(exW - 1, 0))

  val out = Wire(new ScaleMixtureGaussianPostMulArgsOutput(spec, fracW))
  out.lhs  := z0man(fracW, 0)
  out.rhs  := Cat(1.U(1.W), io.xsgmA2Man)
  out.z0ex := z0ex
  io.mul.valid := io.valid
  io.mul.bits  := out
}

// -------------------------------------------------------------------------
//   post-process
// -------------------------------------------------------------------------

class ScaleMixtureGaussianPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val info = Flipped(Valid(new ScaleMixtureGaussianInfo(spec)))
    val z0ex = Flipped(Valid(UInt(exW.W)))
    val mul  = Flipped(Valid(new PostProcMultiplierOutput(manW)))
    val z    = Valid(UInt(spec.W.W))
  })

  val info = io.info.bits
  val z0ex = io.z0ex.bits
  val mul  = io.mul.bits

  val zsgn = ~info.xsgn

  // zExNoBias = xsgmA2.ex + z0Ex + exInc - exBias (the xsgmA2 ex and z0 ex are
  // both biased, and we want a single bias in the result).
  val zex00 = info.xsgmA2Ex +& z0ex +& mul.exInc
  val zex0  = Mux(zex00 < exBias.U, 0.U, zex00 - exBias.U)
  val zex   = Mux(zex0 > (spec.exMax + exBias).U,
                  (spec.exMax + 1 + exBias).U(exW.W),
                  zex0(exW - 1, 0))

  // mantissa: zero out the mantissa for under/overflow.
  val zman = Mux(zex === 0.U || zex === (spec.exMax + 1 + exBias).U, 0.U(manW.W), mul.zman)

  val zNormal = Cat(zsgn, zex, zman)

  // ---- specials ----
  //   x = NaN  -> NaN
  //   x = +/-Inf -> -/+Inf
  //   x = +/-0   -> +/-0  (smg's natural answer for 0 is 0, no sign flip)
  val zNaN  = Cat(0.U(1.W), Fill(exW, 1.U(1.W)), 1.U(1.W), 0.U((manW - 1).W))
  val zInf  = Cat(zsgn, Fill(exW, 1.U(1.W)), 0.U(manW.W))
  val zZero = 0.U(spec.W.W)

  val z = Mux(info.xnan,  zNaN,
          Mux(info.xinf,  zInf,
          Mux(info.xzero, zZero,
                          zNormal)))

  io.z := Pipe(io.info.valid && io.z0ex.valid && io.mul.valid, z, stage)
}
