//% @file SinCos.scala
//
// sin(x) and cos(x): cos is computed as sin(π/2 - x), so the two share the
// entire datapath except for one sign/fold conditional in the preprocess.
//
// Argument reduction (in the preprocess):
//   1. Multiply |x| by 1/π via the shared PreProcMultiplier; the multiplier
//      output `prod` carries x/π at full precision (no rounding) so we can
//      subtract integer multiples of 1/2 without losing low bits.
//   2. Shift the quotient into [0, 2), then fold via the 2 MSBs:
//         00 -> sin(y) [direct]   01 -> sin(1-y)
//         10 -> -sin(1-y)         11 -> -sin(y)
//      Cos uses (π/2 - x); equivalently a different fold per quadrant.
//   3. The folded result is y ∈ [0, 1/2); store it as yex/yman.
//
// The polynomial table approximates z(2y) = sin(π·y)/y in [1, π/2).
// The post-mult then computes z * y = sin(πy) = ±sin(x).
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import scala.math.{max => smax}

import spire.math.Real
import spire.implicits._

import chisel3._
import chisel3.util._

import rial.approx._
import rial.arith.{RealSpec, FloatChiselUtil}
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._

// -------------------------------------------------------------------------
//
// constants for the pre-process multiplier (`1/π`)
//

object SinCosPreMulArgs {

  /** Minimum lhs width that the pre-mult must have for sin/cos. The trailing
   *  pad matches the SW reference's default `oneOverPiPad = max(manW+1, 12)`
   *  so software and hardware produce bit-identical x/π.
   */
  def lhsW(spec: RealSpec): Int = {
    1 + spec.manW + smax(spec.manW + 1, 12)
  }

  /** The 1/π constant as a UInt at width `w` (≥ lhsW(spec)). When the
   *  multiplier was sized for a different function (e.g. exp) we pad the
   *  LSB side with zeros so the high bits still carry the same fraction. */
  def lhs(spec: RealSpec, w: Int): UInt = {
    val base  = lhsW(spec)
    val diffW = w - base
    assert(diffW >= 0, f"PreMul lhsW($w) must be >= sincos lhsW(${base})")
    val coef = Real.one / Real.pi
    val v    = coef(base + 1).toBigInt
    if (diffW == 0) v.U(base.W) else Cat(v.U(base.W), 0.U(diffW.W))
  }

  /** Width of the slice of the multiplier output that sincos actually
   *  consumes (the high bits). */
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

/** Information emitted by [[SinCosPreProcess]] for [[SinCosPostProcess]]
 *  (everything that doesn't go through the polynomial path). */
class SinCosPreProcInfo(val spec: RealSpec) extends Bundle {
  val znan = Bool()
  val zone = Bool()
  val zsgn = UInt(1.W)
  val yex  = UInt(spec.exW.W)
  val yman = UInt(spec.manW.W)
}

class SinCosPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = bitW.addrW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  val io = IO(new Bundle {
    val isSin = Input(Bool())
    val x     = Flipped(Valid(new DecomposedReal(spec)))
    val prod  = Flipped(Valid(UInt(SinCosPreMulArgs.prodW(spec).W)))
    val adr   = Valid(UInt(adrW.W))
    val dx    = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    val info  = Valid(new SinCosPreProcInfo(spec))
  })

  val xsgn = io.x.bits.sgn
  val xex  = io.x.bits.ex
  val xman = io.x.bits.man

  // --------------------------------------------------------------------------
  // x/π from the pre-mult output

  val oneOverPiW = SinCosPreMulArgs.lhsW(spec)

  val xOverPiProd          = Mux(xex === 0.U, 0.U, io.prod.bits)
  val xOverPiProdMoreThan2 = xOverPiProd((1 + manW) + oneOverPiW - 1)
  val xOverPiEx            = xex -& 2.U + xOverPiProdMoreThan2
  val xOverPiIsZero        = xOverPiEx(exW)
  // x.ex > xOverPiEx always (xOverPiProdMoreThan2 ≤ 1); overflow here just
  // means x/π is tiny, not infinite.

  val xOverPi      = Mux(xOverPiProdMoreThan2, xOverPiProd,
                                                Cat(xOverPiProd.tail(1), 0.U(1.W)))
  val xOverPiFracW = (1 + manW) + oneOverPiW - 1
  assert(xOverPi.getWidth == xOverPiFracW + 1)

  // --------------------------------------------------------------------------
  // fold x/π into [0, 2)

  val xOverPiExMoreThan1 = (xOverPiEx >= exBias.U)
  val xOverPiExNobiasPos = (xOverPiEx - exBias.U)(log2Up(xOverPiFracW) - 1, 0)
  val xOverPiExNobiasNeg = (exBias.U - xOverPiEx)(log2Up(xOverPiFracW) - 1, 0)

  val xOverPiAlignPos = (xOverPi << xOverPiExNobiasPos)(xOverPiFracW, 0)
  val xOverPiAlignNeg = (xOverPi >> xOverPiExNobiasNeg)(xOverPiFracW, 0)
  val xOverPiAligned  = Mux(xOverPiExMoreThan1, xOverPiAlignPos, xOverPiAlignNeg)

  // --------------------------------------------------------------------------
  // sign of result, derived from the 2 MSBs of (|x|/π mod 2)
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

  val xOverPiAligned2MSBs = xOverPiAligned.head(2)

  val zsgn = Mux(io.isSin, xOverPiAligned2MSBs(1) ^ xsgn,
                           xOverPiAligned2MSBs(1) ^ xOverPiAligned2MSBs(0))

  val znan = xex.andR

  // --------------------------------------------------------------------------
  // fold sin/cos in [0, 2) into sin [0, 1/2)
  //
  // sin(x)                 cos(x)              MSB(0) | isSin  |  used
  // 2MSBs 00 01 10 11      2MSBs 00 01 10 11   -------+--------+--------
  //      |  :  :  :  :          |  :  :  :  :     0   | sin(1) | ymanPos
  //      | .-. :  :  :        .-|-.:  :  :.-:     1   | sin(1) | ymanNeg
  // _____|'_:_': _:__:     __.__|__:__:__:__:     0   | cos(0) | ymanNeg
  //   _.'|  :  :'._.':          |  :'.:.':  :     1   | cos(0) | ymanPos
  //      |  :  :  :  :          |  :  :  :  :
  //         '--'  '--'          '--'  '--'
  //         pi-x                pi/2-x

  val ymanPos =  xOverPiAligned.tail(2)
  val ymanNeg = ~xOverPiAligned.tail(2)
  val ymanAligned = Mux(xOverPiAligned2MSBs(0) ^ io.isSin, ymanPos, ymanNeg)

  // For 0 < x < π/2 we use x/π itself (un-aligned) to keep precision.
  val ymanAsIs = xOverPi(xOverPi.getWidth - 1, 2)
  val yexAsIs  = xOverPiEx
  val nonAlign = !xOverPiExMoreThan1 && xOverPiAligned2MSBs === 0.U && io.isSin

  assert(ymanAsIs.getWidth == ymanPos.getWidth)
  val yman0         = Mux(nonAlign, ymanAsIs, ymanAligned)
  val yex0          = Mux(nonAlign, yexAsIs,  (exBias - 2).U(exW.W))
  val ymanIsNonZero = yman0.orR

  // round/normalize yman0 to a manW-bit mantissa, like the old code.
  val yman0W        = yman0.getWidth
  val yman0Shift0   = PriorityEncoder(Reverse(yman0))
  val yman0ShiftW   = log2Up(yman0W)
  val yman0Shift    = yman0Shift0(yman0ShiftW - 1, 0)
  val yman0Shifted  = (yman0 << yman0Shift)(yman0W - 1, 0)
  val yman0RoundBit = yman0W - manW - 1
  val yman0Rounded  = yman0Shifted(yman0W - 2, yman0RoundBit) +&
                      yman0Shifted(yman0RoundBit - 1)
  val yman0MoreThan2 = yman0Rounded(manW)

  val ymanNonZero = yman0Rounded(manW - 1, 0)
  val yexNonZero  = yex0 - yman0Shift + yman0MoreThan2
  val yex  = Mux(!ymanIsNonZero, 0.U, yexNonZero)
  val yman = Mux(!ymanIsNonZero, 0.U, ymanNonZero)

  val zone = (yex === (exBias - 1).U) && yman === 0.U

  // --------------------------------------------------------------------------
  // address / dx for the polynomial

  val yaligned = Cat(1.U(1.W), yman) >> ((exBias - 1).U(exW.W) - yex)
  val adr      = yaligned(manW - 1, manW - adrW)

  // --------------------------------------------------------------------------
  // outputs

  val infoOut = Wire(new SinCosPreProcInfo(spec))
  infoOut.znan := znan
  infoOut.zone := zone
  infoOut.zsgn := zsgn
  infoOut.yex  := yex
  infoOut.yman := yman

  // valid is the AND of x.valid and prod.valid (we need both)
  val valid = io.x.valid && io.prod.valid

  io.info := Pipe(valid, infoOut, stage)
  io.adr  := Pipe(valid, adr,     stage)
  if (order != 0) {
    val dx = Cat(~yaligned(manW - adrW - 1), yaligned(manW - adrW - 2, 0))
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

class SinCosTableCoeff(
  val spec:  RealSpec,
  val coeff: PolynomialCoefficients,
  val stage: Int,
) extends Module {

  val adrW  = coeff.bitW.addrW
  val order = coeff.bitW.nOrder

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
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

/** Combines [[SinCosPreProcess]]'s info and the post-mult output into the
 *  final sin(x) / cos(x) bit pattern. */
class SinCosPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val pre   = Flipped(Valid(new SinCosPreProcInfo(spec)))
    val mul   = Flipped(Valid(new PostProcMultiplierOutput(manW)))
    val z     = Valid(UInt(spec.W.W))
  })

  val info  = io.pre.bits
  val zman0 = io.mul.bits.zman
  val zex0  = info.yex + 1.U + io.mul.bits.exInc

  val znan  = info.znan
  val zone  = info.zone
  val zzero = info.yex === 0.U
  val zsgn  = info.zsgn

  val zex = Mux(zzero, 0.U,
            Mux(znan,  Fill(exW, 1.U(1.W)),
            Mux(zone,  exBias.U(exW.W),
                       zex0)))
  val zman = Mux(zzero, 0.U,
             Mux(znan || zone, Cat(znan, 0.U((manW - 1).W)),
                               zman0))

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := Pipe(io.pre.valid && io.mul.valid, z, stage)
}
