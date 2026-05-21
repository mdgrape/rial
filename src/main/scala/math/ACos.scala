//% @file ACos.scala
//
// acos(x) hardware modules. Two-phase computation:
//
//   phase 1:  x -> y = sqrt(1 - |x|)   -- reuses the Sqrt sub-table & non-table
//   phase 2:  y -> acos(1 - y^2)       -- own table; uses the PostProcMultiplier
//
// Phase-1 shares the SqrtTableCoeff (its `adr` is Mux'd with Sqrt's) and the
// SqrtNonTable (its input is Mux'd with xTab). Phase-2 has its own table.
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math._

import spire.math.Real
import spire.implicits._

import chisel3._
import chisel3.util._

import rial.approx._
import rial.arith.{RealSpec, RealGeneric}
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._

// -------------------------------------------------------------------------
//   shared definitions
// -------------------------------------------------------------------------

object ACosSpecialValue {
  val width     = 2
  val xNaNInf   = 0.U(width.W) // x == nan/inf -> nan
  val xZero     = 1.U(width.W) // x ==  0      -> pi/2
  val xLargePos = 2.U(width.W) // x >=  1      -> +0
  val xLargeNeg = 3.U(width.W) // x <= -1      -> pi
}

class ACosFlags extends Bundle {
  val isSpecial = Bool()
  val xsgn      = UInt(1.W)
  val flag      = UInt(ACosSpecialValue.width.W)
}

// =========================================================================
//   Phase 1
// =========================================================================

class ACosPhase1PreProcess(
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
    val x   = Flipped(Valid(new DecomposedReal(spec)))
    // 1-bit Sqrt ex(0) selector + adrW; same layout as SqrtPreProcess.io.adr
    val adr = Valid(UInt((1 + adrW).W))
    val dx  = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    // y = 1-|x| as a DecomposedReal, fed into SqrtNonTable (Mux'd with xTab)
    val y       = Valid(new DecomposedReal(spec))
    val flags   = Valid(new ACosFlags)
  })

  val xsgn = io.x.bits.sgn
  val xex  = io.x.bits.ex
  val xman = io.x.bits.man

  val xLargerThan1 = xex >= exBias.U

  // ---- compute 1 - |x| ----
  val xShiftMax = (1 + manW + 2).U
  val xShift0   = exBias.U(exW.W) - xex
  val xShift    = Mux(xShiftMax < xShift0, xShiftMax, xShift0)

  val xShifted    = (Cat(1.U(1.W), xman, 0.U(2.W)) >> xShift)
  val xSubtracted = Mux(xLargerThan1, 0.U, ~(xShifted(manW + 2 - 1, 0)) +& 1.U)

  val xNormalizeShift = PriorityEncoder(Reverse(xSubtracted))
  val xNormalized     = (xSubtracted << xNormalizeShift)(1 + manW + 2 - 1, 0)
  assert(xNormalized(1 + manW + 2 - 1) === 1.U || xNormalized === 0.U)

  val xRounded          = xNormalized(manW + 1, 2) +& xNormalized(1)
  val xRoundedMoreThan2 = xRounded(manW)
  val xConvertedEx0     = exBias.U(exW.W) - xNormalizeShift + xRoundedMoreThan2
  val xConvertedMan0    = xRounded(manW - 1, 0)

  val xSmallEnough = io.x.bits.zero || xConvertedEx0 === exBias.U // special: y ~ 1, i.e. x ~ 0

  val xConvertedEx  = Mux(xSmallEnough, exBias.U(exW.W), Mux(xLargerThan1, 0.U, xConvertedEx0))
  val xConvertedMan = Mux(xSmallEnough, 0.U(manW.W),     Mux(xLargerThan1, 0.U, xConvertedMan0))

  // ---- table address: same encoding as SqrtPreProcess.adr (ex(0) | top mantissa) ----
  val adr = Cat(xConvertedEx(0), xConvertedMan(manW - 1, dxW))
  io.adr := Pipe(io.x.valid, adr, stage)

  if (order != 0) {
    val dx = Cat(~xConvertedMan(dxW - 1), xConvertedMan(dxW - 2, 0))
    io.dx.get := Pipe(io.x.valid, dx, stage)
  }

  // ---- y = 1-|x| as a DecomposedReal ----
  val yBits = Wire(new DecomposedReal(spec))
  yBits.sgn  := 0.U
  yBits.ex   := xConvertedEx
  yBits.man  := xConvertedMan
  yBits.zero := xLargerThan1
  yBits.inf  := io.x.bits.inf
  yBits.nan  := io.x.bits.nan
  io.y := Pipe(io.x.valid, yBits, stage)

  // ---- special-value flags ----
  val isSpecial = io.x.bits.inf || io.x.bits.nan || xSmallEnough || xLargerThan1
  val flag = MuxCase(0.U(ACosSpecialValue.width.W), Seq(
    (io.x.bits.inf || io.x.bits.nan) -> ACosSpecialValue.xNaNInf,
     xSmallEnough                    -> ACosSpecialValue.xZero,
    (xLargerThan1 && xsgn === 1.U)   -> ACosSpecialValue.xLargeNeg,
    (xLargerThan1 && xsgn === 0.U)   -> ACosSpecialValue.xLargePos,
  ))

  val flagsBits = Wire(new ACosFlags)
  flagsBits.isSpecial := isSpecial
  flagsBits.xsgn      := xsgn
  flagsBits.flag      := flag
  io.flags := Pipe(io.x.valid, flagsBits, stage)
}

class ACosPhase1PostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW       = spec.exW
  val manW      = spec.manW
  val extraBits = bitW.outBitW - manW

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info  = new SqrtPostProcInfo(spec) // from SqrtNonTable (run on y = 1-|x|)
      val zres  = UInt(bitW.outBitW.W)        // polynomial result
      val flags = new ACosFlags
    }))
    val z = Valid(UInt(spec.W.W))
  })

  val info  = io.in.bits.info
  val zres  = io.in.bits.zres
  val flags = io.in.bits.flags

  // Normal sqrt result (mirrors SqrtPostProcess)
  val zexSqrt  = info.zex
  val zmanSqrt = WireDefault(0.U(manW.W))
  if (extraBits == 0) {
    zmanSqrt := zres
  } else {
    val zman0 = (zres >> extraBits) +& zres(extraBits - 1)
    val polynomialOvf = zman0(manW) === 1.U
    zmanSqrt := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW - 1, 0))
  }

  // Special-value encoding (NaN-boxed for phase 2 to decode).
  // For an ordinary NaN/Inf input we emit a quiet NaN with mantissa MSB only
  // set (matches RealGeneric.nan); for the other cases we encode the flag
  // value in the mantissa.
  val zsgnSpecial = flags.xsgn
  val zexSpecial  = Fill(exW, 1.U(1.W))
  val zmanSpecial = Mux(flags.flag === ACosSpecialValue.xNaNInf,
                        Cat(1.U(1.W), 0.U((manW - 1).W)),
                        flags.flag)

  val zsgn = zsgnSpecial // sqrt is always positive; sign carries x's
  val zex  = Mux(flags.isSpecial, zexSpecial,  zexSqrt)
  val zman = Mux(flags.isSpecial, zmanSpecial, zmanSqrt)

  io.z := Pipe(io.in.valid, Cat(zsgn, zex, zman), stage)
}

// =========================================================================
//   Phase 2  (ACosPhase2Table is defined in ACosSim.scala)
// =========================================================================

class ACosPhase2PreProcess(
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
    val x     = Flipped(Valid(new DecomposedReal(spec)))
    val adr   = Valid(UInt(adrW.W))
    val dx    = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    val flags = Valid(new ACosFlags)
    // x's mantissa (with hidden 1) -- the post-mult rhs at the right pipe stage
    val xManW1 = Valid(UInt((1 + manW).W))
    val xex    = Valid(UInt(exW.W))
  })

  val xex  = io.x.bits.ex
  val xman = io.x.bits.man

  val xShift   = exBias.U - xex
  val xAligned = Cat(1.U(1.W), xman) >> xShift

  val adr = xAligned(manW - 1, dxW)
  io.adr := Pipe(io.x.valid, adr, stage)

  if (order != 0) {
    val dx = Cat(~xAligned(dxW - 1), xAligned(dxW - 2, 0))
    io.dx.get := Pipe(io.x.valid, dx, stage)
  }

  // Decode phase-1's NaN-boxed flags (carried in x's exponent / mantissa).
  // Phase 1 emits xMan = 1/2/3 for the encoded special cases and a quiet NaN
  // mantissa (top bit only) for NaN/Inf -- any other value is treated as
  // xNaNInf.
  val xManMatches3 = (xman <= 3.U) && xman =/= 0.U
  val flags = Wire(new ACosFlags)
  flags.isSpecial := xex === Fill(exW, 1.U(1.W))
  flags.xsgn      := io.x.bits.sgn
  flags.flag      := Mux(xManMatches3, xman(ACosSpecialValue.width - 1, 0),
                                       ACosSpecialValue.xNaNInf)
  io.flags := Pipe(io.x.valid, flags, stage)

  // pass-through for the post-mult
  io.xManW1 := Pipe(io.x.valid, Cat(1.U(1.W), xman), stage)
  io.xex    := Pipe(io.x.valid, xex,                  stage)
}

class ACosPhase2TableCoeff(
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

class ACosPhase2PostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val realPi     = new RealGeneric(spec, Pi)
  val realHalfPi = new RealGeneric(spec, Pi * 0.5)

  val io = IO(new Bundle {
    val flags = Flipped(Valid(new ACosFlags))
    val xex   = Flipped(Valid(UInt(exW.W))) // x.ex pipelined from phase-2 pre
    val mul   = Flipped(Valid(new PostProcMultiplierOutput(manW)))
    val z     = Valid(UInt(spec.W.W))
  })

  val zMan0 = io.mul.bits.zman
  val zEx0  = io.xex.bits + io.mul.bits.exInc // mul lhs has ex==exBias unbiased=0; total = x.ex + exInc

  // ---- pi - z (for x < 0) ----
  val piW = 2 + manW + 2
  val pi  = (realPi.manW1 << 3).toBigInt.U(piW.W)

  val zAligned     = Cat(1.U(2.W), zMan0, 0.U(2.W)) >> (exBias.U - zEx0)
  val zNeg         = pi - zAligned
  val zNegMoreThan2 = zNeg(piW - 1)
  val zNegShifted   = Mux(zNegMoreThan2 === 1.U, zNeg(piW - 1, 3), zNeg(piW - 2, 2))
  val zNegRounded   = zNegShifted + Mux(zNegMoreThan2 === 1.U, zNeg(2), zNeg(1))
  val zNegEx0       = exBias.U(exW.W) + zNegMoreThan2
  val zNegMan0      = zNegRounded(manW - 1, 0)

  // ---- decode special values from flags ----
  val flags = io.flags.bits
  val xNaNInf   = flags.isSpecial && flags.flag === ACosSpecialValue.xNaNInf
  val xZero     = flags.isSpecial && flags.flag === ACosSpecialValue.xZero
  val xLargePos = flags.isSpecial && flags.flag === ACosSpecialValue.xLargePos
  val xLargeNeg = flags.isSpecial && flags.flag === ACosSpecialValue.xLargeNeg

  val zman = Mux(xNaNInf,   Cat(1.U(1.W), 0.U((manW - 1).W)),
             Mux(xZero,     realHalfPi.man.toBigInt.U(manW.W),
             Mux(xLargePos, 0.U(manW.W),
             Mux(xLargeNeg, realPi.man.toBigInt.U(manW.W),
             Mux(flags.xsgn === 0.U, zMan0, zNegMan0)))))

  val zex = Mux(xNaNInf,   Fill(exW, 1.U(1.W)),
            Mux(xZero,     realHalfPi.ex.U(exW.W),
            Mux(xLargePos, 0.U(exW.W),
            Mux(xLargeNeg, realPi.ex.U(exW.W),
            Mux(flags.xsgn === 0.U, zEx0, zNegEx0)))))

  val z = Cat(0.U(1.W), zex, zman)
  io.z := Pipe(io.flags.valid && io.mul.valid && io.xex.valid, z, stage)
}

// Combinational helper: drives the PostProcMultiplier inputs (lhs = polynomial
// with hidden 1, rhs = x.manW1 = original mantissa with hidden 1).
class ACosPhase2PostMultArgsOutput(val spec: RealSpec, val fracW: Int) extends Bundle {
  val lhs = UInt((1 + fracW).W)
  val rhs = UInt((1 + spec.manW).W)
}

class ACosPhase2PostMultArgs(
  val spec: RealSpec,
  val bitW: PolynomialBitWidth,
) extends Module {
  val fracW = bitW.outBitW

  val io = IO(new Bundle {
    val zres   = Input(UInt(fracW.W))
    val xManW1 = Flipped(Valid(UInt((1 + spec.manW).W)))
    val mul    = Valid(new ACosPhase2PostMultArgsOutput(spec, fracW))
  })

  val out = Wire(new ACosPhase2PostMultArgsOutput(spec, fracW))
  out.lhs := Cat(1.U(1.W), io.zres)                // polynomial value with hidden 1
  out.rhs := io.xManW1.bits                         // x's mantissa with hidden 1

  io.mul.valid := io.xManW1.valid
  io.mul.bits  := out
}
