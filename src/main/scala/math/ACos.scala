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

import rial.table._
import rial.arith.{RealSpec, RealGeneric}
import rial.util._
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
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = bitW.addrW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ACosPhase1PreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x   = Flipped(Valid(new DecomposedReal(spec)))
    // 1-bit Sqrt ex(0) selector + adrW; same layout as SqrtPreProcess.io.adr
    val adr = Valid(UInt((1 + adrW).W))
    val dx  = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    // y = 1-|x| as a DecomposedReal, fed into SqrtNonTable (Mux'd with xTab)
    val y       = Valid(new DecomposedReal(spec))
    val flags   = Valid(new ACosFlags)
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  val xex  = x.ex
  val xman = x.man

  val xLargerThan1Comb = xex >= exBias.U

  // ---- compute 1 - |x| ----
  val xShiftMax = (1 + manW + 2).U
  val xShift0   = exBias.U(exW.W) - xex
  val xShift    = Mux(xShiftMax < xShift0, xShiftMax, xShift0)

  val xShifted        = (Cat(1.U(1.W), xman, 0.U(2.W)) >> xShift)
  val xSubtractedComb = Mux(xLargerThan1Comb, 0.U, ~(xShifted(manW + 2 - 1, 0)) +& 1.U)

  // --------------------------------------------------------------------------
  // stage "sub": the difference 1 - |x|

  private val cutSub = new SharedPipe("sub", cutIn.valid, stages.at("sub"))

  val xSubtracted  = cutSub(xSubtractedComb)
  val xLargerThan1 = cutSub(xLargerThan1Comb)
  val xSub         = cutSub(x)

  val xNormalizeShiftComb = PriorityEncoder(Reverse(xSubtracted))
  val xNormalizedComb     = (xSubtracted << xNormalizeShiftComb)(1 + manW + 2 - 1, 0)

  // --------------------------------------------------------------------------
  // stage "norm": the normalized difference

  private val cutNorm = new SharedPipe("norm", cutSub.valid, stages.at("norm"))

  val xNormalized      = cutNorm(xNormalizedComb)
  val xNormalizeShift  = cutNorm(xNormalizeShiftComb)
  val xLargerThan1Norm = cutNorm(xLargerThan1)
  val xNorm            = cutNorm(xSub)

  assert(!cutNorm.valid || xNormalized(1 + manW + 2 - 1) === 1.U ||
         xNormalized === 0.U)

  val xRounded          = xNormalized(manW + 1, 2) +& xNormalized(1)
  val xRoundedMoreThan2 = xRounded(manW)
  val xConvertedEx0     = exBias.U(exW.W) - xNormalizeShift + xRoundedMoreThan2
  val xConvertedMan0    = xRounded(manW - 1, 0)

  // special: y ~ 1, i.e. x ~ 0
  val xSmallEnoughComb = xNorm.zero || xConvertedEx0 === exBias.U

  val xConvertedExComb  = Mux(xSmallEnoughComb, exBias.U(exW.W),
                          Mux(xLargerThan1Norm, 0.U, xConvertedEx0))
  val xConvertedManComb = Mux(xSmallEnoughComb, 0.U(manW.W),
                          Mux(xLargerThan1Norm, 0.U, xConvertedMan0))

  // --------------------------------------------------------------------------
  // stage "round": the rounded difference

  private val cutRound = new SharedPipe("round", cutNorm.valid, stages.at("round"))

  val xConvertedEx  = cutRound(xConvertedExComb)
  val xConvertedMan = cutRound(xConvertedManComb)
  val xSmallEnough  = cutRound(xSmallEnoughComb)
  val xLargerThan1Round = cutRound(xLargerThan1Norm)
  val xRound        = cutRound(xNorm)

  // --------------------------------------------------------------------------
  // stage "out"

  // ---- table address: same encoding as SqrtPreProcess.adr (ex(0) | top mantissa) ----
  val adr = Cat(xConvertedEx(0), TableAddress.adr(xConvertedMan, adrW))
  io.adr := Pipe(cutRound.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = TableAddress.dx(xConvertedMan, adrW, dxW)
    io.dx.get := Pipe(cutRound.valid, dx, stages.at("out"))
  }

  // ---- y = 1-|x| as a DecomposedReal ----
  val yBits = Wire(new DecomposedReal(spec))
  yBits.sgn  := 0.U
  yBits.ex   := xConvertedEx
  yBits.man  := xConvertedMan
  yBits.zero := xLargerThan1Round
  yBits.inf  := xRound.inf
  yBits.nan  := xRound.nan
  io.y := Pipe(cutRound.valid, yBits, stages.at("out"))

  // ---- special-value flags ----
  val isSpecial = xRound.inf || xRound.nan || xSmallEnough || xLargerThan1Round
  val flag = MuxCase(0.U(ACosSpecialValue.width.W), Seq(
    (xRound.inf || xRound.nan)                -> ACosSpecialValue.xNaNInf,
     xSmallEnough                             -> ACosSpecialValue.xZero,
    (xLargerThan1Round && xRound.sgn === 1.U) -> ACosSpecialValue.xLargeNeg,
    (xLargerThan1Round && xRound.sgn === 0.U) -> ACosSpecialValue.xLargePos,
  ))

  val flagsBits = Wire(new ACosFlags)
  flagsBits.isSpecial := isSpecial
  flagsBits.xsgn      := xRound.sgn
  flagsBits.flag      := flag
  io.flags := Pipe(cutRound.valid, flagsBits, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ACosPhase1PreProcess]]. */
object ACosPhase1PreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"sub"`   the difference `1-|x|`
   *  - `"norm"`  the normalized difference
   *  - `"round"` the rounded difference
   *  - `"out"`   the table address, `1-|x|` and the flags
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("sub",   3.0),
    StagePoint("norm",  3.0),
    StagePoint("round", 2.0),
    StagePoint("out",   1.0))
}

class ACosPhase1PostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW       = spec.exW
  val manW      = spec.manW
  val extraBits = bitW.outBitW - manW

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ACosPhase1PostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info  = new SqrtPostProcInfo(spec) // from SqrtNonTable (run on y = 1-|x|)
      val zres  = UInt(bitW.outBitW.W)        // polynomial result
      val flags = new ACosFlags
    }))
    val z = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val in   = cutIn(io.in.bits)
  val zres = in.zres

  // Normal sqrt result (mirrors SqrtPostProcess)
  val zmanSqrtComb = WireDefault(0.U(manW.W))
  if (extraBits == 0) {
    zmanSqrtComb := zres
  } else {
    val zman0 = (zres >> extraBits) +& zres(extraBits - 1)
    val polynomialOvf = zman0(manW) === 1.U
    zmanSqrtComb := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW - 1, 0))
  }

  // --------------------------------------------------------------------------
  // stage "round": the rounded polynomial result

  private val cutRound = new SharedPipe("round", cutIn.valid, stages.at("round"))

  val zmanSqrt = cutRound(zmanSqrtComb)
  val info     = cutRound(in.info)
  val flags    = cutRound(in.flags)

  val zexSqrt = info.zex

  // --------------------------------------------------------------------------
  // stage "out"

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

  io.z := Pipe(cutRound.valid, Cat(zsgn, zex, zman), stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ACosPhase1PostProcess]]. */
object ACosPhase1PostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"round"` the rounded polynomial result
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("round", 2.0),
    StagePoint("out",   1.0))
}

// =========================================================================
//   Phase 2  (ACosPhase2Table is defined in ACosSim.scala)
// =========================================================================

class ACosPhase2PreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = bitW.addrW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ACosPhase2PreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x     = Flipped(Valid(new DecomposedReal(spec)))
    val adr   = Valid(UInt(adrW.W))
    val dx    = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    val flags = Valid(new ACosFlags)
    // x's mantissa (with hidden 1) -- the post-mult rhs at the right pipe stage
    val xManW1 = Valid(UInt((1 + manW).W))
    val xex    = Valid(UInt(exW.W))
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val xIn = cutIn(io.x.bits)

  val xShift       = exBias.U - xIn.ex
  val xAlignedComb = Cat(1.U(1.W), xIn.man) >> xShift

  // --------------------------------------------------------------------------
  // stage "align": the argument aligned at the binary point

  private val cutAlign = new SharedPipe("align", cutIn.valid, stages.at("align"))

  val xAligned = cutAlign(xAlignedComb)
  val x        = cutAlign(xIn)

  val xex  = x.ex
  val xman = x.man

  // --------------------------------------------------------------------------
  // stage "out"

  val adr = TableAddress.adr(xAligned(manW - 1, 0), adrW)
  io.adr := Pipe(cutAlign.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = TableAddress.dx(xAligned(manW - 1, 0), adrW, dxW)
    io.dx.get := Pipe(cutAlign.valid, dx, stages.at("out"))
  }

  // Decode phase-1's NaN-boxed flags (carried in x's exponent / mantissa).
  // Phase 1 emits xMan = 1/2/3 for the encoded special cases and a quiet NaN
  // mantissa (top bit only) for NaN/Inf -- any other value is treated as
  // xNaNInf.
  val xManMatches3 = (xman <= 3.U) && xman =/= 0.U
  val flags = Wire(new ACosFlags)
  flags.isSpecial := xex === Fill(exW, 1.U(1.W))
  flags.xsgn      := x.sgn
  flags.flag      := Mux(xManMatches3, xman(ACosSpecialValue.width - 1, 0),
                                       ACosSpecialValue.xNaNInf)
  io.flags := Pipe(cutAlign.valid, flags, stages.at("out"))

  // pass-through for the post-mult
  io.xManW1 := Pipe(cutAlign.valid, Cat(1.U(1.W), xman), stages.at("out"))
  io.xex    := Pipe(cutAlign.valid, xex,                 stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ACosPhase2PreProcess]]. */
object ACosPhase2PreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"align"` the argument aligned at the binary point
   *  - `"out"`   the table address and the flags
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("align", 3.0),
    StagePoint("out",   1.0))
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
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val realPi     = new RealGeneric(spec, Pi)
  val realHalfPi = new RealGeneric(spec, Pi * 0.5)

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ACosPhase2PostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val flags = Flipped(Valid(new ACosFlags))
    val xex   = Flipped(Valid(UInt(exW.W))) // x.ex pipelined from phase-2 pre
    val mul   = Flipped(Valid(new PostProcMultiplierOutput(manW)))
    val z     = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in",
    io.flags.valid && io.mul.valid && io.xex.valid, stages.at("in"))

  val mul     = cutIn(io.mul.bits)
  val xexIn   = cutIn(io.xex.bits)
  val flagsIn = cutIn(io.flags.bits)

  val zMan0Comb = mul.zman
  val zEx0Comb  = xexIn + mul.exInc // mul lhs has ex==exBias unbiased=0; total = x.ex + exInc

  // ---- pi - z (for x < 0) ----
  val piW = 2 + manW + 2
  val pi  = (realPi.manW1 << 3).toBigInt.U(piW.W)

  val zAlignedComb = Cat(1.U(2.W), zMan0Comb, 0.U(2.W)) >> (exBias.U - zEx0Comb)

  // --------------------------------------------------------------------------
  // stage "align": the product aligned with pi

  private val cutAlign = new SharedPipe("align", cutIn.valid, stages.at("align"))

  val zAligned   = cutAlign(zAlignedComb)
  val zMan0      = cutAlign(zMan0Comb)
  val zEx0       = cutAlign(zEx0Comb)
  val flagsAlign = cutAlign(flagsIn)

  val zNegComb = pi - zAligned

  // --------------------------------------------------------------------------
  // stage "neg": the difference pi - z

  private val cutNeg = new SharedPipe("neg", cutAlign.valid, stages.at("neg"))

  val zNeg       = cutNeg(zNegComb)
  val zMan0Neg   = cutNeg(zMan0)
  val zEx0Neg    = cutNeg(zEx0)
  val flagsNeg   = cutNeg(flagsAlign)

  val zNegMoreThan2 = zNeg(piW - 1)
  val zNegShifted   = Mux(zNegMoreThan2 === 1.U, zNeg(piW - 1, 3), zNeg(piW - 2, 2))
  val zNegRoundedComb = zNegShifted + Mux(zNegMoreThan2 === 1.U, zNeg(2), zNeg(1))
  val zNegEx0Comb     = exBias.U(exW.W) + zNegMoreThan2

  // --------------------------------------------------------------------------
  // stage "round": the rounded difference

  private val cutRound = new SharedPipe("round", cutNeg.valid, stages.at("round"))

  val zNegRounded = cutRound(zNegRoundedComb)
  val zNegEx0     = cutRound(zNegEx0Comb)
  val zMan0Round  = cutRound(zMan0Neg)
  val zEx0Round   = cutRound(zEx0Neg)
  val flags       = cutRound(flagsNeg)

  val zNegMan0 = zNegRounded(manW - 1, 0)

  // --------------------------------------------------------------------------
  // stage "out"

  // ---- decode special values from flags ----
  val xNaNInf   = flags.isSpecial && flags.flag === ACosSpecialValue.xNaNInf
  val xZero     = flags.isSpecial && flags.flag === ACosSpecialValue.xZero
  val xLargePos = flags.isSpecial && flags.flag === ACosSpecialValue.xLargePos
  val xLargeNeg = flags.isSpecial && flags.flag === ACosSpecialValue.xLargeNeg

  val zman = Mux(xNaNInf,   Cat(1.U(1.W), 0.U((manW - 1).W)),
             Mux(xZero,     realHalfPi.man.toBigInt.U(manW.W),
             Mux(xLargePos, 0.U(manW.W),
             Mux(xLargeNeg, realPi.man.toBigInt.U(manW.W),
             Mux(flags.xsgn === 0.U, zMan0Round, zNegMan0)))))

  val zex = Mux(xNaNInf,   Fill(exW, 1.U(1.W)),
            Mux(xZero,     realHalfPi.ex.U(exW.W),
            Mux(xLargePos, 0.U(exW.W),
            Mux(xLargeNeg, realPi.ex.U(exW.W),
            Mux(flags.xsgn === 0.U, zEx0Round, zNegEx0)))))

  val z = Cat(0.U(1.W), zex, zman)
  io.z := Pipe(cutRound.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ACosPhase2PostProcess]]. */
object ACosPhase2PostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"align"` the product aligned with pi
   *  - `"neg"`   the difference pi - z
   *  - `"round"` the rounded difference
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("align", 3.0),
    StagePoint("neg",   2.0),
    StagePoint("round", 2.0),
    StagePoint("out",   1.0))
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
