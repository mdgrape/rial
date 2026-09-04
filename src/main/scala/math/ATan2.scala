//% @file ATan2.scala
//
// atan2(y, x) hardware modules. Two-phase computation:
//
//   phase 1: (y, x) -> t = min(|x|,|y|) / max(|x|,|y|)
//   phase 2:   t    -> atan2(y, x)
//
// Phase 1 shares the ReciprocalTableCoeff
// Both phases need the PostProcMultiplier:
//   Phase 1 multiplies 2/(1+max) * (1+min.man) to get the ratio (1+t.man).
//   Phase 2 multiplies atan(t)/t * t to get atan(t).
//
// Phase 1's special-value flags & ratio's exponent/status are NaN-boxed into the
// output so Phase 2 can decode them with the standard sgn/ex/man triplet.
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

object ATan2Status {
  val width = 2
  // status encoding mirrors ATan2Phase1SimGeneric: status = (yIsLarger << 1) | xsgn
  val xIsPosIsLarger  = 0.U(width.W) // ysgn *        atan(|y/x|)   .. x>0, |x|>|y|
  val xIsNegIsLarger  = 1.U(width.W) // ysgn * (pi  - atan(|y/x|))  .. x<0, |x|>|y|
  val xIsPosIsSmaller = 2.U(width.W) // ysgn * (pi/2- atan(|x/y|))  .. x>0, |x|<|y|
  val xIsNegIsSmaller = 3.U(width.W) // ysgn * (pi/2+ atan(|x/y|))  .. x<0, |x|<|y|
}

// Mantissa values written by phase 1 when isSpecial is set; phase 2 decodes by
// exact match.  Anything else (e.g. the NaN top-bit-only encoding) is treated
// as zNaN.
object ATan2SpecialValue {
  val width        = 3
  val zNormal      = 0.U(width.W) // normal atan2
  val zZero        = 1.U(width.W) // |x| >> |y|, x > 0 : atan2(y,x) = (+/-)0
  val zPi          = 2.U(width.W) // |x| >> |y|, x < 0 : atan2(y,x) = (+/-)pi
  val zHalfPi      = 3.U(width.W) // |y| >> |x|        : atan2(y,x) = (+/-)pi/2
  val zQuarterPi   = 4.U(width.W) // |x| == |y|, x > 0 : atan2(y,x) = (+/-)pi/4
  val z3QuarterPi  = 5.U(width.W) // |x| == |y|, x < 0 : atan2(y,x) = (+/-)3pi/4
  val zNaN         = 7.U(width.W) // x or y is NaN
}

class ATan2Flags extends Bundle {
  val isSpecial = Bool()
  val special   = UInt(ATan2SpecialValue.width.W)
  val status    = UInt(ATan2Status.width.W)
  val ysgn      = UInt(1.W)
}

// =========================================================================
//   Phase 1
// =========================================================================

class ATan2Phase1Info(val spec: RealSpec) extends Bundle {
  val zex       = UInt(spec.exW.W) // exponent of the |min| / |max| ratio
  val xySameMan = Bool()
  val maxXYMan0 = Bool()
  val minMan    = UInt(spec.manW.W)
}

class ATan2Phase1PreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  // We pack status1 into ex(exW-1); the sign of NaN encodes ysgn.  This relies
  // on the bias having its MSB clear (e.g. 0b0111_1111) so that "ex(exW-1) == 0"
  // is true for any |ratio| < 1 -- if it weren't, ratio's natural MSB would
  // collide with the status bit.
  assert(spec.exBias == maskI(spec.exW - 1), """
    |ATan2Phase1 calculates |min(x,y)| / |max(x,y)| and encodes information
    |like xsgn, ysgn, |x| >? |y| into the sign, msb of exponent, lsb of mantissa.
    |It depends on |min(x,y)| / |max(x,y)| <= 1 and in case of == 1 we use
    |special value flags, so we don't use the msb of the exponent.
    |It means that, if x:float < 1, x.ex.msb must be zero.""".stripMargin)

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = bitW.addrW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ATan2Phase1PreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x = Flipped(Valid(new DecomposedReal(spec)))
    val y = Flipped(Valid(new DecomposedReal(spec)))

    // shared ReciprocalTableCoeff address / interpolation residual:
    val adr = Valid(UInt(adrW.W))
    val dx  = if (order != 0) { Some(Valid(UInt(dxW.W))) } else { None }

    val flags = Valid(new ATan2Flags)
    val info  = Valid(new ATan2Phase1Info(spec))
  })

  // both operands are needed
  private val cutIn = new SharedPipe("in", io.x.valid && io.y.valid,
    stages.at("in"))

  val xIn = cutIn(io.x.bits)
  val yIn = cutIn(io.y.bits)

  // ---- compare |x|, |y| (sign-free since both are non-negative bit-patterns
  // when interpreted as UInts assembled from ex ## man) ----
  val xExMan = Cat(xIn.ex, xIn.man)
  val yExMan = Cat(yIn.ex, yIn.man)
  val yIsLargerComb = xExMan < yExMan

  // --------------------------------------------------------------------------
  // stage "cmp": which of |x| and |y| is the larger

  private val cutCmp = new SharedPipe("cmp", cutIn.valid, stages.at("cmp"))

  val yIsLarger = cutCmp(yIsLargerComb)
  val x         = cutCmp(xIn)
  val y         = cutCmp(yIn)

  // ---- exponent of min/max ratio ----
  val minex     = Mux(yIsLarger, x.ex, y.ex)
  val maxex     = Mux(yIsLarger, y.ex, x.ex)
  val diffexDec = Mux(yIsLarger, y.man > x.man, x.man > y.man)
  val zeroed    = (minex +& exBias.U(exW.W)) <= (maxex +& diffexDec.asUInt)

  val zex0 = Mux(x.inf && y.inf, exBias.U(exW.W),
             Mux((x.inf && !y.inf) || (!x.inf && y.inf) || (x.zero) || (y.zero), 0.U(exW.W),
                 ((minex +& exBias.U) - maxex) - diffexDec.asUInt))
  val zexComb = Mux(zeroed, 0.U(exW.W), zex0(exW - 1, 0))

  // ---- specials (NaN-boxed) ----
  val xpos      = x.sgn === 0.U
  val xneg      = x.sgn === 1.U
  val tooLargeX = (zeroed && !yIsLarger) || ( x.inf && !y.inf) || (!x.zero &&  y.zero)
  val tooLargeY = (zeroed &&  yIsLarger) || (!x.inf &&  y.inf) || ( x.zero && !y.zero)
  val samexy    = x.ex === y.ex && x.man === y.man && !x.nan && !y.nan
  val zerozero  = x.zero && y.zero

  val znan      =  x.nan || y.nan
  val zzero     = (tooLargeX && xpos) || (zerozero && xpos)
  val zpi       = (tooLargeX && xneg) || (zerozero && xneg)
  val zhalfpi   =  tooLargeY
  val z1piover4 =  samexy && xpos
  val z3piover4 =  samexy && xneg

  val flagsComb = WireDefault(0.U.asTypeOf(new ATan2Flags))
  flagsComb.isSpecial := znan || zzero || zpi || zhalfpi || z1piover4 || z3piover4
  flagsComb.special   := MuxCase(ATan2SpecialValue.zNormal, Seq(
    znan      -> ATan2SpecialValue.zNaN,
    zzero     -> ATan2SpecialValue.zZero,
    zpi       -> ATan2SpecialValue.zPi,
    zhalfpi   -> ATan2SpecialValue.zHalfPi,
    z1piover4 -> ATan2SpecialValue.zQuarterPi,
    z3piover4 -> ATan2SpecialValue.z3QuarterPi,
  ))
  flagsComb.status := Cat(yIsLarger, x.sgn)
  flagsComb.ysgn   := y.sgn

  // ---- table address: ReciprocalPreProcess pattern (top adrW bits of max.man) ----
  val maxManComb = Mux(yIsLarger, y.man, x.man)
  val minManComb = Mux(yIsLarger, x.man, y.man)

  val xySameManComb = x.man === y.man

  // --------------------------------------------------------------------------
  // stage "ex": the ratio exponent and the flags

  private val cutEx = new SharedPipe("ex", cutCmp.valid, stages.at("ex"))

  val zex       = cutEx(zexComb)
  val flags     = cutEx(flagsComb)
  val maxMan    = cutEx(maxManComb)
  val minMan    = cutEx(minManComb)
  val xySameMan = cutEx(xySameManComb)

  // --------------------------------------------------------------------------
  // stage "out"

  io.flags := Pipe(cutEx.valid, flags, stages.at("out"))

  val adr = maxMan(manW - 1, dxW)
  io.adr := Pipe(cutEx.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = Cat(~maxMan(dxW - 1), maxMan(dxW - 2, 0))
    io.dx.get := Pipe(cutEx.valid, dx, stages.at("out"))
  }

  // ---- info bundle to post-process ----
  val info = Wire(new ATan2Phase1Info(spec))
  info.zex       := zex
  info.xySameMan := xySameMan
  info.maxXYMan0 := !maxMan.orR
  info.minMan    := minMan
  io.info := Pipe(cutEx.valid, info, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ATan2Phase1PreProcess]]. */
object ATan2Phase1PreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"cmp"` which of `|x|` and `|y|` is the larger
   *  - `"ex"`  the ratio exponent and the flags
   *  - `"out"` the table address, the info and the flags
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("cmp", 2.0),
    StagePoint("ex",  2.0),
    StagePoint("out", 1.0))
}

// Combinational driver of the post-process multiplier inputs.
//   lhs = (1 << fracW) + polyResult  -- builds (2 / (1+max_frac)) * 2^fracW with
//                                       hidden 1, gated to 0 when max.man == 0
//                                       (then min/max == min so the multiplier
//                                       short-circuits to (1+min.man))
//   rhs = Cat(1, minMan)             -- (1 + min.man)
class ATan2Phase1PostMultArgsOutput(val spec: RealSpec, val fracW: Int) extends Bundle {
  val lhs = UInt((1 + fracW).W)
  val rhs = UInt((1 + spec.manW).W)
}

class ATan2Phase1PostMultArgs(
  val spec: RealSpec,
  val bitW: PolynomialBitWidth,
) extends Module {
  val fracW = bitW.outBitW
  val manW  = spec.manW

  val io = IO(new Bundle {
    val zres      = Input(UInt(fracW.W))
    val maxXYMan0 = Input(Bool())
    val minMan    = Input(UInt(manW.W))
    val valid     = Input(Bool())
    val mul       = Valid(new ATan2Phase1PostMultArgsOutput(spec, fracW))
  })

  val zFrac = Mux(io.maxXYMan0, 0.U(fracW.W), io.zres)
  val out   = Wire(new ATan2Phase1PostMultArgsOutput(spec, fracW))
  out.lhs := Cat(1.U(1.W), zFrac)
  out.rhs := Cat(1.U(1.W), io.minMan)

  io.mul.valid := io.valid
  io.mul.bits  := out
}

class ATan2Phase1PostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ATan2Phase1PostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val info  = Flipped(Valid(new ATan2Phase1Info(spec)))
    val flags = Flipped(Valid(new ATan2Flags))
    val mul   = Flipped(Valid(new PostProcMultiplierOutput(manW)))
    val z     = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in",
    io.info.valid && io.flags.valid && io.mul.valid, stages.at("in"))

  val infoIn  = cutIn(io.info.bits)
  val flagsIn = cutIn(io.flags.bits)
  val mul     = cutIn(io.mul.bits)

  // The multiplier handles the ratio's mantissa with full rounding.  When the
  // mul output's exInc is asserted the ratio >= 2, but our exponent bookkeeping
  // already accounts for that by adding 1 to zex.  For phase 1 ratios <= 1 so
  // exInc should rarely fire; we still propagate it.
  val zManMul = mul.zman

  // ---- assemble the normal-case ratio ----
  val zManNormalComb = Mux(!infoIn.zex.orR || infoIn.xySameMan, 0.U(manW.W),
                       Mux(infoIn.maxXYMan0,                    infoIn.minMan,
                                                                zManMul))

  // --------------------------------------------------------------------------
  // stage "man": the ratio mantissa

  private val cutMan = new SharedPipe("man", cutIn.valid, stages.at("man"))

  val zManNormal = cutMan(zManNormalComb)
  val info       = cutMan(infoIn)
  val flags      = cutMan(flagsIn)

  // --------------------------------------------------------------------------
  // stage "out"

  val zexNoStat = info.zex  // exW bits, but always has MSB = 0 (we asserted that exBias's MSB is 0)
  val zSgnEnc = flags.ysgn
  val zExEnc  = Cat(flags.status(1), zexNoStat(exW - 2, 0))
  val zManEnc = Cat(zManNormal(manW - 1, 1), flags.status(0))

  // ---- special-value encoding ----
  val zSgnSpecial = flags.ysgn
  val zExSpecial  = Fill(exW, 1.U(1.W))
  val zManSpecial = Mux(flags.special === ATan2SpecialValue.zNaN,
                        Cat(1.U(1.W), 0.U((manW - 1).W)),
                        flags.special)

  val zSpecial = Cat(zSgnSpecial, zExSpecial, zManSpecial)
  val zNormal  = Cat(zSgnEnc,     zExEnc,     zManEnc)

  val z = Mux(flags.isSpecial, zSpecial, zNormal)
  io.z := Pipe(cutMan.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ATan2Phase1PostProcess]]. */
object ATan2Phase1PostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"man"` the ratio mantissa
   *  - `"out"` the encoded ratio
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("man", 1.0),
    StagePoint("out", 1.0))
}

// =========================================================================
//   Phase 2
// =========================================================================

class ATan2Phase2PreProcess(
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
  val stagePoints = ATan2Phase2PreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x      = Flipped(Valid(new DecomposedReal(spec)))
    val adr    = Valid(UInt(adrW.W))
    val dx     = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    val flags  = Valid(new ATan2Flags)
    // pipelined helpers consumed by the post-mult args / post-process
    val xManW1 = Valid(UInt((1 + manW).W))
    val xex    = Valid(UInt(exW.W))
    val isLinear = Valid(Bool())
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  // ---- decode NaN-boxed flags from phase-1 output ----
  val isSpecial = x.ex === Fill(exW, 1.U(1.W))
  // status0 is the LSB of x.man; status1 is the MSB of x.ex.
  val status = Cat(x.ex(exW - 1), x.man(0))

  val flagsComb = Wire(new ATan2Flags)
  flagsComb.isSpecial := isSpecial
  flagsComb.ysgn      := x.sgn
  flagsComb.status    := status
  flagsComb.special   := x.man(ATan2SpecialValue.width - 1, 0)

  // ---- decode the ratio t (strip status1 from ex, status0 from man) ----
  // If isSpecial, force xex to 0 so subsequent shifts don't underflow.
  val xexDecodedComb = Mux(isSpecial, 0.U(exW.W), Cat(0.U(1.W), x.ex(exW - 2, 0)))
  // xmanW1 := (1 + decoded_man) | 1  -- the OR-with-1 matches SimGeneric's
  // `(xmanW1 | SafeLong(1))` (avoids hitting an exact-zero address corner case).
  val xManW1 = Cat(1.U(1.W), x.man(manW - 1, 1), 1.U(1.W))

  // ---- shift x.manW1 into a fixed-point address ----
  val shiftW    = log2Up(2 + manW)
  val xShift0   = exBias.U(exW.W) - xexDecodedComb
  val xShiftOut = xShift0(exW - 1, shiftW).orR
  val xShift    = Mux(xShiftOut, (manW + 1).U(shiftW.W),
                      Mux(xShift0(shiftW - 1, 0) > (manW + 1).U(shiftW.W),
                          (manW + 1).U(shiftW.W),
                          xShift0(shiftW - 1, 0)))

  val xFixedComb = (xManW1 >> xShift)(manW, 0)

  // ---- linear-region detection: atan(t)/t ~ 1 when t is small enough ----
  val linearThreshold = (ATan2Phase2SimGeneric.calcLinearThreshold(manW) + exBias)
  val isLinearComb    = xexDecodedComb <= linearThreshold.U(exW.W)

  // ---- pipe x's mantissa-with-hidden-1 (used by post-mult) and decoded ex ----
  // post-mult rhs is (1+xman) where xman = x.man with status bit stripped (the
  // bottom bit cleared); we deliberately do NOT use the OR-with-1 form here.
  val postRhsComb = Cat(1.U(1.W), x.man(manW - 1, 1), 0.U(1.W))

  // --------------------------------------------------------------------------
  // stage "shift": the ratio aligned at the binary point

  private val cutShift = new SharedPipe("shift", cutIn.valid, stages.at("shift"))

  val xFixed     = cutShift(xFixedComb)
  val xexDecoded = cutShift(xexDecodedComb)
  val isLinear   = cutShift(isLinearComb)
  val postRhs    = cutShift(postRhsComb)
  val flags      = cutShift(flagsComb)

  // --------------------------------------------------------------------------
  // stage "out"

  io.flags := Pipe(cutShift.valid, flags, stages.at("out"))

  val adr = xFixed(manW - 1, dxW)
  io.adr := Pipe(cutShift.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = Cat(~xFixed(dxW - 1), xFixed(dxW - 2, 0))
    io.dx.get := Pipe(cutShift.valid, dx, stages.at("out"))
  }

  io.isLinear := Pipe(cutShift.valid, isLinear, stages.at("out"))

  io.xManW1 := Pipe(cutShift.valid, postRhs,    stages.at("out"))
  io.xex    := Pipe(cutShift.valid, xexDecoded, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ATan2Phase2PreProcess]]. */
object ATan2Phase2PreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"shift"` the ratio aligned at the binary point
   *  - `"out"`   the table address and the decoded ratio
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("shift", 3.0),
    StagePoint("out",   1.0))
}

class ATan2Phase2TableCoeff(
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

// Combinational driver of the post-process multiplier inputs for phase 2.
//   lhs = polyResult << 1 (with the table's top bit becoming the hidden 1)
//         or, in the linear region, 1 << fracW (representing atan(t)/t ~ 1).
//   rhs = (1 + decoded x.man).
class ATan2Phase2PostMultArgsOutput(val spec: RealSpec, val fracW: Int) extends Bundle {
  val lhs = UInt((1 + fracW).W)
  val rhs = UInt((1 + spec.manW).W)
}

class ATan2Phase2PostMultArgs(
  val spec: RealSpec,
  val bitW: PolynomialBitWidth,
) extends Module {
  val fracW = bitW.outBitW
  val manW  = spec.manW

  val io = IO(new Bundle {
    val zres     = Input(UInt(fracW.W))
    val isLinear = Input(Bool())
    val xManW1   = Flipped(Valid(UInt((1 + manW).W)))
    val mul      = Valid(new ATan2Phase2PostMultArgsOutput(spec, fracW))
  })

  // table value is always >= 0.5 (atan(t)/t ranges from 0.785 to 1-eps over
  // t in [0,1)), so the polynomial result's MSB is always 1; shifting left
  // by one bit puts that 1 into the hidden-bit position.
  val lhsTable  = Cat(io.zres, 0.U(1.W))
  val lhsLinear = Cat(1.U(1.W), 0.U(fracW.W))
  val out       = Wire(new ATan2Phase2PostMultArgsOutput(spec, fracW))
  out.lhs := Mux(io.isLinear, lhsLinear, lhsTable)
  out.rhs := io.xManW1.bits

  io.mul.valid := io.xManW1.valid
  io.mul.bits  := out
}

class ATan2Phase2PostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ATan2Phase2PostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val flags    = Flipped(Valid(new ATan2Flags))
    val xex      = Flipped(Valid(UInt(exW.W)))
    val isLinear = Flipped(Valid(Bool()))
    val mul      = Flipped(Valid(new PostProcMultiplierOutput(manW)))
    val z        = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in",
    io.flags.valid && io.xex.valid && io.isLinear.valid && io.mul.valid,
    stages.at("in"))

  val mul        = cutIn(io.mul.bits)
  val flagsIn    = cutIn(io.flags.bits)
  val xex        = cutIn(io.xex.bits)
  val isLinearIn = cutIn(io.isLinear.bits)
  val atanMan0   = mul.zman
  val atanInc    = mul.exInc

  // ---- atan(t) = (atan(t)/t) * t ----
  // SimGeneric: atanEx = xex + resEx + atanProdMoreThan2 + atanRoundMoreThan2
  //             where resEx = 0 (linear) or -1 (table)
  val atanEx0 = (xex +& atanInc) - Mux(isLinearIn, 0.U(1.W), 1.U(1.W))
  val atanExComb  = atanEx0(exW - 1, 0)
  val atanManComb = atanMan0

  // --------------------------------------------------------------------------
  // stage "ex": the exponent and mantissa of atan(t)

  private val cutEx = new SharedPipe("ex", cutIn.valid, stages.at("ex"))

  val atanEx   = cutEx(atanExComb)
  val atanMan  = cutEx(atanManComb)
  val flagsEx  = cutEx(flagsIn)

  // ---- align atan to the pi / pi-2 layout (manW + 4 bits: "1m.mmm...m mm") ----
  // pi   = 1x.xxxx....x000   (ex 1)
  // pi/2 = 01.xxxx....xx00   (ex 0)
  // atan = 01.xxxx....xx00   (ex 0, but shifted down by `exBias - atanEx`)
  val piManW1     = Real.pi(manW + 2).toBigInt.U((manW + 4).W)
  val halfPiManW1 = (Real.pi / Real.two)(manW + 2).toBigInt.U((manW + 4).W)
  val atanManW1   = Cat(1.U(2.W), atanMan, 0.U(2.W))

  assert(piManW1.getWidth     == manW + 4)
  assert(halfPiManW1.getWidth == manW + 4)
  assert(atanManW1.getWidth   == manW + 4)

  val shiftW      = log2Up(manW + 5)
  val atanShift0  = exBias.U(exW.W) - atanEx
  val atanShift   = Mux(atanShift0(exW - 1, shiftW).orR,
                        (manW + 1 + 3).U(shiftW.W),
                        atanShift0(shiftW - 1, 0))
  val atanAlignedComb = Mux(atanShift > (manW + 1 + 3).U(shiftW.W), 0.U((manW + 4).W),
                            atanManW1 >> atanShift)

  // --------------------------------------------------------------------------
  // stage "align": atan(t) aligned with pi

  private val cutAlign = new SharedPipe("align", cutEx.valid, stages.at("align"))

  val atanAligned  = cutAlign(atanAlignedComb)
  val atanExAlign  = cutAlign(atanEx)
  val atanManAlign = cutAlign(atanMan)
  val flagsAlign   = cutAlign(flagsEx)

  // ---- status 1: ysgn * (pi - atan(|y|/|x|))  (x<0, |x|>|y|) ----
  // pi - atan in (3pi/4, pi] ~ (2.35, 3.14], so the result has ex = 1.
  val piMinusATanComb = piManW1 - atanAligned

  // ---- status 2: ysgn * (pi/2 - atan(|x|/|y|)) (x>0, |x|<|y|) ----
  // pi/2 - atan in (pi/4, pi/2] ~ (0.78, 1.57], so ex is exBias-1 or exBias.
  val halfPiMinusATanComb = halfPiManW1 - atanAligned

  // ---- status 3: ysgn * (pi/2 + atan(|x|/|y|)) (x<0, |x|<|y|) ----
  // pi/2 + atan in (pi/2, 3pi/4] ~ (1.57, 2.35], so ex is exBias or exBias+1.
  val halfPiPlusATanComb = halfPiManW1 + atanAligned

  // --------------------------------------------------------------------------
  // stage "diff": pi - atan, pi/2 - atan and pi/2 + atan

  private val cutDiff = new SharedPipe("diff", cutAlign.valid, stages.at("diff"))

  val piMinusATan     = cutDiff(piMinusATanComb)
  val halfPiMinusATan = cutDiff(halfPiMinusATanComb)
  val halfPiPlusATan  = cutDiff(halfPiPlusATanComb)
  val atanExDiff      = cutDiff(atanExAlign)
  val atanManDiff     = cutDiff(atanManAlign)
  val flagsDiff       = cutDiff(flagsAlign)

  // ---- status 0: ysgn * atan(|y|/|x|)         (x>0, |x|>|y|) ----
  val zExPL  = atanExDiff
  val zManPL = atanManDiff(manW - 1, 0)

  val zExNL  = (exBias + 1).U(exW.W)
  val zManNL = piMinusATan((manW + 4) - 2, 3) + piMinusATan(2)

  val halfPiMinusATanMoreThan1 = halfPiMinusATan((manW + 4) - 2)
  val zExPS  = (exBias - 1).U(exW.W) + halfPiMinusATanMoreThan1
  val zManPS = Mux(halfPiMinusATanMoreThan1,
                   halfPiMinusATan((manW + 4) - 3, 2) + halfPiMinusATan(1),
                   halfPiMinusATan((manW + 4) - 4, 1) + halfPiMinusATan(0))

  val halfPiPlusATanRounded  = dropLSB(2, halfPiPlusATan) + halfPiPlusATan(1)
  val halfPiPlusATanMoreThan2 = halfPiPlusATanRounded(manW + 1)
  val zExNS  = exBias.U(exW.W) + halfPiPlusATanMoreThan2
  val zManNS = Mux(halfPiPlusATanMoreThan2,
                   halfPiPlusATanRounded(manW,     1),
                   halfPiPlusATanRounded(manW - 1, 0))

  val zExNormalComb = MuxCase(0.U(exW.W), Seq(
    (flagsDiff.status === ATan2Status.xIsPosIsLarger ) -> zExPL,
    (flagsDiff.status === ATan2Status.xIsNegIsLarger ) -> zExNL,
    (flagsDiff.status === ATan2Status.xIsPosIsSmaller) -> zExPS,
    (flagsDiff.status === ATan2Status.xIsNegIsSmaller) -> zExNS,
  ))
  val zManNormalComb = MuxCase(0.U(manW.W), Seq(
    (flagsDiff.status === ATan2Status.xIsPosIsLarger ) -> zManPL(manW - 1, 0),
    (flagsDiff.status === ATan2Status.xIsNegIsLarger ) -> zManNL(manW - 1, 0),
    (flagsDiff.status === ATan2Status.xIsPosIsSmaller) -> zManPS(manW - 1, 0),
    (flagsDiff.status === ATan2Status.xIsNegIsSmaller) -> zManNS(manW - 1, 0),
  ))

  // --------------------------------------------------------------------------
  // stage "sel": the selected normal-case exponent and mantissa

  private val cutSel = new SharedPipe("sel", cutDiff.valid, stages.at("sel"))

  val zExNormal  = cutSel(zExNormalComb)
  val zManNormal = cutSel(zManNormalComb)
  val flagsSel   = cutSel(flagsDiff)

  // --------------------------------------------------------------------------
  // stage "out"

  // ---- special values ----
  val realPi         = new RealGeneric(spec, Pi)
  val realHalfPi     = new RealGeneric(spec, Pi * 0.5)
  val realQuarterPi  = new RealGeneric(spec, Pi * 0.25)
  val realQuarter3Pi = new RealGeneric(spec, Pi * 0.75)

  val zExSpecial = MuxCase(0.U(exW.W), Seq(
    (flagsSel.special === ATan2SpecialValue.zNaN       ) -> Fill(exW, 1.U(1.W)),
    (flagsSel.special === ATan2SpecialValue.zZero      ) -> 0.U(exW.W),
    (flagsSel.special === ATan2SpecialValue.zPi        ) -> realPi        .ex.U(exW.W),
    (flagsSel.special === ATan2SpecialValue.zHalfPi    ) -> realHalfPi    .ex.U(exW.W),
    (flagsSel.special === ATan2SpecialValue.zQuarterPi ) -> realQuarterPi .ex.U(exW.W),
    (flagsSel.special === ATan2SpecialValue.z3QuarterPi) -> realQuarter3Pi.ex.U(exW.W),
  ))
  val zManSpecial = MuxCase(0.U(manW.W), Seq(
    (flagsSel.special === ATan2SpecialValue.zNaN       ) -> Cat(1.U(1.W), 0.U((manW - 1).W)),
    (flagsSel.special === ATan2SpecialValue.zZero      ) -> 0.U(manW.W),
    (flagsSel.special === ATan2SpecialValue.zPi        ) -> realPi        .man.toBigInt.U(manW.W),
    (flagsSel.special === ATan2SpecialValue.zHalfPi    ) -> realHalfPi    .man.toBigInt.U(manW.W),
    (flagsSel.special === ATan2SpecialValue.zQuarterPi ) -> realQuarterPi .man.toBigInt.U(manW.W),
    (flagsSel.special === ATan2SpecialValue.z3QuarterPi) -> realQuarter3Pi.man.toBigInt.U(manW.W),
  ))

  // ---- assemble ----
  // NaN's sign is 0 (per RealGeneric.nan); all other paths carry ysgn.
  val zSgnSpecial = Mux(flagsSel.special === ATan2SpecialValue.zNaN, 0.U(1.W),
                        flagsSel.ysgn)

  val zSgn = Mux(flagsSel.isSpecial, zSgnSpecial, flagsSel.ysgn)
  val zEx  = Mux(flagsSel.isSpecial, zExSpecial,  zExNormal)
  val zMan = Mux(flagsSel.isSpecial, zManSpecial, zManNormal)

  val z = Cat(zSgn, zEx, zMan)
  io.z := Pipe(cutSel.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ATan2Phase2PostProcess]]. */
object ATan2Phase2PostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"ex"`    the exponent and mantissa of atan(t)
   *  - `"align"` atan(t) aligned with pi
   *  - `"diff"`  pi - atan, pi/2 - atan and pi/2 + atan
   *  - `"sel"`   the selected normal-case exponent and mantissa
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("ex",    2.0),
    StagePoint("align", 2.0),
    StagePoint("diff",  2.0),
    StagePoint("sel",   2.0),
    StagePoint("out",   1.0))
}
