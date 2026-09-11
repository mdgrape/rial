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
package rial.math

import scala.math.{max => smax, min => smin}

import spire.math.Real
import spire.implicits._

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{RealSpec, FloatChiselUtil}
import rial.util._
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

  /** The product of the mantissa by the constant one, in the field [[prod]]
   *  occupies: what the pre-process reads for an argument the caller has
   *  already divided by π. The pre-process takes the exponent of such a
   *  product from the argument itself, see `scaled` of
   *  [[rial.math.SinCosPreProcess]].
   *
   *  @param spec floating point spec of the argument
   *  @param man  mantissa of the argument, without its hidden one
   */
  def prodOfOne(spec: RealSpec, man: UInt): UInt = {
    val product = Cat(0.U(1.W), 1.U(1.W), man, 0.U((lhsW(spec) - 1).W))
    assert(product.getWidth == prodW(spec),
      f"prodOfOne(${product.getWidth}) must be prodW(${prodW(spec)}) wide")
    product
  }
}

// -------------------------------------------------------------------------
//
// the range of the argument that the reduction covers
//

/** The range of |x| that the sin/cos argument reduction is built for.
 *
 *  The reduction folds |x|/π into [0, 2) by shifting it by its own
 *  exponent, so the range of |x| it has to cover sets the width of that
 *  shift. `maxEx` is the largest unbiased exponent of |x| the reduction
 *  covers: it is built for every finite |x| < 2^(maxEx+1), down to the
 *  argument at which |x|/π underflows (below that sin(x) is ±0 and
 *  cos(x) is ±1). For an |x| outside the range the result is a
 *  meaningless, but finite, value in [-1, 1].
 *
 *  A narrower range shortens the segments of the reduction that carry the
 *  fold and the normalisation of what it leaves: `"pre.sinCos.norm" ->
 *  "pre.sinCos.fold"`, which aligns |x|/π over fewer bits,
 *  `"pre.sinCos.fold" -> "pre.sinCos.lead"`, which searches a shorter
 *  window for the leading one, and `"pre.sinCos.lead" ->
 *  "pre.sinCos.shift"`, which shifts by a narrower amount; it changes
 *  neither the boundaries nor the latency of
 *  [[rial.math.SinCosPreProcess]].
 *
 *  {{{
 *  val cfg = MathFuncConfig.standardFP32.restrictSinCos(
 *    Some(SinCosDomain.twoPi(RealSpec.Float32Spec)))
 *  }}}
 *
 *  @param maxEx the largest unbiased exponent of |x| to cover
 */
case class SinCosDomain(maxEx: Int) {

  /** The largest number of bits the fold shifts |x|/π up by.
   *
   *  |x| < 2^(maxEx+1) and π > 2, so |x|/π < 2^maxEx.
   */
  def maxFoldShift: Int = {
    smax(0, maxEx - 1)
  }

  /** The largest shift the normalisation of the folded argument needs, when
   *  the range is narrow enough for that shift to be known.
   *
   *  The fold cancels the leading bits of |x|/π when x is close to a multiple
   *  of π/2, and the normalisation shifts what is left back up. How far the
   *  cancellation reaches is a property of the binary expansion of the 1/π
   *  constant over the arguments the range admits, so this is a checked
   *  bound, not a derived one: every mantissa of the two formats named here
   *  has been walked, at every fold shift a range up to the second period
   *  admits, and none cancels more than `manW + 4` bits. A range that admits
   *  more, or a format that has not been walked, gives `None` and the whole
   *  folded argument is searched.
   */
  def maxNormShiftOf(spec: RealSpec): Option[Int] = {
    val walked = Seq(RealSpec.Float32Spec, RealSpec.BFloat16Spec)
    if (maxEx <= 2 && walked.exists(_ == spec)) {
      Some(spec.manW + 4)
    } else {
      None
    }
  }
}

/** Factory for [[rial.math.SinCosDomain]]. */
object SinCosDomain {

  /** Every finite argument of the format. */
  def full(spec: RealSpec): SinCosDomain = {
    SinCosDomain(spec.exMax)
  }

  /** |x| up to the largest power of two above 2π. */
  def twoPi(spec: RealSpec): SinCosDomain = {
    SinCosDomain(smin(2, spec.exMax))
  }

  /** |x| up to the largest power of two above π. */
  def pi(spec: RealSpec): SinCosDomain = {
    SinCosDomain(smin(1, spec.exMax))
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
  val spec:   RealSpec,
  val bitW:   PolynomialBitWidth,
  val stage:  PipelineStageConfig,
  val domain: Option[SinCosDomain] = None,
  // whether the caller has already divided the argument by π, in which case
  // `io.prod` carries the product of the mantissa by the constant one
  val scaled: Boolean = false,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  /** The range of the argument the reduction is built for. */
  val range = domain.getOrElse(SinCosDomain.full(spec))

  val adrW   = bitW.addrW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = SinCosPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val isSin = Input(Bool())
    val x     = Flipped(Valid(new DecomposedReal(spec)))
    val prod  = Flipped(Valid(UInt(SinCosPreMulArgs.prodW(spec).W)))
    val adr   = Valid(UInt(adrW.W))
    val dx    = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    val info  = Valid(new SinCosPreProcInfo(spec))
  })

  // --------------------------------------------------------------------------
  // stage "in": both the argument and the product are needed

  private val cutIn = new SharedPipe("in", io.x.valid && io.prod.valid,
    stages.at("in"))

  val isSinIn = cutIn(io.isSin)
  val x       = cutIn(io.x.bits)
  val prodIn  = cutIn(io.prod.bits)

  val xsgn = x.sgn
  val xex  = x.ex
  val xman = x.man

  val znanIn = xex.andR

  // --------------------------------------------------------------------------
  // x/π from the pre-mult output

  val oneOverPiW = SinCosPreMulArgs.lhsW(spec)

  val xOverPiProd          = Mux(xex === 0.U, 0.U, prodIn)
  val xOverPiProdMoreThan2 = xOverPiProd((1 + manW) + oneOverPiW - 1)
  // With the constant one the product is the mantissa itself and x/π is x,
  // so the exponent is the one of x.
  val xOverPiExComb = if (scaled) {
    Cat(0.U(1.W), xex)
  } else {
    xex -& 2.U + xOverPiProdMoreThan2
  }
  // x.ex > xOverPiEx always (xOverPiProdMoreThan2 ≤ 1); overflow here just
  // means x/π is tiny, not infinite.

  val xOverPiComb  = Mux(xOverPiProdMoreThan2, xOverPiProd,
                                               Cat(xOverPiProd.tail(1), 0.U(1.W)))
  val xOverPiFracW = (1 + manW) + oneOverPiW - 1

  // --------------------------------------------------------------------------
  // stage "norm": the normalized x/π

  private val cutNorm = new SharedPipe("norm", cutIn.valid, stages.at("norm"))

  val xOverPi   = cutNorm(xOverPiComb)
  val xOverPiEx = cutNorm(xOverPiExComb)
  val xsgnNorm  = cutNorm(xsgn)
  val isSinNorm = cutNorm(isSinIn)
  val znanNorm  = cutNorm(znanIn)

  // |x|/π underflows: `xOverPiEx` is a biased exponent one bit wider than
  // the format's, so the argument is below the smallest normal |x|/π when
  // the subtraction above borrowed or landed on the zero exponent.
  val xOverPiIsZero = xOverPiEx(exW) || xOverPiEx === 0.U
  assert(xOverPi.getWidth == xOverPiFracW + 1)

  // --------------------------------------------------------------------------
  // fold x/π into [0, 2)

  // The up-shift is bounded by the range the reduction covers, the down-shift
  // by the point past which every bit the rounding below reads is zero: both
  // shift amounts are held at that width and saturate there.
  val foldShiftW = smin(log2Ceil(range.maxFoldShift + 1), log2Ceil(xOverPiFracW))
  val downShiftW = smin(log2Ceil(manW + 5), log2Ceil(xOverPiFracW))
  assert(downShiftW <= exW, f"downShiftW(${downShiftW}) must be <= exW(${exW})")

  val xOverPiExMoreThan1 = (xOverPiEx >= exBias.U)

  val xOverPiExNobiasNegFull = (exBias.U - xOverPiEx)
  val xOverPiExNobiasNeg     = Mux(
    xOverPiExNobiasNegFull(exW, downShiftW).orR,
    Fill(downShiftW, 1.U(1.W)),
    xOverPiExNobiasNegFull(downShiftW - 1, 0))

  val xOverPiAlignPos = if (foldShiftW == 0) {
    xOverPi
  } else {
    val xOverPiExNobiasPos = (xOverPiEx - exBias.U)(foldShiftW - 1, 0)
    (xOverPi << xOverPiExNobiasPos)(xOverPiFracW, 0)
  }
  val xOverPiAlignNeg = (xOverPi >> xOverPiExNobiasNeg)(xOverPiFracW, 0)
  val xOverPiAligned  = Mux(xOverPiIsZero, 0.U((xOverPiFracW + 1).W),
                        Mux(xOverPiExMoreThan1, xOverPiAlignPos, xOverPiAlignNeg))

  // A right shift by two or more clears the two top bits, so the branch the
  // exponent below takes needs a two-bit window of `xOverPi` rather than the
  // whole shifter; this is `xOverPiAlignNeg.head(2)`.
  val xOverPiAlignNeg2MSBs =
    Mux(xOverPiExNobiasNeg === 0.U, xOverPi.head(2),
    Mux(xOverPiExNobiasNeg === 1.U, Cat(0.U(1.W), xOverPi.head(1)), 0.U(2.W)))

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

  val zsgnComb = Mux(isSinNorm, xOverPiAligned2MSBs(1) ^ xsgnNorm,
                                xOverPiAligned2MSBs(1) ^ xOverPiAligned2MSBs(0))

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
  val ymanAligned = Mux(xOverPiAligned2MSBs(0) ^ isSinNorm, ymanPos, ymanNeg)

  // For 0 < x < π/2 we use x/π itself (un-aligned) to keep precision.
  val ymanAsIs = xOverPi(xOverPi.getWidth - 1, 2)
  val yexAsIs  = xOverPiEx
  val nonAlign = !xOverPiIsZero && !xOverPiExMoreThan1 &&
                 xOverPiAlignNeg2MSBs === 0.U && isSinNorm

  assert(ymanAsIs.getWidth == ymanPos.getWidth)
  val yman0Comb = Mux(nonAlign, ymanAsIs, ymanAligned)
  val yex0Comb  = Mux(nonAlign, yexAsIs,  (exBias - 2).U(exW.W))

  // --------------------------------------------------------------------------
  // stage "fold": the folded argument

  private val cutFold = new SharedPipe("fold", cutNorm.valid, stages.at("fold"))

  val yman0    = cutFold(yman0Comb)
  val yex0     = cutFold(yex0Comb)
  val zsgnFold = cutFold(zsgnComb)
  val znanFold = cutFold(znanNorm)

  val ymanIsNonZeroComb = yman0.orR

  // normalize yman0: shift its leading one up, and take the shift out of the
  // exponent.
  val yman0W       = yman0.getWidth
  val maxNormShift = range.maxNormShiftOf(spec).getOrElse(yman0W - 1)
  val yman0ShiftW  = log2Ceil(maxNormShift + 1)
  val yman0Shift0  = PriorityEncoder(Reverse(yman0.head(maxNormShift + 1)))
  val yman0Shift   = yman0Shift0(yman0ShiftW - 1, 0)

  // --------------------------------------------------------------------------
  // stage "lead": the position of the leading one

  private val cutLead = new SharedPipe("lead", cutFold.valid, stages.at("lead"))

  val yman0ShiftLead    = cutLead(yman0Shift)
  val yman0Lead         = cutLead(yman0)
  val yex0Lead          = cutLead(yex0)
  val ymanIsNonZeroLead = cutLead(ymanIsNonZeroComb)
  val zsgnLead          = cutLead(zsgnFold)
  val znanLead          = cutLead(znanFold)

  val yman0ShiftedComb = (yman0Lead << yman0ShiftLead)(yman0W - 1, 0)
  val yexShiftedComb   = yex0Lead - yman0ShiftLead

  // --------------------------------------------------------------------------
  // stage "shift": the normalized argument

  private val cutShift = new SharedPipe("shift", cutLead.valid, stages.at("shift"))

  val yman0Shifted  = cutShift(yman0ShiftedComb)
  val yexShifted    = cutShift(yexShiftedComb)
  val ymanIsNonZero = cutShift(ymanIsNonZeroLead)
  val zsgnShift     = cutShift(zsgnLead)
  val znanShift     = cutShift(znanLead)

  // round it to a manW-bit mantissa, like the old code.
  val yman0RoundBit = yman0W - manW - 1
  val yman0Rounded  = yman0Shifted(yman0W - 2, yman0RoundBit) +&
                      yman0Shifted(yman0RoundBit - 1)
  val yman0MoreThan2 = yman0Rounded(manW)

  val ymanNonZero = yman0Rounded(manW - 1, 0)
  val yexNonZero  = yexShifted + yman0MoreThan2
  val yexComb  = Mux(!ymanIsNonZero, 0.U, yexNonZero)
  val ymanComb = Mux(!ymanIsNonZero, 0.U, ymanNonZero)

  val zoneComb = (yexComb === (exBias - 1).U) && ymanComb === 0.U

  // --------------------------------------------------------------------------
  // stage "round": the rounded argument

  private val cutRound = new SharedPipe("round", cutShift.valid, stages.at("round"))

  val yex       = cutRound(yexComb)
  val yman      = cutRound(ymanComb)
  val zone      = cutRound(zoneComb)
  val zsgnRound = cutRound(zsgnShift)
  val znanRound = cutRound(znanShift)

  // --------------------------------------------------------------------------
  // address / dx for the polynomial

  val yaligned = Cat(1.U(1.W), yman) >> ((exBias - 1).U(exW.W) - yex)
  val adr      = yaligned(manW - 1, manW - adrW)

  // --------------------------------------------------------------------------
  // stage "out"

  val infoOut = Wire(new SinCosPreProcInfo(spec))
  infoOut.znan := znanRound
  infoOut.zone := zone
  infoOut.zsgn := zsgnRound
  infoOut.yex  := yex
  infoOut.yman := yman

  io.info := Pipe(cutRound.valid, infoOut, stages.at("out"))
  io.adr  := Pipe(cutRound.valid, adr,     stages.at("out"))
  if (order != 0) {
    val dx = Cat(~yaligned(manW - adrW - 1), yaligned(manW - adrW - 2, 0))
    io.dx.get := Pipe(cutRound.valid, dx, stages.at("out"))
  }
}

/** Pipeline stage points of [[rial.math.SinCosPreProcess]]. */
object SinCosPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"norm"`  the normalized x/π
   *  - `"fold"`  the folded argument
   *  - `"lead"`  the position of the leading one of the folded argument
   *  - `"shift"` the normalized argument
   *  - `"round"` the rounded argument
   *  - `"out"`   the table address and the info
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("norm",  3.0),
    StagePoint("fold",  4.0),
    StagePoint("lead",  1.0),
    StagePoint("shift", 1.0),
    StagePoint("round", 2.0),
    StagePoint("out",   2.0))
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
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = SinCosPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val pre   = Flipped(Valid(new SinCosPreProcInfo(spec)))
    val mul   = Flipped(Valid(new PostProcMultiplierOutput(manW)))
    val z     = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.pre.valid && io.mul.valid,
    stages.at("in"))

  val info  = cutIn(io.pre.bits)
  val mul   = cutIn(io.mul.bits)

  val zman0Comb = mul.zman
  val zex0Comb  = info.yex + 1.U + mul.exInc

  val znanComb  = info.znan
  val zoneComb  = info.zone
  val zzeroComb = info.yex === 0.U
  val zsgnComb  = info.zsgn

  // --------------------------------------------------------------------------
  // stage "ex": the result exponent

  private val cutEx = new SharedPipe("ex", cutIn.valid, stages.at("ex"))

  val zex0  = cutEx(zex0Comb)
  val zman0 = cutEx(zman0Comb)
  val znan  = cutEx(znanComb)
  val zone  = cutEx(zoneComb)
  val zzero = cutEx(zzeroComb)
  val zsgn  = cutEx(zsgnComb)

  // --------------------------------------------------------------------------
  // stage "out"

  val zex = Mux(zzero, 0.U,
            Mux(znan,  Fill(exW, 1.U(1.W)),
            Mux(zone,  exBias.U(exW.W),
                       zex0)))
  val zman = Mux(zzero, 0.U,
             Mux(znan || zone, Cat(znan, 0.U((manW - 1).W)),
                               zman0))

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := Pipe(cutEx.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.SinCosPostProcess]]. */
object SinCosPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"ex"`  the result exponent
   *  - `"out"` the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("ex",  1.0),
    StagePoint("out", 1.0))
}
