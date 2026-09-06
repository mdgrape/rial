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
import rial.util._
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
  val stage: PipelineStageConfig,
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

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ExpPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val prod = Flipped(Valid(UInt(ExpPreMulArgs.prodW(spec).W)))
    val adr  = Valid(UInt(adrW.W))
    val dx   = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    val info = Valid(new ExpPreProcInfo(spec, padding))
  })

  // both the argument and the product are needed
  private val cutIn = new SharedPipe("in", io.x.valid && io.prod.valid,
    stages.at("in"))

  val x = cutIn(io.x.bits)

  val xsgn = x.sgn
  val xex0 = x.ex
  val xman = x.man

  // ------------------------------------------------------------------------
  // round the pre-product to (manW + extraMan) bits

  val log2eManW          = ExpPreMulArgs.lhsW(spec) - 1
  val xprod              = cutIn(io.prod.bits)
  val xprodW             = ExpPreMulArgs.prodW(spec)
  val xprodMoreThan2Comb = xprod(xprodW - 1) === 1.U

  assert(log2eManW - extraMan > 0)

  val xprodbp        = manW + log2eManW
  val xprodRoundBits = xprodbp - (manW + extraMan)

  val xprodSticky  = (xprodMoreThan2Comb & xprod(xprodRoundBits - 1)) |
                     xprod(xprodRoundBits - 2, 0).orR
  val xprodRound   = Mux(xprodMoreThan2Comb, xprod(xprodRoundBits),
                                             xprod(xprodRoundBits - 1))
  val xprodShifted = Mux(xprodMoreThan2Comb, xprod(xprodbp,     xprodRoundBits + 1),
                                             xprod(xprodbp - 1, xprodRoundBits))
  val xprodLSB     = xprodShifted(0)
  val xprodInc     = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
                                                    xprodLSB, xprodRound, xprodSticky)
  val xprodRoundedComb = xprodShifted +& xprodInc

  // ------------------------------------------------------------------------
  // stage "round": the rounded product x * log2(e)

  private val cutRound = new SharedPipe("round", cutIn.valid, stages.at("round"))

  val xprodRounded   = cutRound(xprodRoundedComb)
  val xprodMoreThan2 = cutRound(xprodMoreThan2Comb)
  val xexRound       = cutRound(xex0)
  val xsgnRound      = cutRound(xsgn)

  assert(xprodRounded.getWidth == manW + 1 + extraMan)
  val xprodMoreThan2AfterRounded = xprodRounded(manW + extraMan)

  val xman2 = xprodRounded(manW + extraMan - 1, 0)
  val xex   = xexRound + xprodMoreThan2AfterRounded + xprodMoreThan2

  assert(!cutRound.valid || xprodMoreThan2AfterRounded +& xprodMoreThan2 =/= 2.U)

  val xexdComb = xprodMoreThan2AfterRounded + xprodMoreThan2
  val xVal     = Cat(1.U(1.W), xman2)
  assert(xVal.getWidth == xValW)

  // ------------------------------------------------------------------------
  // split into integer/fractional parts

  val xshift0     = (xIntW + exBias).U(exW.W) - xex
  val xshift      = xshift0(xIntW - 1, 0)
  val xValShifted = Mux(xshift0 >= (xValW + 1).U, 0.U, xVal >> (xshift - 1.U(1.W)))
  val xValRoundedComb = (xValShifted >> 1) + xValShifted(0)

  // ------------------------------------------------------------------------
  // stage "shift": the argument aligned at the binary point

  private val cutShift = new SharedPipe("shift", cutRound.valid, stages.at("shift"))

  val xValRounded = cutShift(xValRoundedComb)
  val xexdShift   = cutShift(xexdComb)
  val xsgnShift   = cutShift(xsgnRound)

  val xint0  = xValRounded(xIntW + xFracW - 1, xFracW)
  val xfrac0 = xValRounded(xFracW - 1, 0)

  // (sign handling: ExpNonTable converts xint to a biased exponent.)
  val xintComb = xint0 +& (xsgnShift.asBool && (xfrac0 =/= 0.U)).asUInt

  val xfracNeg  = (BigInt(1) << xFracW).U((xFracW + 1).W) - xfrac0
  val xfracComb = Mux(xsgnShift === 0.U, xfrac0, xfracNeg(xFracW - 1, 0))

  // ------------------------------------------------------------------------
  // stage "neg": the fraction, negated for a negative argument

  private val cutNeg = new SharedPipe("neg", cutShift.valid, stages.at("neg"))

  val xint  = cutNeg(xintComb)
  val xfrac = cutNeg(xfracComb)
  val xexd  = cutNeg(xexdShift)

  // ------------------------------------------------------------------------
  // outputs

  val adr = xfrac(xFracW - 1, (xFracW - 1) - adrW + 1)

  val infoOut = Wire(new ExpPreProcInfo(spec, padding))
  infoOut.xint := xint
  infoOut.xexd := xexd
  if (padding != 0) {
    infoOut.xfracLSBs.get := xfrac(padding - 1, 0)
  }

  io.info := Pipe(cutNeg.valid, infoOut, stages.at("out"))
  io.adr  := Pipe(cutNeg.valid, adr,     stages.at("out"))
  if (order != 0) {
    val dx = Cat(~xfrac(xFracW - 1 - adrW), xfrac(xFracW - 1 - adrW - 1, padding))
    io.dx.get := Pipe(cutNeg.valid, dx, stages.at("out"))
  }
}

/** Pipeline stage points of [[rial.math.ExpPreProcess]]. */
object ExpPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"round"` the rounded product x * log2(e)
   *  - `"shift"` the argument aligned at the binary point
   *  - `"neg"`   the fraction, negated for a negative argument
   *  - `"out"`   the table address and the info
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("round", 3.0),
    StagePoint("shift", 3.0),
    StagePoint("neg",   2.0),
    StagePoint("out",   1.0))
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

  // The polynomial approximates 2^t-1 on [0, 1) with `padding` bits below the
  // mantissa LSB, and ExpPostProcess rounds that result to manW bits.  Should
  // the rounding carry into the hidden bit, the correction term in
  // ExpPostProcess would be built from the wrong mantissa, so require that the
  // table cannot reach the value at which the rounding carries.
  private val padding = coeff.bitW.outBitW - spec.manW

  if (padding > 0) {
    val maxPoly   = ExpTableCoeff.maxPolynomial(coeff)
    val threshold = (BigInt(1) << coeff.bitW.outBitW) - (BigInt(1) << (padding - 1))

    require(maxPoly < threshold,
      f"ExpTableCoeff: the polynomial reaches $maxPoly, at or above the " +
      f"$threshold at which rounding to ${spec.manW} bits carries into the " +
      f"hidden bit (nOrder=${coeff.bitW.nOrder}, addrW=$adrW, " +
      f"inBitW=${coeff.bitW.inBitW}, outBitW=${coeff.bitW.outBitW}, " +
      f"manW=${spec.manW}, padding=$padding). " +
      "Widen outBitW or refine the table so that the approximation of 2^t-1 " +
      "stays below the threshold.")
  }

  val io = IO(new Bundle {
    val adr = Flipped(Valid(UInt(adrW.W)))
    val cs  = Valid(UInt(coeff.bitW.coefBits.reduce(_+_).W))
  })

  val packed = VecInit(PolynomialGenerator.packTable(coeff).map(c => {
    c.U(coeff.bitW.coefBits.reduce(_+_).W)
  }))
  io.cs := Pipe(io.adr.valid, packed(io.adr.bits), stage)
}

/** Table checks of [[rial.math.ExpTableCoeff]]. */
object ExpTableCoeff {

  /** The largest value the polynomial of `coeff` can evaluate to.
   *
   *  It evaluates the table with [[rial.table.PolynomialEval.eval]], the
   *  software model of the hardware evaluator, over every `dx` of the last
   *  interval and the two extreme `dx` of every other interval.  The table
   *  approximates the increasing 2^t-1, so its maximum lives in the last
   *  interval; the endpoints of the others guard that assumption.
   *
   *  @param coeff the coefficients of the exp table
   *  @return the maximum of the integer polynomial result
   */
  def maxPolynomial(coeff: PolynomialCoefficients): BigInt = {
    val nInterval = coeff.coeffs.length
    val dxMax     = (BigInt(1) << coeff.bitW.inBitW) - 1

    val lastMax = (BigInt(0) to dxMax).map(dx => {
      PolynomialEval.eval(nInterval - 1, dx, coeff)
    }).max

    (0 until (nInterval - 1)).foldLeft(lastMax)( (acc, addr) => {
      acc.max(PolynomialEval.eval(addr, BigInt(0), coeff))
         .max(PolynomialEval.eval(addr, dxMax,     coeff))
    })
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

class ExpPostProcInfo(val spec: RealSpec, val padding: Int) extends Bundle {
  val zex         = UInt(spec.exW.W)
  val znan        = Bool()
  val zinf        = Bool()
  val zzero       = Bool()
  val zIsNonTable = Bool()
  val xfracLSBs   = if (padding != 0) Some(UInt(padding.W)) else None
}

class ExpNonTable(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val padding = bitW.outBitW - manW

  val zCorrTmpW = scala.math.min(6, 1 + manW + padding)
  val coeffW    = scala.math.min(8, manW + 1)

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ExpNonTable.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val info = Flipped(Valid(new ExpPreProcInfo(spec, padding)))
    val z    = Valid(new ExpPostProcInfo(spec, padding))
  })

  private val cutIn = new SharedPipe("in", io.x.valid && io.info.valid,
    stages.at("in"))

  val xIn    = cutIn(io.x.bits)
  val infoIn = cutIn(io.info.bits)

  // ------------------------------------------------------------------------
  // stage "mul": the correction coefficient

  private val cutMul = new SharedPipe("mul", cutIn.valid, stages.at("mul"))

  val zout = Wire(new ExpPostProcInfo(spec, padding))
  if (padding != 0) {
    zout.xfracLSBs.get := cutMul(infoIn.xfracLSBs.get)
  }

  val x    = cutMul(xIn)
  val info = cutMul(infoIn)

  // ------------------------------------------------------------------------
  // exponent + special-flag computation

  private val log2 = (a: Double) => log(a) / log(2.0)
  val xExOvfLimit = scala.math.ceil(log2(spec.exMax.toDouble)).toInt
  val xExUdfLimit = scala.math.ceil(log2((-spec.exMin).toDouble)).toInt

  val xExOvfLimBiased = (xExOvfLimit + exBias).U(exW.W)
  val xExUdfLimBiased = (xExUdfLimit + exBias).U(exW.W)

  val xint = info.xint
  val xexd = info.xexd

  val zexPos = xint +& exBias.U(exW.W)
  val zexNeg = exBias.U(exW.W) -& xint
  val zex0   = Mux(x.sgn === 0.U, zexPos(exW - 1, 0), zexNeg(exW - 1, 0))

  val xex   = x.ex + xexd
  val znan  = x.nan
  val zinf  = (x.sgn === 0.U) &&
              ((xExOvfLimBiased <= xex) || (  spec.exMax.U  < xint) || x.inf)
  val zzero = (x.sgn === 1.U) &&
              ((xExUdfLimBiased <= xex) || ((-spec.exMin).U < xint) || x.inf)

  zout.zex         := Mux(znan || zinf, Fill(exW, 1.U(1.W)),
                      Mux(zzero,        0.U(exW.W),
                                        zex0))
  zout.znan        := znan
  zout.zinf        := zinf
  zout.zzero       := zzero
  zout.zIsNonTable := znan || zinf || zzero

  io.z := Pipe(cutMul.valid, zout, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ExpNonTable]]. */
object ExpNonTable {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"mul"` the correction coefficient
   *  - `"out"` the exponent and the special values
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("mul", 3.0),
    StagePoint("out", 3.0))
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
  val stage: PipelineStageConfig,
) extends Module {

  val exW       = spec.exW
  val manW      = spec.manW
  val padding   = bitW.outBitW - manW
  val extraBits = padding

  val zCorrTmpW = scala.math.min(6, 1 + manW + padding)
  val coeffW    = scala.math.min(8, manW + 1)

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ExpPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info = new ExpPostProcInfo(spec, padding)
      val zres = UInt(bitW.outBitW.W)
    }))
    val z = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val in     = cutIn(io.in.bits)
  val infoIn = in.info
  val zres   = in.zres

  val zman0Comb = if (extraBits > 0) {
    (zres >> extraBits) +& zres(extraBits - 1)
  } else {
    zres
  }

  // ------------------------------------------------------------------------
  // stage "round": the rounded polynomial result

  private val cutRound = new SharedPipe("round", cutIn.valid, stages.at("round"))

  val zman0     = cutRound(zman0Comb)
  val infoRound = cutRound(infoIn)

  val zCorrectionComb = Wire(UInt(1.W))
  if (padding != 0) {
    // The correction asks whether (2^manW + zman0) * coef(k) reaches
    // 2^(manW + zCorrTmpW - 1), where k is the padding-wide xfracLSBs and
    // coef(k) = (ln2 * k) >> (coeffW + padding - zCorrTmpW).  k takes only
    // 2^padding values, so each one is a constant compare on zman0.
    val ln2       = (Real.log(Real.two))(coeffW).toBigInt
    val threshold = BigInt(1) << (manW + zCorrTmpW - 1)
    val hidden    = BigInt(1) << manW
    val corrCases = (0 until (1 << padding)).flatMap(k => {
      val coef = (ln2 * k) >> (coeffW + padding - zCorrTmpW)
      if (coef == 0) {
        None
      } else {
        val lo = (threshold + coef - 1) / coef - hidden
        if (lo >= hidden) {
          None
        } else if (lo <= 0) {
          Some(k.U(padding.W) -> true.B)
        } else {
          Some(k.U(padding.W) -> (zman0(manW - 1, 0) >= lo.U))
        }
      }
    })
    zCorrectionComb := MuxLookup(infoRound.xfracLSBs.get, false.B)(corrCases).asUInt
  } else {
    zCorrectionComb := 0.U(1.W)
  }

  // ------------------------------------------------------------------------
  // stage "corr": the correction term

  private val cutCorr = new SharedPipe("corr", cutRound.valid, stages.at("corr"))

  val zCorrection = cutCorr(zCorrectionComb)
  val zman0Corr   = cutCorr(zman0)
  val infoCorr    = cutCorr(infoRound)

  val zmanCorrectedComb = zman0Corr +& zCorrection
  val zexCorrection     = zmanCorrectedComb(manW) & !infoCorr.zIsNonTable
  val zExComb           = infoCorr.zex + zexCorrection

  // ------------------------------------------------------------------------
  // stage "add": the corrected mantissa and exponent

  private val cutAdd = new SharedPipe("add", cutCorr.valid, stages.at("add"))

  val zmanCorrected = cutAdd(zmanCorrectedComb)
  val zEx           = cutAdd(zExComb)
  val info          = cutAdd(infoCorr)

  // ------------------------------------------------------------------------
  // stage "out"

  val zMan = Mux(info.zIsNonTable, Cat(info.znan, 0.U((manW - 1).W)),
                 zmanCorrected(manW - 1, 0))
  val z = Cat(0.U(1.W), zEx, zMan)
  assert(z.getWidth == spec.W)

  io.z := Pipe(cutAdd.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ExpPostProcess]]. */
object ExpPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"round"` the rounded polynomial result
   *  - `"corr"`  the correction term
   *  - `"add"`   the corrected mantissa and exponent
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("round", 2.0),
    StagePoint("corr",  2.0),
    StagePoint("add",   1.0),
    StagePoint("out",   1.0))
}
