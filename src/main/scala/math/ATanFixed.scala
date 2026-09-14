//% @file ATanFixed.scala
//
// atan(t) for |t| <= 1, without either shared multiplier.
//
// The caller hands over the ratio itself: the min/max, the reciprocal and the
// quadrant of an atan2 stay on its side, so neither the pre-process multiplier
// nor a first phase takes part (see ATan2.scala for the pipeline that does
// them in hardware). What differs from that pipeline's second phase is what
// the table holds: atan(t) itself, not atan(t)/t, so the polynomial output is
// already the result as a fixed point value in [0, 1). Turning it into a float
// takes a leading-zero count and a shift instead of the post-process
// multiplier, which is what the modules below do.
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math.Pi

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{FloatChiselUtil, RealGeneric, RealSpec, RoundSpec}
import rial.util._

/** What [[rial.math.ATanFixedPreProcess]] leaves for
 *  [[rial.math.ATanFixedPostProcess]] beside the table address. */
class ATanFixedPostProcInfo(val spec: RealSpec) extends Bundle {
  val tsgn       = UInt(1.W)
  val tex        = UInt(spec.exW.W)
  val tman       = UInt(spec.manW.W)
  val isLinear   = Bool() // atan(t) is t to the working precision
  val znan       = Bool() // t is NaN
  val zzero      = Bool() // t is zero
  val zquarterPi = Bool() // |t| >= 1, where the result is atan(1) = pi/4
}

// -------------------------------------------------------------------------
//                                                          _
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/___/
// |_|            |_|
// -------------------------------------------------------------------------

/** Aligns |t| at the binary point, where it addresses the table, and decides
 *  the cases the table does not answer. */
class ATanFixedPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW  = bitW.addrW
  val dxW   = bitW.inBitW
  val order = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ATanFixedPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val adr  = Valid(UInt(adrW.W))
    val dx   = if (order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
    val info = Valid(new ATanFixedPostProcInfo(spec))
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val t = cutIn(io.x.bits)

  // atan(t) is t to the working precision below this exponent; the accurate
  // second phase of atan2 shortens the same range the same way.
  val linearThreshold = ATan2Phase2SimGeneric.calcLinearThreshold(manW) + exBias

  val tAtLeastOne = t.ex >= exBias.U(exW.W)

  // the four cases the polynomial does not answer, one of them at a time: NaN,
  // |t| >= 1 (where the result is atan(1)), zero, and the linear range.
  val infoComb = Wire(new ATanFixedPostProcInfo(spec))
  infoComb.tsgn       := t.sgn
  infoComb.tex        := t.ex
  infoComb.tman       := t.man
  infoComb.znan       := t.nan
  infoComb.zquarterPi := tAtLeastOne && !t.nan
  infoComb.zzero      := t.zero
  infoComb.isLinear   := !t.zero && (t.ex <= linearThreshold.U(exW.W))

  // ---- |t| as a fixed point value in [0, 1) with manW fraction bits ----
  // The shift saturates so that a tiny t leaves the address at zero and an
  // exponent at or above the bias, whose difference wraps, does not reach it.
  val shiftW    = log2Up(2 + manW)
  val tShift0   = exBias.U(exW.W) - t.ex
  val tShiftOut = tShift0(exW - 1, shiftW).orR
  val tShift    = Mux(tShiftOut, (manW + 1).U(shiftW.W),
                      Mux(tShift0(shiftW - 1, 0) > (manW + 1).U(shiftW.W),
                          (manW + 1).U(shiftW.W),
                          tShift0(shiftW - 1, 0)))

  val tFixedComb = (Cat(1.U(1.W), t.man) >> tShift)(manW, 0)

  // --------------------------------------------------------------------------
  // stage "shift": the argument aligned at the binary point

  private val cutShift = new SharedPipe("shift", cutIn.valid, stages.at("shift"))

  val tFixed = cutShift(tFixedComb)
  val info   = cutShift(infoComb)

  // --------------------------------------------------------------------------
  // stage "out"

  val adr = TableAddress.adr(tFixed(manW - 1, 0), adrW)
  io.adr := Pipe(cutShift.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = TableAddress.dx(tFixed(manW - 1, 0), adrW, dxW)
    io.dx.get := Pipe(cutShift.valid, dx, stages.at("out"))
  }

  io.info := Pipe(cutShift.valid, info, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ATanFixedPreProcess]]. */
object ATanFixedPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"shift"` the argument aligned at the binary point
   *  - `"out"`   the table address and the info
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("shift", 3.0),
    StagePoint("out",   1.0))
}

// -------------------------------------------------------------------------
//   table coefficient ROM
// -------------------------------------------------------------------------

class ATanFixedTableCoeff(
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
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

/** Turns [[ATanFixedPreProcess]]'s info and the polynomial output, a fixed
 *  point value in [0, 1), into the final atan(t) bit pattern. */
class ATanFixedPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val fracW = bitW.outBitW

  /** bits of the normalized polynomial output below the mantissa. */
  val extraBits = fracW - 1 - manW

  require(extraBits >= 0,
    f"outBitW(${fracW}) must hold a hidden one above a manW(${manW}) mantissa")

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ATanFixedPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val pre  = Flipped(Valid(new ATanFixedPostProcInfo(spec)))
    val zres = Flipped(Valid(UInt(fracW.W)))
    val z    = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.pre.valid && io.zres.valid,
    stages.at("in"))

  val info = cutIn(io.pre.bits)
  val zres = cutIn(io.zres.bits)

  val countComb = PriorityEncoder(Reverse(zres))

  // --------------------------------------------------------------------------
  // stage "lead": the position of the leading one

  private val cutLead = new SharedPipe("lead", cutIn.valid, stages.at("lead"))

  val count    = cutLead(countComb)
  val zresLead = cutLead(zres)
  val infoLead = cutLead(info)

  // a value in [0, 1) with `count` leading zeros lies in
  // [2^-(count+1), 2^-count), so its exponent is -(count + 1).
  val zshiftedComb = (zresLead << count)(fracW - 1, 0)
  val zex0Comb     = (exBias - 1).U(exW.W) - count

  // --------------------------------------------------------------------------
  // stage "shift": the normalized polynomial output

  private val cutShift = new SharedPipe("shift", cutLead.valid, stages.at("shift"))

  val zshifted  = cutShift(zshiftedComb)
  val zex0      = cutShift(zex0Comb)
  val infoShift = cutShift(infoLead)

  val zmanRoundedComb = WireDefault(0.U(manW.W))
  val zexRoundedComb  = WireDefault(0.U(exW.W))
  if (extraBits == 0) {
    zmanRoundedComb := zshifted(fracW - 2, 0)
    zexRoundedComb  := zex0
  } else {
    val lsb    = zshifted(extraBits)
    val round  = zshifted(extraBits - 1)
    val sticky = if (extraBits == 1) {
      false.B
    } else {
      zshifted(extraBits - 2, 0).orR
    }
    val roundInc = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
      lsb, round, sticky)
    // the carry of the rounding leaves a zero mantissa one exponent higher
    val zman0 = zshifted(fracW - 2, extraBits) +& roundInc
    zmanRoundedComb := zman0(manW - 1, 0)
    zexRoundedComb  := zex0 + zman0(manW)
  }

  // --------------------------------------------------------------------------
  // stage "round": the rounded mantissa and the exponent it leaves

  private val cutRound = new SharedPipe("round", cutShift.valid, stages.at("round"))

  val zmanRounded = cutRound(zmanRoundedComb)
  val zexRounded  = cutRound(zexRoundedComb)
  val infoRound   = cutRound(infoShift)

  // --------------------------------------------------------------------------
  // stage "out"

  val znan       = infoRound.znan
  val zquarterPi = infoRound.zquarterPi
  val zzero      = infoRound.zzero
  val isLinear   = infoRound.isLinear

  val quarterPi = new RealGeneric(spec, Pi * 0.25)

  val zex = Mux(znan,       Fill(exW, 1.U(1.W)),
            Mux(zquarterPi, quarterPi.ex.U(exW.W),
            Mux(zzero,      0.U(exW.W),
            Mux(isLinear,   infoRound.tex,
                            zexRounded))))
  val zman = Mux(znan,       Cat(1.U(1.W), 0.U((manW - 1).W)),
             Mux(zquarterPi, quarterPi.man.toBigInt.U(manW.W),
             Mux(zzero,      0.U(manW.W),
             Mux(isLinear,   infoRound.tman,
                             zmanRounded))))

  // a NaN carries no sign, as [[rial.arith.RealGeneric.nan]] makes it
  val zsgn = Mux(znan, 0.U(1.W), infoRound.tsgn)

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := Pipe(cutRound.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ATanFixedPostProcess]]. */
object ATanFixedPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"lead"`  the position of the leading one of the polynomial output
   *  - `"shift"` the normalized polynomial output and its exponent
   *  - `"round"` the rounded mantissa
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("lead",  1.0),
    StagePoint("shift", 1.0),
    StagePoint("round", 1.0),
    StagePoint("out",   1.0))
}
