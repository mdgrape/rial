//% @file ATanFixed.scala
//
// atan(t) for |t| <= 1, without either shared multiplier.
//
// The caller hands over the ratio itself: the min/max, the reciprocal and the
// quadrant of an atan2 stay on its side, so neither the pre-process multiplier
// nor a first phase takes part (see ATan2.scala for the pipeline that does
// them in hardware). What differs from that pipeline's second phase is what
// the table holds: atan(t) itself, not atan(t)/t, so the polynomial output is
// already the result as a fixed point value. Turning it into a float takes a
// leading-zero count and a shift instead of the post-process multiplier.
//
// atan(t) vanishes with t, so the table gives a region of its own to each
// binade of |t| and stores its values scaled by a power of two (see
// GradedTable.scala): whatever the binade, the polynomial output fills the
// same range, and the result keeps its precision down to the smallest
// argument the table answers. Below that binade atan(t) is t itself, and at
// or above one, where the caller has left the domain, it is atan(1); neither
// needs the table. The modules below form the address of that table and
// assemble the result from what it leaves.
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

/** Picks the region of the table that answers |t|: the binade of |t|, halved
 *  by the top bit of its mantissa over the widest binades. It forms the
 *  address of that region and the position inside its sub-interval, and
 *  decides the cases the table does not answer. */
class ATanFixedPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val dxW   = bitW.inBitW
  val order = bitW.nOrder

  /** The regions of the table and the address that selects one. */
  val layout = ATanFixedTable.layout(order, bitW.addrW, manW)

  val subW    = layout.subW
  val regionW = layout.regionW
  val adrW    = layout.adrW

  /** width of the binade the region is looked up on. */
  val binadeW = log2Ceil(layout.maxBinade + 1)

  require(manW >= layout.splitMaxW + subW,
    f"a mantissa of ${manW}%d bits must hold the " +
    f"${layout.splitMaxW}%d bits that choose a region within a binade and " +
    f"the ${subW}%d that choose a sub-interval of it")

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ATanFixedPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val adr  = Valid(UInt(adrW.W))
    val dx   = if (order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
    val info = Valid(new ATanFixedPostProcInfo(spec))
  })

  val t = io.x.bits

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

  // |t| in [2^-e, 2^-(e-1)) has the binade e; a binade the table holds no
  // region of its own for reads the nearest one it does
  val binade0    = exBias.U(exW.W) - t.ex
  val binadeComb = Mux(binade0 > layout.maxBinade.U(exW.W),
                       layout.maxBinade.U(binadeW.W),
                   Mux(binade0 < layout.minBinade.U(exW.W),
                       layout.minBinade.U(binadeW.W),
                       binade0(binadeW - 1, 0)))

  val baseOf = VecInit(layout.baseOf.map(b => b.U(regionW.W)))

  // the region the binade holds, and the fraction the sub-interval index and
  // the position inside it are taken from
  val (regionComb, keyComb) = if (layout.splitMaxW == 0) {
    (baseOf(binadeComb), t.man)
  } else {
    // how many top bits of the mantissa choose the region is itself a
    // lookup on the binade, so it is this wide
    val splitWW = log2Ceil(layout.splitMaxW + 1)
    val splitOf = VecInit(layout.splitWOf.map(k => k.U(splitWW.W)))
    val splitW  = splitOf(binadeComb)
    val top     = t.man(manW - 1, manW - layout.splitMaxW) >>
                  (layout.splitMaxW.U(splitWW.W) - splitW)
    ((baseOf(binadeComb) +& top)(regionW - 1, 0),
     (t.man << splitW)(manW - 1, 0))
  }

  // --------------------------------------------------------------------------
  // stage "region": the region of the table and the fraction inside it

  private val cutRegion = new SharedPipe("region", io.x.valid,
    stages.at("region"))

  val region = cutRegion(regionComb)
  val key    = cutRegion(keyComb)
  val info   = cutRegion(infoComb)

  // --------------------------------------------------------------------------
  // stage "out"

  val adr = Cat(region, TableAddress.adr(key, subW))
  io.adr := Pipe(cutRegion.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = TableAddress.dx(key, subW, dxW)
    io.dx.get := Pipe(cutRegion.valid, dx, stages.at("out"))
  }

  io.info := Pipe(cutRegion.valid, info, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ATanFixedPreProcess]]. */
object ATanFixedPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"region"` the region of the table and the fraction inside it
   *  - `"out"`    the table address, the position inside its sub-interval and
   *               the info
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("region", 2.0),
    StagePoint("out",    1.0))
}

// -------------------------------------------------------------------------
//   table coefficient ROM
// -------------------------------------------------------------------------

class ATanFixedTableCoeff(
  val spec:  RealSpec,
  // the regions of ATanFixedTable, concatenated in address order
  val coeff: PolynomialCoefficients,
  val stage: Int,
) extends Module {

  val coefBW = coeff.bitW.coefBits.reduce(_+_)

  // the concatenated table is already in the order the address selects
  val all = PolynomialGenerator.packTable(coeff)

  val adrW = log2Ceil(all.length)

  require(all.length == (1 << adrW),
    s"expected a power of two of atan table entries, got ${all.length}")

  val io = IO(new Bundle {
    val adr = Flipped(Valid(UInt(adrW.W)))
    val cs  = Valid(UInt(coefBW.W))
  })

  val packed = VecInit(all.map(c => c.U(coefBW.W)))

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

/** Turns [[ATanFixedPreProcess]]'s info and the polynomial output, atan(t)
 *  of the argument scaled into [1/4, 1), into the final atan(t) bit
 *  pattern. */
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

  // every region holds at least a quarter, so its value has at most two
  // leading zeros: a leading-zero detector would be wasted here
  val leadComb = Mux(zres(fracW - 1), 0.U(2.W),
                 Mux(zres(fracW - 2), 1.U(2.W),
                                      2.U(2.W)))

  // --------------------------------------------------------------------------
  // stage "lead": the leading zeros of the value

  private val cutLead = new SharedPipe("lead", cutIn.valid, stages.at("lead"))

  val lead     = cutLead(leadComb)
  val zresLead = cutLead(zres)
  val infoLead = cutLead(info)

  // the region of the binade e holds atan(t) scaled up by 2^(e-1), and its
  // value with `lead` leading zeros is a fraction in [1/2, 1) times 2^-lead,
  // so the result is that fraction times 2^(1 - e - lead); e is exBias - tex
  val zshiftedComb = (zresLead << lead)(fracW - 1, 0)
  val zex0Comb     = infoLead.tex - lead

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
   *  - `"lead"`  the leading zeros of the polynomial output
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
