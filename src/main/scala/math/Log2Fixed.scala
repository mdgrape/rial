//% @file Log2Fixed.scala
//
// log2(x) without either shared multiplier.
//
// log2 vanishes at x = 1, and a table indexed by the mantissa of x alone
// loses the arguments that approach it: what the result keeps there is the
// few leading bits of a fraction that carries a fixed number of them. What
// the table holds instead is |log2(x)| itself over the binades of the
// distance of x from one -- x - 1 on one side of it, 1 - x on the other,
// both exact in fixed point -- each binade scaled by the power of two that
// brings it into the same range as every other (see GradedTable.scala), so
// the result keeps its precision however close x comes to one. Turning the
// polynomial output into a float takes a leading-zero count and a shift
// instead of the post-process multiplier.
//
// Away from one the result is the unbiased exponent E of x plus log2(1 + t)
// of its mantissa t, and |E + log2(1 + t)| is at least one. The regions above
// x = 1 hold log2(1 + d) of the mantissa d, which is that same t, so the same
// table answers there too: its value shifted back by the scale of the binade
// it was read at is the fraction of the result, and the exponent carries the
// integer part. The modules below form the address of that table and assemble
// the result from what it leaves.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math.max

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{FloatChiselUtil, RealSpec, RoundSpec}
import rial.util._

/** What [[rial.math.Log2FixedPreProcess]] leaves for
 *  [[rial.math.Log2FixedPostProcess]] beside the table address. */
class Log2FixedPostProcInfo(
  val spec: RealSpec,
  val bitW: PolynomialBitWidth,
) extends Bundle {

  /** width of the binade the region was read at. */
  val binadeW = Log2FixedTable.binadeW(spec.manW)

  val zsgn  = UInt(1.W)
  val znan  = Bool()
  val zinf  = Bool()
  val zzero = Bool()

  /** x lies in (1/2, 2), where the table holds |log2 x| itself */
  val nearOne = Bool()

  /** x is a power of two, where the table has nothing to add */
  val whole = Bool()

  /** the magnitude of the integer part of the result */
  val zint = UInt(spec.exW.W)

  /** the binade the region was read at */
  val binade = UInt(binadeW.W)
}

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

/** Picks the region of the table that answers x -- the binade of the distance
 *  of x from one, on the side of one it lies on -- forms the address and the
 *  position inside its sub-interval, and decides the cases the table does not
 *  answer. */
class Log2FixedPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val dxW   = bitW.inBitW
  val order = bitW.nOrder

  /** The regions above x = 1 and the ones below it. */
  val above = Log2FixedTable.above(order, bitW.addrW, manW)
  val below = Log2FixedTable.below(order, bitW.addrW, manW)

  val subW    = above.subW
  val regionW = above.regionW

  /** Width of the whole address: the side of x = 1 above the region. */
  val adrW = 1 + regionW + subW

  /** width of the binade the region is looked up on. */
  val binadeW = Log2FixedTable.binadeW(manW)

  /** The most bits of the mantissa any binade of either side is split by. */
  val splitMaxW = max(above.splitMaxW, below.splitMaxW)

  require(above.regionW == below.regionW,
    f"the two sides of x = 1 share the region index of the address, so they " +
    f"must be the same width: ${above.regionW}%d above and " +
    f"${below.regionW}%d below")

  require(manW >= splitMaxW + subW,
    f"a mantissa of ${manW}%d bits must hold the " +
    f"${splitMaxW}%d bits that choose a region within a binade and " +
    f"the ${subW}%d that choose a sub-interval of it")

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = Log2FixedPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val adr  = Valid(UInt(adrW.W))
    val dx   = if (order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
    val info = Valid(new Log2FixedPostProcInfo(spec, bitW))
  })

  /** The binade a side looks its region up at: the nearest one that side
   *  holds a region for, as wide as its lookup is deep. */
  private def binadeIn(layout: GradedLayout, e: UInt): UInt = {
    val nearest = Mux(e > layout.maxBinade.U(binadeW.W),
                      layout.maxBinade.U(binadeW.W),
                  Mux(e < layout.minBinade.U(binadeW.W),
                      layout.minBinade.U(binadeW.W),
                      e))
    nearest(log2Ceil(layout.baseOf.length) - 1, 0)
  }

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  // ---- the distance of x from one ----
  // Below one the difference 1 - x is exact in fixed point (Sterbenz), and is
  // the complement of the mantissa there; above one, and away from one
  // altogether, the mantissa is the distance itself.
  val belowOneComb = x.ex === (exBias - 1).U(exW.W) && x.man.orR

  val diffComb = Mux(belowOneComb, (~x.man + 1.U)(manW - 1, 0), x.man)

  // --------------------------------------------------------------------------
  // stage "diff": the distance from one

  private val cutDiff = new SharedPipe("diff", cutIn.valid, stages.at("diff"))

  val diff  = cutDiff(diffComb)
  val xDiff = cutDiff(x)

  // a distance in [2^-e, 2^-(e-1)) has e - 1 leading zeros in a fraction of
  // manW bits, so the binade it lies in is that count plus one
  val leadComb = PriorityEncoder(Reverse(diff))

  // --------------------------------------------------------------------------
  // stage "lead": the leading zeros of the distance

  private val cutLead = new SharedPipe("lead", cutDiff.valid, stages.at("lead"))

  val lead     = cutLead(leadComb)
  val diffLead = cutLead(diff)
  val xLead    = cutLead(xDiff)

  val belowOne = xLead.ex === (exBias - 1).U(exW.W) && xLead.man.orR

  // below one the distance is that fraction halved, so its binade is one more
  val binadeComb = (lead +& 1.U + belowOne.asUInt)(binadeW - 1, 0)

  // the mantissa of the distance, its leading one dropped
  val manComb = (diffLead << (lead +& 1.U))(manW - 1, 0)

  // both sides are looked up and the one x lies on is chosen, so each of them
  // takes the binade at the nearest one it holds a region for
  val binadeAbove = binadeIn(above, binadeComb)
  val binadeBelow = binadeIn(below, binadeComb)

  val baseOfAbove = VecInit(above.baseOf.map(b => b.U(regionW.W)))
  val baseOfBelow = VecInit(below.baseOf.map(b => b.U(regionW.W)))

  val base = Mux(belowOne, baseOfBelow(binadeBelow), baseOfAbove(binadeAbove))

  // the region the binade holds, and the fraction the sub-interval index and
  // the position inside it are taken from
  val (regionComb, keyComb) = if (splitMaxW == 0) {
    (base, manComb)
  } else {
    // how many top bits of the mantissa choose the region within the binade
    // is itself a lookup on the binade, so it is this wide
    val splitWW      = log2Ceil(splitMaxW + 1)
    val splitOfAbove = VecInit(above.splitWOf.map(k => k.U(splitWW.W)))
    val splitOfBelow = VecInit(below.splitWOf.map(k => k.U(splitWW.W)))
    val splitW = Mux(belowOne, splitOfBelow(binadeBelow),
                               splitOfAbove(binadeAbove))
    val top    = manComb(manW - 1, manW - splitMaxW) >>
                 (splitMaxW.U(splitWW.W) - splitW)
    ((base +& top)(regionW - 1, 0), (manComb << splitW)(manW - 1, 0))
  }

  // --------------------------------------------------------------------------
  // stage "shift": the region of the table and the fraction inside it

  private val cutShift = new SharedPipe("shift", cutLead.valid, stages.at("shift"))

  val region = cutShift(regionComb)
  val key    = cutShift(keyComb)
  val side   = cutShift(belowOne)
  val binade = cutShift(binadeComb)
  val xShift = cutShift(xLead)

  // --------------------------------------------------------------------------
  // stage "out"

  val adr = Cat(side, region, TableAddress.adr(key, subW))
  assert(adr.getWidth == adrW)

  io.adr := Pipe(cutShift.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = TableAddress.dx(key, subW, dxW)
    io.dx.get := Pipe(cutShift.valid, dx, stages.at("out"))
  }

  // ---- special values ----
  //   log2(nan) -> nan      log2(inf) -> inf
  //   log2(0)   -> -inf     log2(neg) -> nan      log2(1) -> +0
  //
  // a subnormal argument is a zero here, as it is everywhere else, so it
  // gives -inf whatever its sign is.
  //
  // ---- integer part of the result ----
  //  x >= 1 (E >= 0):  log2(x) =    E      + frac
  //  x <  1 (E <  0):  log2(x) = -((|E|-1) + (1 - frac))
  val zintAbove = xShift.ex - exBias.U(exW.W)        // E,       for x >= 1
  val zintBelow = (exBias - 1).U(exW.W) - xShift.ex  // |E| - 1, for x <  1

  val info = Wire(new Log2FixedPostProcInfo(spec, bitW))
  info.zsgn    := xShift.ex < exBias.U
  info.znan    := xShift.nan || (xShift.sgn === 1.U && !xShift.zero)
  info.zinf    := xShift.inf || xShift.zero
  info.zzero   := xShift.ex === exBias.U && !xShift.man.orR
  info.nearOne := (xShift.ex === exBias.U ||
                   xShift.ex === (exBias - 1).U(exW.W)) && xShift.man.orR
  info.whole   := !xShift.man.orR
  info.zint    := Mux(xShift.ex >= exBias.U, zintAbove, zintBelow)
  info.binade  := binade

  io.info := Pipe(cutShift.valid, info, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.Log2FixedPreProcess]]. */
object Log2FixedPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"diff"`  the distance of the argument from one
   *  - `"lead"`  the leading zeros of the distance
   *  - `"shift"` the region of the table and the fraction inside it
   *  - `"out"`   the table address, the position inside its sub-interval and
   *              the info
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("diff",  2.0),
    StagePoint("lead",  1.0),
    StagePoint("shift", 2.0),
    StagePoint("out",   1.0))
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

/** Turns [[Log2FixedPreProcess]]'s info and the polynomial output into the
 *  final log2(x) bit pattern.
 *
 *  Near one the output is |log2(x)| scaled into [1/4, 1), which a
 *  leading-zero count and a shift normalise on their own. Away from one it is
 *  log2(1 + t) at the scale of the binade of the mantissa t, and the integer
 *  part of the result carries it: the two together are a fixed point
 *  magnitude of `exW + outBitW` bits. */
class Log2FixedPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val fracW = bitW.outBitW

  /** width of the fixed point magnitude of the result. */
  val magW = exW + fracW

  /** bits of the normalized magnitude below the mantissa. */
  val extraBits = magW - 1 - manW

  require(extraBits >= 1,
    f"exW(${exW}) + outBitW(${fracW}) must hold a hidden one and a round bit " +
    f"above a manW(${manW}) mantissa")

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = Log2FixedPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val pre  = Flipped(Valid(new Log2FixedPostProcInfo(spec, bitW)))
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

  // near one the region of the binade e holds |log2 x| scaled up by 2^(e-2),
  // and its value with `lead` leading zeros is a fraction in [1/2, 1) times
  // 2^-(lead + 1), so the result is that fraction times 2^(2 - e - lead); it
  // is normalized already, and sits at the top of the magnitude
  val znearComb   = Cat((zres << leadComb)(fracW - 1, 0), 0.U(exW.W))
  val zexNearComb = (exBias + 1).U(exW.W) - info.binade - leadComb

  // away from one the same region holds log2(1 + t) of the mantissa t scaled
  // up by 2^(e-2), so the fraction of the result is it shifted back. A power
  // of two leaves the table nothing to add.
  val fracComb = Mux(info.whole, 0.U(fracW.W),
                     (Cat(zres, 0.U(2.W)) >> info.binade)(fracW - 1, 0))

  // --------------------------------------------------------------------------
  // stage "frac": the fraction of the result

  private val cutFrac = new SharedPipe("frac", cutIn.valid, stages.at("frac"))

  val frac     = cutFrac(fracComb)
  val znear    = cutFrac(znearComb)
  val zexNear  = cutFrac(zexNearComb)
  val infoFrac = cutFrac(info)

  // Below one the magnitude is (|E| - 1) + (1 - frac). A frac of zero leaves
  // 1 - frac at one, an increment of the integer part, and the complement
  // below at zero.
  val belowOne = infoFrac.zsgn === 1.U
  val fracMag  = Mux(belowOne, (~frac + 1.U)(fracW - 1, 0), frac)
  val carry    = belowOne && !frac.orR

  val magComb = Cat(infoFrac.zint + carry.asUInt, fracMag)
  assert(magComb.getWidth == magW)

  val countComb = PriorityEncoder(Reverse(magComb))

  // --------------------------------------------------------------------------
  // stage "lead": the position of the leading one of the magnitude

  private val cutLead = new SharedPipe("lead", cutFrac.valid, stages.at("lead"))

  val count       = cutLead(countComb)
  val mag         = cutLead(magComb)
  val znearLead   = cutLead(znear)
  val zexNearLead = cutLead(zexNear)
  val infoLead    = cutLead(infoFrac)

  // the magnitude is an integer of `exW` bits and a fraction of `fracW` ones,
  // so with `count` leading zeros its leading one weighs 2^(exW - 1 - count).
  val zshiftedComb = Mux(infoLead.nearOne, znearLead,
                         (mag << count)(magW - 1, 0))
  val zex0Comb     = Mux(infoLead.nearOne, zexNearLead,
                         (exBias + exW - 1).U(exW.W) - count)

  // --------------------------------------------------------------------------
  // stage "shift": the normalized magnitude and its exponent

  private val cutShift = new SharedPipe("shift", cutLead.valid, stages.at("shift"))

  val zshifted  = cutShift(zshiftedComb)
  val zex0      = cutShift(zex0Comb)
  val infoShift = cutShift(infoLead)

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
  val zman0           = zshifted(magW - 2, extraBits) +& roundInc
  val zmanRoundedComb = zman0(manW - 1, 0)
  val zexRoundedComb  = zex0 + zman0(manW)

  // --------------------------------------------------------------------------
  // stage "round": the rounded mantissa and the exponent it leaves

  private val cutRound = new SharedPipe("round", cutShift.valid, stages.at("round"))

  val zmanRounded = cutRound(zmanRoundedComb)
  val zexRounded  = cutRound(zexRoundedComb)
  val infoRound   = cutRound(infoShift)

  // --------------------------------------------------------------------------
  // stage "out"

  val znan  = infoRound.znan
  val zinf  = infoRound.zinf
  val zzero = infoRound.zzero

  val zex = Mux(znan || zinf, Fill(exW, 1.U(1.W)),
            Mux(zzero,        0.U(exW.W),
                              zexRounded))
  val zman = Mux(znan,          Cat(1.U(1.W), 0.U((manW - 1).W)),
             Mux(zinf || zzero, 0.U(manW.W),
                                zmanRounded))
  // SW reference `RealGeneric.nan` uses sgn=0; force the sign bit to 0 when
  // we are emitting NaN so the test sees identical bit patterns.
  val zsgn = Mux(znan, 0.U(1.W), infoRound.zsgn)

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := Pipe(cutRound.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.Log2FixedPostProcess]]. */
object Log2FixedPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"frac"`  the fraction of the result
   *  - `"lead"`  the position of the leading one of the magnitude
   *  - `"shift"` the normalized magnitude and its exponent
   *  - `"round"` the rounded mantissa
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("frac",  1.0),
    StagePoint("lead",  1.0),
    StagePoint("shift", 1.0),
    StagePoint("round", 1.0),
    StagePoint("out",   1.0))
}
