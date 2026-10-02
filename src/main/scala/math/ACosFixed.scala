//% @file ACosFixed.scala
//
// acos(x) in one pass, without either shared multiplier.
//
// The accurate acos takes two phases and a multiplier: the first maps x to
// y = sqrt(1 - |x|) through the Sqrt table, the second reads acos(1 - y^2)/y
// and multiplies it back by y (see ACos.scala). What the modules below do
// instead is to read the angle itself out of a table indexed by the binade of
// u = 1 - |x|, which is exact in fixed point for the |x| in [1/2, 1) that has
// one, and by |x| itself below a half, where acos is smooth. Each binade
// stores its angle scaled by a power of two, so the polynomial output is the
// angle as a fixed point value in [1/4, 1) whatever the binade: turning it
// into a float takes a leading-zero count and a shift instead of the
// post-process multiplier, and the angle of a negative x, pi - theta, is a
// subtraction at that same fixed point scale.
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math._

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{FloatChiselUtil, RealGeneric, RealSpec, RoundSpec}
import rial.util._

/** What [[rial.math.ACosFixedPreProcess]] leaves for
 *  [[rial.math.ACosFixedPostProcess]] beside the table address. */
class ACosFixedPostProcInfo(
  val spec: RealSpec,
  val bitW: PolynomialBitWidth,
) extends Bundle {

  /** width of the region index the table address carries. */
  val regionW = ACosFixedTable.regionW(bitW.outBitW)

  val xsgn    = UInt(1.W)
  val region  = UInt(regionW.W) // 0 for |x| < 1/2, else the binade of 1 - |x|
  val znan    = Bool() // x is NaN or an infinity
  val zzero   = Bool() // x is +1 or beyond, where acos is +0
  val zpi     = Bool() // x is -1 or beyond, where acos is pi
  val zhalfPi = Bool() // x is zero, where acos is pi/2
}

// -------------------------------------------------------------------------
//                                                          _
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/___/
// |_|            |_|
// -------------------------------------------------------------------------

/** Picks the region of the table that answers x -- the argument itself below
 *  a half, the binade of 1 - |x| above it -- forms the address and the
 *  position inside its sub-interval, and decides the cases the table does not
 *  answer. */
class ACosFixedPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val dxW   = bitW.inBitW
  val order = bitW.nOrder

  /** sub-intervals a region is divided into, as a power of two. */
  val subW = ACosFixedTable.subW(order, bitW.addrW)

  /** width of the region index the address carries above the sub-interval. */
  val regionW = ACosFixedTable.regionW(bitW.outBitW)

  val adrW = regionW + subW

  require(manW >= subW,
    f"a mantissa of ${manW}%d bits must hold a sub-interval index of " +
    f"${subW}%d bits")

  require(ACosFixedTable.nRegion(bitW.outBitW) >= manW + 2,
    f"outBitW(${bitW.outBitW}%d) must leave the table a region for each of " +
    f"the ${manW + 1}%d binades a difference of manW(${manW}%d) bits has")

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ACosFixedPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val adr  = Valid(UInt(adrW.W))
    val dx   = if (order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
    val info = Valid(new ACosFixedPostProcInfo(spec, bitW))
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  // ---- |x| as a fraction of manW + 1 bits ----
  // The shift saturates so that a tiny x leaves the fraction at zero, where
  // the table gives acos(0); |x| in [1/2, 1) is not shifted at all.
  val shiftW      = log2Up(2 + manW)
  val alignBy0    = (exBias - 1).U(exW.W) - x.ex
  val alignTooFar = alignBy0(exW - 1, shiftW).orR
  val alignBy     = Mux(alignTooFar, (manW + 1).U(shiftW.W),
                    Mux(alignBy0(shiftW - 1, 0) > (manW + 1).U(shiftW.W),
                        (manW + 1).U(shiftW.W),
                        alignBy0(shiftW - 1, 0)))

  val xFixedComb = (Cat(1.U(1.W), x.man) >> alignBy)(manW, 0)

  // ---- u = 1 - |x| ----
  // Sterbenz: for the |x| in [1/2, 1) whose binade is read below, the
  // difference of the two fractions is exact.
  val uComb = ((BigInt(1) << (manW + 1)).U((manW + 2).W) -
               Cat(1.U(1.W), x.man))(manW, 0)

  // --------------------------------------------------------------------------
  // stage "diff": the difference 1 - |x|

  private val cutDiff = new SharedPipe("diff", cutIn.valid, stages.at("diff"))

  val xFixedDiff = cutDiff(xFixedComb)
  val uDiff      = cutDiff(uComb)
  val xDiff      = cutDiff(x)

  // u in [2^-e, 2^-(e-1)) has e - 1 leading zeros in a fraction of manW + 1
  // bits, so the binade it lies in is that count plus one
  val leadComb = PriorityEncoder(Reverse(uDiff))

  // --------------------------------------------------------------------------
  // stage "lead": the leading zeros of the difference

  private val cutLead = new SharedPipe("lead", cutDiff.valid, stages.at("lead"))

  val lead       = cutLead(leadComb)
  val uLead      = cutLead(uDiff)
  val xFixedLead = cutLead(xFixedDiff)
  val xLead      = cutLead(xDiff)

  val xBelowHalf = xLead.ex < (exBias - 1).U(exW.W)

  // the mantissa of u, its leading one dropped
  val uManComb = (uLead << lead)(manW - 1, 0)

  // the fraction the table is read at, and the region that holds it
  val keyComb    = Mux(xBelowHalf, xFixedLead(manW - 1, 0), uManComb)
  val regionComb = Mux(xBelowHalf, 0.U(regionW.W), (lead +& 1.U)(regionW - 1, 0))

  // --------------------------------------------------------------------------
  // stage "shift": the fraction the table is read at, and its region

  private val cutShift = new SharedPipe("shift", cutLead.valid, stages.at("shift"))

  val key    = cutShift(keyComb)
  val region = cutShift(regionComb)
  val xShift = cutShift(xLead)

  // --------------------------------------------------------------------------
  // stage "out"

  val adr = Cat(region, TableAddress.adr(key, subW))
  io.adr := Pipe(cutShift.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = TableAddress.dx(key, subW, dxW)
    io.dx.get := Pipe(cutShift.valid, dx, stages.at("out"))
  }

  // ---- the cases the table does not answer ----
  // acos is defined on [-1, 1] alone; a finite argument beyond an end gets
  // the value of that end, so one that overshoots it by rounding needs no
  // clamp. An infinity and a NaN give NaN.
  val xAtLeastOne = xShift.ex >= exBias.U(exW.W) && !xShift.inf && !xShift.nan

  val info = Wire(new ACosFixedPostProcInfo(spec, bitW))
  info.xsgn    := xShift.sgn
  info.region  := region
  info.znan    := xShift.nan || xShift.inf
  info.zzero   := xAtLeastOne && xShift.sgn === 0.U
  info.zpi     := xAtLeastOne && xShift.sgn === 1.U
  info.zhalfPi := xShift.zero
  io.info := Pipe(cutShift.valid, info, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ACosFixedPreProcess]]. */
object ACosFixedPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"diff"`  the difference `1-|x|`
   *  - `"lead"`  the leading zeros of the difference
   *  - `"shift"` the fraction the table is read at, and its region
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

/** Turns [[ACosFixedPreProcess]]'s info and the polynomial output, the angle
 *  of |x| scaled into [1/4, 1), into the final acos(x) bit pattern. */
class ACosFixedPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val fracW = bitW.outBitW

  /** width of the region index of the info. */
  val regionW = ACosFixedTable.regionW(fracW)

  /** guard bits below the polynomial output the subtraction works at. */
  val guardW = ACosFixedTable.guardW

  /** width of an angle in [0, pi] in units of 2^-(fracW + guardW). */
  val angleW = 2 + fracW + guardW

  /** bits of the normalized angle below the mantissa. */
  val extraBits = angleW - 1 - manW

  require(extraBits >= 1,
    f"outBitW(${fracW}) must hold a hidden one and a round bit above a " +
    f"manW(${manW}) mantissa")

  val realPi     = new RealGeneric(spec, Pi)
  val realHalfPi = new RealGeneric(spec, Pi * 0.5)

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ACosFixedPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val pre  = Flipped(Valid(new ACosFixedPostProcInfo(spec, bitW)))
    val zres = Flipped(Valid(UInt(fracW.W)))
    val z    = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.pre.valid && io.zres.valid,
    stages.at("in"))

  val info = cutIn(io.pre.bits)
  val zres = cutIn(io.zres.bits)

  val pi = ACosFixedTable.pi(fracW + guardW).toBigInt.U(angleW.W)

  // every region holds at least a quarter, so its value has at most two
  // leading zeros: a leading-zero detector would be wasted here
  val leadComb = Mux(zres(fracW - 1), 0.U(2.W),
                 Mux(zres(fracW - 2), 1.U(2.W),
                                      2.U(2.W)))

  // the power of two the region stored its angle scaled up by
  val thetaShift = Mux(info.region === 0.U, 1.U(regionW.W),
                       info.region(regionW - 1, 1))

  // the angle of |x| at the scale the subtraction below works at
  val thetaComb = (Cat(zres, 0.U((2 + guardW).W)) >> thetaShift)(angleW - 1, 0)

  // --------------------------------------------------------------------------
  // stage "lead": the leading zeros of the value, and the angle of |x|

  private val cutLead = new SharedPipe("lead", cutIn.valid, stages.at("lead"))

  val lead     = cutLead(leadComb)
  val theta    = cutLead(thetaComb)
  val zresLead = cutLead(zres)
  val infoLead = cutLead(info)

  // the angle of a negative x is pi minus the one of |x|, exact at this scale
  val angleNeg        = pi - theta
  val angleAtLeastTwo = angleNeg(angleW - 1)
  val angleNegShifted = (angleNeg << Mux(angleAtLeastTwo, 0.U, 1.U))(angleW - 1, 0)

  // the value of a region with `lead` leading zeros is a fraction in [1/2, 1)
  // times 2^-(lead + 1), and the angle of |x| that fraction times 2 (region
  // zero) or 4 / 2^floor(e/2) (the binade e)
  val zexPos = Mux(infoLead.region === 0.U,
                   exBias.U(exW.W) - lead,
                   (exBias + 1).U(exW.W) - lead - infoLead.region(regionW - 1, 1))

  val zshiftedComb = Mux(infoLead.xsgn === 0.U,
                         (zresLead << (lead +& (2 + guardW).U))(angleW - 1, 0),
                         angleNegShifted)
  val zex0Comb     = Mux(infoLead.xsgn === 0.U, zexPos,
                         exBias.U(exW.W) + angleAtLeastTwo)

  // --------------------------------------------------------------------------
  // stage "shift": the normalized angle and its exponent

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
  val zman0           = zshifted(angleW - 2, extraBits) +& roundInc
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

  val znan    = infoRound.znan
  val zzero   = infoRound.zzero
  val zpi     = infoRound.zpi
  val zhalfPi = infoRound.zhalfPi

  val zex = Mux(znan,    Fill(exW, 1.U(1.W)),
            Mux(zzero,   0.U(exW.W),
            Mux(zpi,     realPi.ex.U(exW.W),
            Mux(zhalfPi, realHalfPi.ex.U(exW.W),
                         zexRounded))))
  val zman = Mux(znan,    Cat(1.U(1.W), 0.U((manW - 1).W)),
             Mux(zzero,   0.U(manW.W),
             Mux(zpi,     realPi.man.toBigInt.U(manW.W),
             Mux(zhalfPi, realHalfPi.man.toBigInt.U(manW.W),
                          zmanRounded))))

  // acos is an angle in [0, pi]: it carries no sign
  val z = Cat(0.U(1.W), zex, zman)
  assert(z.getWidth == spec.W)

  io.z := Pipe(cutRound.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ACosFixedPostProcess]]. */
object ACosFixedPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"lead"`  the leading zeros of the polynomial output and the angle of
   *              `|x|`
   *  - `"shift"` the normalized angle and its exponent
   *  - `"round"` the rounded mantissa
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("lead",  1.0),
    StagePoint("shift", 1.0),
    StagePoint("round", 1.0),
    StagePoint("out",   1.0))
}
