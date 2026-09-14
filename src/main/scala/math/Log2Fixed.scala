//% @file Log2Fixed.scala
//
// log2(x) without either shared multiplier.
//
// log2(2^e * (1 + t)) = e + log2(1 + t), and the table holds log2(1 + t)
// itself, so the polynomial output is the fraction of the result as a fixed
// point value in [0, 1) and the integer part is the exponent of the argument.
// Turning the two into a float takes a leading-zero count and a shift instead
// of the post-process multiplier, which is what the modules below do.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{FloatChiselUtil, RealSpec, RoundSpec}
import rial.util._

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

/** Reads the table of log2(1 + t) at the mantissa of the argument. */
class Log2FixedPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val manW = spec.manW

  val adrW  = bitW.addrW
  val dxW   = bitW.inBitW
  val order = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = Log2FixedPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x   = Flipped(Valid(new DecomposedReal(spec)))
    val adr = Valid(UInt(adrW.W))
    val dx  = if (order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  // the mantissa is the position t in [0, 1) the table is read at
  val adr = TableAddress.adr(x.man, adrW)
  io.adr := Pipe(cutIn.valid, adr, stages.at("out"))

  if (order != 0) {
    // convert UInt 0000~1111 -> SInt 1000,..1111,0000,..0111
    val dx = TableAddress.dx(x.man, adrW, dxW)
    io.dx.get := Pipe(cutIn.valid, dx, stages.at("out"))
  }
}

/** Pipeline stage points of [[rial.math.Log2FixedPreProcess]]. */
object Log2FixedPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"out"` the table address
   */
  val stagePoints: Seq[StageElem] = Seq(StagePoint("out", 1.0))
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class Log2FixedTableCoeff(
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

  val packed = VecInit(PolynomialGenerator.packTable(coeff).map(c => {
    c.U(coefBW.W)
  }))
  io.cs := Pipe(io.adr.valid, packed(io.adr.bits), stage)
}

// -------------------------------------------------------------------------
//                        _        _     _                    _   _
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_|
//                                               |_|
// -------------------------------------------------------------------------

/** The integer part of log2(x) and its special values, which
 *  [[rial.math.Log2FixedNonTable]] forms beside the table. */
class Log2FixedPostProcInfo(val spec: RealSpec) extends Bundle {
  val zsgn  = UInt(1.W)
  val znan  = Bool()
  val zinf  = Bool()
  val zzero = Bool()
  /** the magnitude of the integer part of the result */
  val zint  = UInt(spec.exW.W)
}

// No pathway other than table interpolation: the integer part of the result
// is the exponent of the argument, and the sign tells the post-process which
// way to read the fraction the table gives.
class Log2FixedNonTable(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = Log2FixedNonTable.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x = Flipped(Valid(new DecomposedReal(spec)))
    val z = Valid(new Log2FixedPostProcInfo(spec))
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  // ---- special values ----
  //   log2(nan) -> nan      log2(inf) -> inf
  //   log2(0)   -> -inf     log2(neg) -> nan      log2(1) -> +0
  //
  // a subnormal argument is a zero here, as it is everywhere else, so it
  // gives -inf whatever its sign is.
  val znanComb  = x.nan || (x.sgn === 1.U && !x.zero)
  val zinfComb  = x.inf || x.zero
  val zzeroComb = x.ex === exBias.U && !x.man.orR
  val zsgnComb  = x.ex < exBias.U

  // ---- integer part of the result ----
  //  x >= 1 (e >= 0):  log2(x) =    e     + frac
  //  x <  1 (e <  0):  log2(x) = -((|e|-1) + (1 - frac))
  val zintAbove = x.ex - exBias.U(exW.W)        // e,       for x >= 1
  val zintBelow = (exBias - 1).U(exW.W) - x.ex  // |e| - 1, for x <  1
  val zintComb  = Mux(x.ex >= exBias.U, zintAbove, zintBelow)

  // --------------------------------------------------------------------------
  // stage "ex": the integer part of the result

  private val cutEx = new SharedPipe("ex", cutIn.valid, stages.at("ex"))

  val zint  = cutEx(zintComb)
  val zsgn  = cutEx(zsgnComb)
  val znan  = cutEx(znanComb)
  val zinf  = cutEx(zinfComb)
  val zzero = cutEx(zzeroComb)

  // --------------------------------------------------------------------------
  // stage "out"

  val out = Wire(new Log2FixedPostProcInfo(spec))
  out.zsgn  := zsgn
  out.znan  := znan
  out.zinf  := zinf
  out.zzero := zzero
  out.zint  := zint

  io.z := Pipe(cutEx.valid, out, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.Log2FixedNonTable]]. */
object Log2FixedNonTable {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"ex"`  the integer part of the result
   *  - `"out"` the integer part, the sign and the special values
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("ex",  2.0),
    StagePoint("out", 2.0))
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

/** Turns [[Log2FixedNonTable]]'s info and the polynomial output, a fixed point
 *  value in [0, 1), into the final log2(x) bit pattern.
 *
 *  The two form the result as a fixed point magnitude with `outBitW` fraction
 *  bits, so the precision of the assembly is absolute: the relative error
 *  grows as the result approaches zero. */
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
    val in = Flipped(Valid(new Bundle {
      val info = new Log2FixedPostProcInfo(spec)
      val zres = UInt(fracW.W)
    }))
    val z = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val in = cutIn(io.in.bits)

  val info = in.info
  val zres = in.zres

  // Below one the magnitude is (|e| - 1) + (1 - frac). A frac of zero leaves
  // 1 - frac at one, an increment of the integer part, and the complement
  // below at zero.
  val belowOne = info.zsgn === 1.U
  val fracMag  = Mux(belowOne, (~zres + 1.U)(fracW - 1, 0), zres)
  val carry    = belowOne && !zres.orR

  val magComb = Cat(info.zint + carry.asUInt, fracMag)
  assert(magComb.getWidth == magW)

  val countComb = PriorityEncoder(Reverse(magComb))

  // --------------------------------------------------------------------------
  // stage "lead": the position of the leading one of the magnitude

  private val cutLead = new SharedPipe("lead", cutIn.valid, stages.at("lead"))

  val count    = cutLead(countComb)
  val mag      = cutLead(magComb)
  val infoLead = cutLead(info)

  // the magnitude is an integer of `exW` bits and a fraction of `fracW` ones,
  // so with `count` leading zeros its leading one weighs 2^(exW - 1 - count).
  val zshiftedComb = (mag << count)(magW - 1, 0)
  val zex0Comb     = (exBias + exW - 1).U(exW.W) - count

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
   *  - `"lead"`  the position of the leading one of the magnitude
   *  - `"shift"` the normalized magnitude and its exponent
   *  - `"round"` the rounded mantissa
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("lead",  1.0),
    StagePoint("shift", 1.0),
    StagePoint("round", 1.0),
    StagePoint("out",   1.0))
}
