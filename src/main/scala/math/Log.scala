//% @file Log.scala
//
// log(x) hardware modules.
//
// - path 0 (x < 0.5 or x >= 2)  ->  log2(1+t)
// - path 1 (0.5 <= x < 1.0)     -> -log((1+t)/2)/(1-(1+t)/2)-1
// - path 2 (1.0 <= x < 2.0)     ->  log(1+t)/t
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
import rial.arith.RealSpec
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class LogPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = bitW.addrW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = LogPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x   = Flipped(Valid(new DecomposedReal(spec)))
    val adr = Valid(UInt((2 + adrW).W)) // 2 MSBs = path selector (exAddr)
    val dx  = if (order != 0) Some(Valid(UInt(dxW.W))) else None
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  // exAddr = 0 : x < 0.5 or x >= 2        (normal/polynomial)
  // exAddr = 1 : 0.5 <= x < 1.0           (smallNeg)
  // exAddr = 2 : 1.0 <= x < 2.0           (smallPos)
  val exAddr = Mux(x.ex === exBias.U,       2.U(2.W),
               Mux(x.ex === (exBias - 1).U, 1.U(2.W),
                                            0.U(2.W)))
  val subAdr = TableAddress.adr(x.man, adrW)
  val adr    = Cat(exAddr, subAdr)
  io.adr := Pipe(cutIn.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = TableAddress.dx(x.man, adrW, dxW)
    io.dx.get := Pipe(cutIn.valid, dx, stages.at("out"))
  }
}

/** Pipeline stage points of [[rial.math.LogPreProcess]]. */
object LogPreProcess {

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

class LogTableCoeff(
  val spec:  RealSpec,
  // concat [N coeffs] ++ [G coeffs] ++ [P coeffs] ++ [zeros]
  val coeff: PolynomialCoefficients,
  val stage: Int,
) extends Module {

  val adrW   = coeff.bitW.addrW
  val coefBW = coeff.bitW.coefBits.reduce(_+_)

  val io = IO(new Bundle {
    val adr = Flipped(Valid(UInt((2 + adrW).W)))
    val cs  = Valid(UInt(coefBW.W))
  })

  // The concatenated table is already in the order the address selects.
  val all = PolynomialGenerator.packTable(coeff)
  require(all.length == 4 * (1 << adrW),
    s"expected 4 * 2^${adrW} log table entries, got ${all.length}")

  val packed = VecInit(all.map(c => c.U(coefBW.W)))

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

class LogNonTableInfo(val spec: RealSpec) extends Bundle {
  val zsgn       = UInt(1.W)
  val znan       = Bool()
  val zinf       = Bool()
  val zzero      = Bool()
  val xHalfToOne  = Bool()  // path 1
  val xOneToTwo  = Bool()  // path 2
  val xOtherwise = Bool()  // path 0
  val constant   = UInt((spec.manW + 1).W) // 1.constant in (1-x)/(x-1)/ln2 form
  val zex0       = UInt(spec.exW.W)        // result exponent before mul correction
  val zint       = UInt(spec.exW.W)        // |xexNobias|; used for path 0
  val zintShift  = UInt(log2Up(spec.exW).W) // leading-zero count of zint
}

/** The path selectors and special-value flags [[LogNonTable]] carries across
 *  its own pipeline boundaries. */
class LogNonTableFlags extends Bundle {
  val zsgn       = UInt(1.W)
  val znan       = Bool()
  val zinf       = Bool()
  val zzero      = Bool()
  val xHalfToOne  = Bool()
  val xOneToTwo  = Bool()
  val xOtherwise = Bool()
}

class LogNonTable(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = LogNonTable.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x = Flipped(Valid(new DecomposedReal(spec)))
    val z = Valid(new LogNonTableInfo(spec))
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  val xex  = x.ex
  val xman = x.man
  val xsgn = x.sgn

  val xHalfToOneComb  = xex === (exBias - 1).U
  val xOneToTwoComb  = xex === exBias.U
  val xOtherwiseComb = !xHalfToOneComb && !xOneToTwoComb

  // ---- the difference the constant is built from ----
  //  path 1 (0.5 <= x < 1.0):  c = 1 - x  = (2^manW - xman) / 2^(manW+1)
  val c0_5to1_0 = ~xman +& 1.U                  // width manW+1
  //  path 2 (1.0 <= x < 2.0):  c = x - 1  = xman / 2^manW
  val c1_0to2_0 = Cat(0.U(1.W), xman)           // width manW+1

  val c0Comb = Mux(xHalfToOneComb, c0_5to1_0, c1_0to2_0)

  // ---- integer part of z (for path 0) ----
  val xexPos   = xex - exBias.U
  val xexNeg   = (exBias - 1).U - xex
  val zintComb = Mux(xex >= exBias.U, xexPos, xexNeg)

  // ---- special values ----
  //   log(nan)  -> nan       log(inf)  -> inf
  //   log(0)    -> -inf      log(neg)  -> nan       log(1) -> 0
  val znanComb  = x.nan || (xsgn === 1.U && !x.zero)
  val zinfComb  = x.inf || x.zero
  val zzeroComb = !xman.orR && xex === exBias.U
  val zsgnComb  = xex < exBias.U

  val flagsComb = Wire(new LogNonTableFlags)
  flagsComb.zsgn       := zsgnComb
  flagsComb.znan       := znanComb
  flagsComb.zinf       := zinfComb
  flagsComb.zzero      := zzeroComb
  flagsComb.xHalfToOne  := xHalfToOneComb
  flagsComb.xOneToTwo  := xOneToTwoComb
  flagsComb.xOtherwise := xOtherwiseComb

  // --------------------------------------------------------------------------
  // stage "diff": the difference 1-x / x-1 and the integer part of the result

  private val cutDiff = new SharedPipe("diff", cutIn.valid, stages.at("diff"))

  val c0Diff    = cutDiff(c0Comb)
  val zintDiff  = cutDiff(zintComb)
  val flagsDiff = cutDiff(flagsComb)

  val c0ShiftComb   = PriorityEncoder(Reverse(c0Diff))
  val zintShiftComb = PriorityEncoder(Reverse(zintDiff))

  // --------------------------------------------------------------------------
  // stage "clz": the leading-zero counts of both

  private val cutClz = new SharedPipe("clz", cutDiff.valid, stages.at("clz"))

  val c0Shift   = cutClz(c0ShiftComb)
  val zintShift = cutClz(zintShiftComb)
  val c0        = cutClz(c0Diff)
  val zint      = cutClz(zintDiff)
  val flags     = cutClz(flagsDiff)

  // --------------------------------------------------------------------------
  // stage "out"

  val c0Shifted = (c0 << c0Shift)(manW, 0)
  assert(!cutClz.valid || (BigInt(1) << manW).U <= c0Shifted || c0 === 0.U)

  val zex0closeto1 = (exBias - 1).U(exW.W) - c0Shift

  // path 0 constant: ln(2)
  val ln2 = Real.log(Real.two)(manW + 1).toBigInt.U((manW + 1).W)
  assert((BigInt(1) << manW).U < ln2 && ln2 < (BigInt(1) << (manW + 1)).U)

  val constant = Mux(flags.xOtherwise, ln2, c0Shifted)

  val zex0farfrom1 = (exBias + exW - 2).U - zintShift
  val zex0 = Mux(flags.xOtherwise, zex0farfrom1, zex0closeto1)

  val zout = Wire(new LogNonTableInfo(spec))
  zout.zsgn       := flags.zsgn
  zout.znan       := flags.znan
  zout.zinf       := flags.zinf
  zout.zzero      := flags.zzero
  zout.xHalfToOne  := flags.xHalfToOne
  zout.xOneToTwo  := flags.xOneToTwo
  zout.xOtherwise := flags.xOtherwise
  zout.constant   := constant
  zout.zex0       := zex0
  zout.zint       := zint
  zout.zintShift  := zintShift

  io.z := Pipe(cutClz.valid, zout, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.LogNonTable]]. */
object LogNonTable {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"diff"` the difference 1-x / x-1 and the integer part of the result
   *  - `"clz"`  the leading-zero counts of both
   *  - `"out"`  the multiplier constant, the exponent and the special values
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("diff", 3.0),
    StagePoint("clz",  4.0),
    StagePoint("out",  4.0))
}

// -------------------------------------------------------------------------
// post-polynomial argument builder. Combinational: consumes (zres + non-
// table info) and emits the (lhs, rhs) the shared PostProcMultiplier needs.
// -------------------------------------------------------------------------

class LogPostMultArgsOutput(val spec: RealSpec, val fracW: Int) extends Bundle {
  val lhs = UInt((1 + fracW).W)
  val rhs = UInt((1 + spec.manW).W)
}

class LogPostMultArgs(
  val spec: RealSpec,
  val bitW: PolynomialBitWidth,
) extends Module {

  val exW   = spec.exW
  val manW  = spec.manW
  val fracW = bitW.outBitW

  val io = IO(new Bundle {
    val info = Flipped(Valid(new LogNonTableInfo(spec)))
    val zres = Input(UInt(fracW.W))
    val mul  = Valid(new LogPostMultArgsOutput(spec, fracW))
  })

  val info = io.info.bits
  val zres = io.zres

  // For path 0 (xOtherwise), build zfull0 = (zint << fracW) + signed-adjusted zres.
  // For path 1 / 2, polynomial output directly forms the mantissa.
  val zfrac0 = Mux(info.zsgn === 0.U, zres,
               Mux(zres === 0.U, Fill(fracW, 1.U(1.W)), ~zres + 1.U))
  val zfull0       = Cat(info.zint, zfrac0)                       // exW + fracW bits
  val zfullShifted = (zfull0 << info.zintShift)(exW + fracW - 1, 0)
  assert(!io.info.valid || zfullShifted(exW + fracW - 1) === 1.U || !info.xOtherwise)

  val ymanW1 = Mux(info.xHalfToOne, Cat(1.U(1.W), zres),                     // path 1
               Mux(info.xOneToTwo, Cat(zres, 0.U(1.W)),                     // path 2
                                   zfullShifted(exW + fracW - 1, exW - 1))) // path 0

  val out = Wire(new LogPostMultArgsOutput(spec, fracW))
  out.lhs := ymanW1
  out.rhs := Cat(1.U(1.W), info.constant)

  io.mul.valid := io.info.valid
  io.mul.bits  := out
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class LogPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = LogPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val info = Flipped(Valid(new LogNonTableInfo(spec)))
    val mul  = Flipped(Valid(new PostProcMultiplierOutput(manW)))
    val z    = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.info.valid && io.mul.valid,
    stages.at("in"))

  val infoIn = cutIn(io.info.bits)
  val mul    = cutIn(io.mul.bits)

  val zman0Comb = mul.zman
  val zex0Comb  = infoIn.zex0 +& mul.exInc

  // --------------------------------------------------------------------------
  // stage "ex": the result exponent

  private val cutEx = new SharedPipe("ex", cutIn.valid, stages.at("ex"))

  val zman0 = cutEx(zman0Comb)
  val zex0  = cutEx(zex0Comb)
  val info  = cutEx(infoIn)

  // --------------------------------------------------------------------------
  // stage "out"

  val znan  = info.znan
  val zinf  = info.zinf || zex0(exW) === 1.U // exponent overflowed
  val zzero = info.zzero
  val zsgn  = info.zsgn

  val zman = Mux(znan || zinf || zzero, Cat(znan.asUInt, Fill(manW - 1, 0.U(1.W))), zman0)
  val zex  = Mux(znan || zinf, Fill(exW, 1.U(1.W)),
             Mux(zzero,        0.U(exW.W),
                               zex0(exW - 1, 0)))
  // SW reference `RealGeneric.nan` uses sgn=0; force the sign bit to 0 when
  // we are emitting NaN so the test sees identical bit patterns.
  val zsgnOut = Mux(znan, 0.U(1.W), zsgn)

  val z = Cat(zsgnOut, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := Pipe(cutEx.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.LogPostProcess]]. */
object LogPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"ex"`  the result exponent
   *  - `"out"` the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("ex",  1.0),
    StagePoint("out", 1.0))
}
