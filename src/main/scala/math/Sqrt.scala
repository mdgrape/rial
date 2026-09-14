//% @file sqrt.scala
//
// square root function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.RealSpec
import rial.util._

// An implementation of sqrt using (normally 2nd-order) polynomial interpolation
// Since FP is composed of exponent and mantissa, the result of sqrt will be
//     sqrt(2^ex * 1.m) = 2^(ex/2) * sqrt(1.m) .
// The fractional part of exponent appears when ex is an odd number, so
//     sqrt(2^ex * 1.m) = 2^(floor(ex/2)) * sqrt(1.m)   if ex is even
//                        2^(floor(ex/2)) * sqrt(2*1.m) otherwise.
// It requires 2 different tables for sqrt(1.m) and sqrt(2*1.m).
// To distinguish them, we use LSB of ex to detect whether ex is even or odd.

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

// sqrt(x): floating => floating
// - if x < 0, returns 0.
class SqrtPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val manW   = spec.manW

  val adrW   = bitW.addrW
  val fracW  = bitW.outBitW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = SqrtPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x   = Flipped(Valid(new DecomposedReal(spec)))
    val adr = Valid(UInt((1+adrW).W)) // we use LSB of x.ex
    val dx  = if(order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  // use the LSB of x.ex to distinguish 1~2 and 2~4
  val adr = Cat(x.ex(0), TableAddress.adr(x.man, adrW))
  io.adr := Pipe(cutIn.valid, adr, stages.at("out"))

  if(order != 0) {
    // convert UInt 0000~1111 -> SInt 1000,..1111,0000,..0111
    val dx = TableAddress.dx(x.man, adrW, dxW)
    io.dx.get := Pipe(cutIn.valid, dx, stages.at("out"))
  }
}

/** Pipeline stage points of [[rial.math.SqrtPreProcess]]. */
object SqrtPreProcess {

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

class SqrtTableCoeff(
  val spec:  RealSpec,
  val coeff: PolynomialCoefficients,
  val stage: Int,
) extends Module {

  val manW   = spec.manW
  val adrW   = coeff.bitW.addrW
  val fracW  = coeff.bitW.outBitW
  val order  = coeff.bitW.nOrder

  val io = IO(new Bundle {
    val adr = Flipped(Valid(UInt((1+adrW).W)))
    val cs  = Valid(UInt(coeff.bitW.coefBits.reduce(_+_).W))
  })

  val packed = VecInit(PolynomialGenerator.packTable(coeff).map(c => {
    c.U(coeff.bitW.coefBits.reduce(_+_).W)
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


// No pathway other than table interpolation. just calculate ex and sgn.
class SqrtNonTable(
  val spec:  RealSpec, // Input / Output floating spec
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = SqrtNonTable.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x = Flipped(Valid(new DecomposedReal(spec)))
    val z = Valid(new SqrtPostProcInfo(spec))
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  val xneg      = if (spec.disableSign) { false.B } else { x.sgn.asBool }
  val znanComb  = x.nan
  val zinfComb  = x.inf
  val zzeroComb = x.zero || xneg

  val xExNobias      = x.ex - exBias.U
  val zExShiftedComb = Cat(xExNobias(exW-1), xExNobias(exW-1, 1)) // arith right shift

  // --------------------------------------------------------------------------
  // stage "ex": the halved exponent

  private val cutEx = new SharedPipe("ex", cutIn.valid, stages.at("ex"))

  val zExShifted = cutEx(zExShiftedComb)
  val znan       = cutEx(znanComb)
  val zinf       = cutEx(zinfComb)
  val zzero      = cutEx(zzeroComb)

  // --------------------------------------------------------------------------
  // stage "out"

  val out = WireDefault(0.U.asTypeOf(new SqrtPostProcInfo(spec)))
  out.zex   := Mux(zinf || znan, Fill(exW, 1.U(1.W)),
               Mux(zzero,        0.U(exW.W),
                                 zExShifted + exBias.U))
  out.znan  := znan
  out.zinf  := zinf
  out.zzero := zzero

  io.z := Pipe(cutEx.valid, out, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.SqrtNonTable]]. */
object SqrtNonTable {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"ex"`  the halved exponent
   *  - `"out"` the exponent and the special values
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

class SqrtPostProcInfo(val spec: RealSpec) extends Bundle {
  val zex   = UInt(spec.exW.W)
  val znan  = Bool()
  val zinf  = Bool()
  val zzero = Bool()
}

class SqrtPostProcess(
  val spec:  RealSpec, // Input / Output floating spec
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = bitW.addrW
  val fracW  = bitW.outBitW
  val order  = bitW.nOrder

  val extraBits = fracW - manW

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = SqrtPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info = new SqrtPostProcInfo(spec)
      val zres = UInt(fracW.W)
    }))
    val z = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val in = cutIn(io.in.bits)

  val zmanRoundedComb = WireDefault(0.U(manW.W))
  if (extraBits == 0) {
    zmanRoundedComb := in.zres
  } else {
    val zman0 = (in.zres >> extraBits) +& in.zres(extraBits-1)
    val polynomialOvf = zman0(manW) === 1.U
    zmanRoundedComb := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW-1,0))
  }

  // --------------------------------------------------------------------------
  // stage "round": the rounded polynomial result

  private val cutRound = new SharedPipe("round", cutIn.valid, stages.at("round"))

  val zmanRounded = cutRound(zmanRoundedComb)
  val info        = cutRound(in.info)

  // --------------------------------------------------------------------------
  // stage "out"

  val zsgn  = 0.U
  val zex   = info.zex
  val znan  = info.znan
  val zinf  = info.zinf
  val zzero = info.zzero

  val zIsNonTable  = znan || zinf || zzero
  val zmanNonTable = Cat(znan, 0.U((manW-1).W))

  val zman = Mux(zIsNonTable, zmanNonTable, zmanRounded)
  val z = Cat(zsgn, zex, zman)

  io.z := Pipe(cutRound.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.SqrtPostProcess]]. */
object SqrtPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"round"` the rounded polynomial result
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("round", 2.0),
    StagePoint("out",   1.0))
}
