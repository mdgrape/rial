//% @file Reciprocal.scala
//
// Reciprocal function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.RealSpec
import rial.util._

// An implementation of reciprocal (1/x) using polynomial interpolation
// Since FP is composed of exponent and mantissa, the result of 1/x will be
//     1/x = 1/(2^ex * 1.m) = 2^(-ex-1) * 2/(1.m)

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

// reciprocal(x): floating => floating
class ReciprocalPreProcess(
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
  val stagePoints = ReciprocalPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x   = Flipped(Valid(new DecomposedReal(spec)))
    val adr = Valid(UInt(adrW.W))
    val dx  = if(order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  val adr = TableAddress.adr(x.man, adrW)
  io.adr := Pipe(cutIn.valid, adr, stages.at("out"))

  if(order != 0) {
    // convert UInt 0000~1111 -> SInt 1000,..1111,0000,..0111
    val dx = TableAddress.dx(x.man, adrW, dxW)
    io.dx.get := Pipe(cutIn.valid, dx, stages.at("out"))
  }
}

/** Pipeline stage points of [[rial.math.ReciprocalPreProcess]]. */
object ReciprocalPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"out"` the table address
   */
  val stagePoints: Seq[StageElem] = Seq(StagePoint("out", 1.0))
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
class ReciprocalNonTable(
  val spec:  RealSpec, // Input / Output floating spec
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = ReciprocalNonTable.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x = Flipped(Valid(new DecomposedReal(spec)))
    val z = Valid(new ReciprocalPostProcInfo(spec))
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  val znanComb  = x.nan
  val zinfComb  = x.zero
  val zzeroComb = x.inf

  val xmanNonZero = x.man.orR
  val zex0Comb    = (exBias<<1).U((exW+1).W) - x.ex - xmanNonZero.asUInt
  val zexMSBComb  = zex0Comb(exW)

  // Check overflow/underflow
  // 1/x = 2^(-e-1) * 2/1.m
  val invExMax = -spec.exMin - 1
  val invExMin = -spec.exMax - 1
  val (zinfOvfComb, zzeroUnfComb) = if (spec.exMax < invExMax) { // check overflow
    (zexMSBComb, zex0Comb === 0.U)
  } else if (invExMin < spec.exMin) { // check underflow
    (false.B, zexMSBComb || zex0Comb === 0.U)
  } else { // nothing happens.
    (false.B, zex0Comb === 0.U)
  }

  val zman0Comb = (x.man === 0.U) && (!x.zero)
  // if x.man == 0, that means x = 2^N. then z = 2^-N, so zman = 0.

  // --------------------------------------------------------------------------
  // stage "ex": the negated exponent and the over/underflow flags

  private val cutEx = new SharedPipe("ex", cutIn.valid, stages.at("ex"))

  val zex0     = cutEx(zex0Comb)
  val znan     = cutEx(znanComb)
  val zinf     = cutEx(zinfComb)
  val zzero    = cutEx(zzeroComb)
  val zinfOvf  = cutEx(zinfOvfComb)
  val zzeroUnf = cutEx(zzeroUnfComb)
  val zman0    = cutEx(zman0Comb)
  val xsgn     = cutEx(x.sgn)

  // --------------------------------------------------------------------------
  // stage "out"

  val out = WireDefault(0.U.asTypeOf(new ReciprocalPostProcInfo(spec)))
  out.zsgn  := Mux(znan || zzero || zinf || zinfOvf || zzeroUnf, 0.U, xsgn)
  out.zex   := Mux(znan || zinf || zinfOvf, Fill(exW, 1.U(1.W)),
               Mux(zzero || zzeroUnf,       0.U(exW.W),
                                            zex0(exW-1,0)))
  out.znan  := znan
  out.zinf  := zinf  || zinfOvf
  out.zzero := zzero || zzeroUnf
  out.zman0 := zman0

  io.z := Pipe(cutEx.valid, out, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ReciprocalNonTable]]. */
object ReciprocalNonTable {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"ex"`  the negated exponent and the over/underflow flags
   *  - `"out"` the exponent, the sign and the special values
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("ex",  3.0),
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

class ReciprocalPostProcInfo(val spec: RealSpec) extends Bundle {
  val zsgn        = UInt(1.W)
  val zex         = UInt(spec.exW.W)
  val zinf        = Bool()
  val znan        = Bool()
  val zzero       = Bool()
  val zman0       = Bool()
}

// It just rounds the polynomial result.
// sqrt, invsqrt, reciprocal requires the same postprocess.
//
class ReciprocalPostProcess(
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
  val stagePoints = ReciprocalPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info = new ReciprocalPostProcInfo(spec)
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
    val zres0 = (in.zres >> extraBits) +& in.zres(extraBits-1)
    val polynomialOvf = zres0(manW)
    zmanRoundedComb := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zres0(manW-1,0))
  }

  // --------------------------------------------------------------------------
  // stage "round": the rounded polynomial result

  private val cutRound = new SharedPipe("round", cutIn.valid, stages.at("round"))

  val zmanRounded = cutRound(zmanRoundedComb)
  val info        = cutRound(in.info)

  // --------------------------------------------------------------------------
  // stage "out"

  val zsgn  = info.zsgn
  val zex   = info.zex
  val znan  = info.znan
  val zinf  = info.zinf
  val zzero = info.zzero
  val zman0 = info.zman0

  val zIsNonTable  = zinf || znan || zzero || zman0
  val zmanNonTable = Cat(znan, 0.U((manW-1).W))

  val zman = Mux(zIsNonTable, zmanNonTable, zmanRounded)
  val z = Cat(zsgn, zex, zman)

  io.z := Pipe(cutRound.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.ReciprocalPostProcess]]. */
object ReciprocalPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"round"` the rounded polynomial result
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("round", 2.0),
    StagePoint("out",   1.0))
}

