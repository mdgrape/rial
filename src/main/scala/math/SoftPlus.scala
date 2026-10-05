//% @file SoftPlus.scala
//
// SoftPlus (log(1 + exp(x))) hardware modules.
//
//   pre  : xIn  .. xTab    -- decode |x| into xScaled; compute flags & table addr
//                             (adr includes a 1-bit selector = xsgn that the
//                             coefficient ROM reads as the pos/neg sub-table)
//   tab  : xTab .. xCalc   -- the coefficient ROM, read at that address
//   poly : xCalc .. xPost  -- shared PolynomialEval evaluates the chosen sub-table
//   post : xPost .. xOut   -- normalize zres, apply xsgn-dependent exponent
//                             shift, fold in special-value overrides
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.RealSpec
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._

// -------------------------------------------------------------------------
//   info bundle piped from pre to post
// -------------------------------------------------------------------------

class SoftPlusPostProcInfo(val spec: RealSpec) extends Bundle {
  val xsgn   = UInt(1.W)
  val znan   = Bool()         // x is NaN
  val zinf   = Bool()         // x = +Inf
  val zzero  = Bool()         // x = -Inf, or x < -rangeMaxEx
  val xlarge = Bool()         // |x| > rangeMaxEx and x is finite (and x.sgn = 0)
  val xVal   = UInt(spec.W.W) // original x, used for the xlarge passthrough case
}

// -------------------------------------------------------------------------
//  pre-process
// -------------------------------------------------------------------------

class SoftPlusPreProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = bitW.addrW
  val fracW  = bitW.outBitW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = SoftPlusPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    // 1-bit pos/neg selector (= xsgn) prepended to the natural addrW bits
    val adr  = Valid(UInt((1 + adrW).W))
    val dx   = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    val info = Valid(new SoftPlusPostProcInfo(spec))
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  val rangeMaxLog2 = SoftPlusTable.rangeMaxLog2(fracW)
  val rangeMaxEx   = exBias + rangeMaxLog2 - 1

  val xExLargeEnough = rangeMaxEx.U(exW.W) < x.ex
  val xinf  = x.inf && !x.nan
  val znan  = x.nan
  val zinf  = xinf && (x.sgn === 0.U)
  val zzero = (xinf || xExLargeEnough) && (x.sgn === 1.U)
  val xlarge = xExLargeEnough && (x.sgn === 0.U) && !xinf && !znan

  val xScaleShift = (rangeMaxEx + 1).U(exW.W) - x.ex
  val xScaledComb = Cat(1.U(1.W), x.man) >> xScaleShift

  val infoComb = Wire(new SoftPlusPostProcInfo(spec))
  infoComb.xsgn   := x.sgn
  infoComb.znan   := znan
  infoComb.zinf   := zinf
  infoComb.zzero  := zzero
  infoComb.xlarge := xlarge
  infoComb.xVal   := Cat(x.sgn, x.ex, x.man)

  // --------------------------------------------------------------------------
  // stage "scale": the argument scaled into the table domain

  private val cutScale = new SharedPipe("scale", cutIn.valid, stages.at("scale"))

  val xScaled = cutScale(xScaledComb)
  val info    = cutScale(infoComb)

  // --------------------------------------------------------------------------
  // stage "out"

  // pos/neg selector goes in the MSB of the adr; it picks the sub-table.
  val adr = Cat(info.xsgn, TableAddress.adr(xScaled(manW - 1, 0), adrW))
  io.adr := Pipe(cutScale.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = TableAddress.dx(xScaled(manW - 1, 0), adrW, dxW)
    io.dx.get := Pipe(cutScale.valid, dx, stages.at("out"))
  }

  io.info := Pipe(cutScale.valid, info, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.SoftPlusPreProcess]]. */
object SoftPlusPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"scale"` the argument scaled into the table domain
   *  - `"out"`   the table address and the info
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("scale", 3.0),
    StagePoint("out",   1.0))
}

// -------------------------------------------------------------------------
//   post-process
// -------------------------------------------------------------------------

class SoftPlusPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val fracW     = bitW.outBitW
  val calcW     = fracW
  val extraBits = fracW - manW

  val rangeMaxLog2 = SoftPlusTable.rangeMaxLog2(fracW)

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = SoftPlusPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info = new SoftPlusPostProcInfo(spec)
      val zres = UInt(fracW.W)
    }))
    val z = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val in     = cutIn(io.in.bits)
  val infoIn = in.info
  val zres   = in.zres

  // ---- normalize: leftshift until top bit reaches position calcW ----
  // zShift = (calcW + 1) - binaryWidth(zres)
  //        = leading zeros count of zres (within fracW bits) + 1
  // Use +& so the +1 doesn't overflow when zres == 0 (the all-zero branch
  // returns `zres.getWidth.U` whose width is log2Up(fracW+1); for fracW that
  // is a power of two + 1, that width is just log2Up(fracW), and a plain
  // `+ 1.U` would wrap to 0).
  val zShiftComb = Mux(zres === 0.U, zres.getWidth.U,
                       PriorityEncoder(Reverse(zres))) +& 1.U
  // the shift is the leading zero count plus one, so the leading one lands on
  // bit calcW and nothing is shifted past it.
  val zmanW1Comb = (zres << zShiftComb)(calcW, 0)

  // --------------------------------------------------------------------------
  // stage "norm": the normalized mantissa

  private val cutNorm = new SharedPipe("norm", cutIn.valid, stages.at("norm"))

  val zmanW1   = cutNorm(zmanW1Comb)
  val zShift   = cutNorm(zShiftComb)
  val infoNorm = cutNorm(infoIn)

  // +& : the carry out of the round is what zmanRoundedMoreThan1 tests
  val zmanW1Rounded = if (extraBits > 0) {
    (zmanW1 >> extraBits) +& zmanW1(extraBits - 1)
  } else {
    zmanW1
  }

  val zmanRoundedMoreThan1 = zmanW1Rounded >= (BigInt(1) << (manW + 1)).U
  val zmanComb = Mux(zmanRoundedMoreThan1, Fill(manW, 1.U(1.W)),
                     zmanW1Rounded(manW - 1, 0))

  // ---- zEx: shifted bias, plus the pos branch's hidden (rangeMaxLog2+1) factor ----
  val zexComb = (exBias.U(exW.W) - zShift) +
                Mux(infoNorm.xsgn === 0.U, (rangeMaxLog2 + 1).U(exW.W), 0.U(exW.W))

  // --------------------------------------------------------------------------
  // stage "round": the rounded mantissa and the exponent

  private val cutRound = new SharedPipe("round", cutNorm.valid, stages.at("round"))

  val zman = cutRound(zmanComb)
  val zex  = cutRound(zexComb)
  val info = cutRound(infoNorm)

  // --------------------------------------------------------------------------
  // stage "out"

  // ---- special-value overrides ----
  val zsgn = 0.U(1.W)
  val zNanOrInf = Cat(zsgn, Fill(exW, 1.U(1.W)), info.znan.asUInt, 0.U((manW - 1).W))
  val zPassthrough = Cat(zsgn, info.xVal(spec.W - 2, 0)) // strip x's sign (always positive)

  val z = Mux(info.znan || info.zinf, zNanOrInf,
          Mux(info.zzero,             0.U(spec.W.W),
          Mux(info.xlarge,            zPassthrough,
                                      Cat(zsgn, zex, zman))))

  io.z := Pipe(cutRound.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.SoftPlusPostProcess]]. */
object SoftPlusPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"norm"`  the normalized mantissa
   *  - `"round"` the rounded mantissa and the exponent
   *  - `"out"`   the assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(
    StagePoint("norm",  2.0),
    StagePoint("round", 2.0),
    StagePoint("out",   1.0))
}
