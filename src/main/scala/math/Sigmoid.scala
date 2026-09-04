//% @file Sigmoid.scala
//
// Sigmoid (1 / (1 + exp(-x))) hardware modules.
//
//   pre  : xIn  .. xTab    -- decode |x| into xScaled, compute flags & table addr
//   tab  : xTab .. xCalc   -- standard polynomial coefficient ROM
//   poly : xCalc .. xPost  -- shared PolynomialEval evaluates `2 - 2/(1+exp(-|x|))`
//   post : xPost .. xOut   -- fold polynomial into the final mantissa, apply
//                             special-value (NaN / saturation / near-zero) overrides
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

class SigmoidPostProcInfo(val spec: RealSpec) extends Bundle {
  val xsgn  = UInt(1.W)
  val xzero = Bool() // xScaled <= 1, sigmoid result is ~0.5
  val znan  = Bool() // input is NaN
  val zone  = Bool() // x >>  rangeMax  -> sigmoid = 1
  val zzero = Bool() // x << -rangeMax  -> sigmoid = 0
}

// -------------------------------------------------------------------------
//  pre-process
// -------------------------------------------------------------------------

class SigmoidPreProcess(
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
  val stagePoints = SigmoidPreProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    val adr  = Valid(UInt(adrW.W))
    val dx   = if (order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
    val info = Valid(new SigmoidPostProcInfo(spec))
  })

  private val cutIn = new SharedPipe("in", io.x.valid, stages.at("in"))

  val x = cutIn(io.x.bits)

  // rangeMaxLog2 must match SigmoidSimGeneric -- both derive it from fracW.
  val rangeMaxLog2 = SigmoidTable.rangeMaxLog2(fracW)
  val rangeMaxEx   = exBias + rangeMaxLog2 - 1

  val xExLargeEnough = rangeMaxEx.U(exW.W) < x.ex

  val znan  =  x.nan
  val zone  = (x.sgn === 0.U) && xExLargeEnough && !znan
  val zzero = (x.sgn === 1.U) && xExLargeEnough && !znan

  val xScaleShift = (rangeMaxEx + 1).U(exW.W) - x.ex
  val xScaled     = Cat(1.U(1.W), x.man) >> xScaleShift

  val adr = xScaled(manW - 1, dxW)
  io.adr := Pipe(cutIn.valid, adr, stages.at("out"))

  if (order != 0) {
    val dx = Cat(~xScaled(dxW - 1), xScaled(dxW - 2, 0))
    io.dx.get := Pipe(cutIn.valid, dx, stages.at("out"))
  }

  val info = Wire(new SigmoidPostProcInfo(spec))
  info.xsgn  := x.sgn
  info.xzero := xScaled === 0.U || xScaled === 1.U
  info.znan  := znan
  info.zone  := zone
  info.zzero := zzero
  io.info := Pipe(cutIn.valid, info, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.SigmoidPreProcess]]. */
object SigmoidPreProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"out"` the table address and the info
   */
  val stagePoints: Seq[StageElem] = Seq(StagePoint("out", 4.0))
}

// -------------------------------------------------------------------------
//   table coefficient ROM
// -------------------------------------------------------------------------

class SigmoidTableCoeff(
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
//   post-process
// -------------------------------------------------------------------------

class SigmoidPostProcess(
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

  /** Boundaries at which pipeline registers can be placed. */
  val stagePoints = SigmoidPostProcess.stagePoints

  val stages = stage.resolve(stagePoints)

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info = new SigmoidPostProcInfo(spec)
      val zres = UInt(fracW.W)
    }))
    val z = Valid(UInt(spec.W.W))
  })

  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val in   = cutIn(io.in.bits)
  val info = in.info
  val zres = in.zres

  // ---- fold:  zPos = (1 << (calcW+1)) - zres   for x >= 0
  //             zNeg =                      zres   for x <  0
  val zPos = (BigInt(1) << (calcW + 1)).U((calcW + 2).W) - zres
  val zNeg = zres
  val zRes = Mux(info.xsgn === 0.U, zPos, zNeg)

  // ---- normalize: leftshift until top bit is at position (calcW + 1) ----
  val zShift = Mux(zRes === 0.U, zRes.getWidth.U, PriorityEncoder(Reverse(zRes)))
  val zmanW1 = zRes << zShift
  val zmanW1Rounded = (zmanW1 >> (extraBits + 1)) + zmanW1(extraBits)

  val zmanRoundedMoreThan1 = zmanW1Rounded >= (BigInt(1) << (manW + 1)).U
  val zman = Mux(zmanRoundedMoreThan1, Fill(manW, 1.U(1.W)), zmanW1Rounded(manW - 1, 0))

  val zex = exBias.U(exW.W) - zShift

  // ---- special-value overrides ----
  val znan  = info.znan
  val zone  = info.zone
  val zzero = info.zzero
  val xzero = info.xzero

  val zsgn = 0.U(1.W)
  val z = Mux(znan,  Cat(zsgn, Fill(exW, 1.U(1.W)), 1.U(1.W), 0.U((manW - 1).W)),
          Mux(zone,  Cat(zsgn, exBias.U(exW.W),     0.U(manW.W)),
          Mux(zzero, 0.U(spec.W.W),
          Mux(xzero, Cat(zsgn, (exBias - 1).U(exW.W), 0.U(manW.W)),
                     Cat(zsgn, zex, zman)))))

  io.z := Pipe(cutIn.valid, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.math.SigmoidPostProcess]]. */
object SigmoidPostProcess {

  /** Boundaries at which pipeline registers can be placed.
   *
   *  - `"out"` the folded and assembled result
   */
  val stagePoints: Seq[StageElem] = Seq(StagePoint("out", 3.0))
}
