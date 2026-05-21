//% @file SoftPlus.scala
//
// SoftPlus (log(1 + exp(x))) hardware modules.
//
//   pre  : xIn  .. xTab    -- decode |x| into xScaled; compute flags & table addr
//                             (adr includes a 1-bit selector = xsgn that the
//                             TableCoeff uses to pick the pos/neg sub-table)
//   tab  : xTab .. xCalc   -- concat'd pos+neg coefficient ROM, Mux'd by xsgn
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
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = bitW.addrW
  val fracW  = bitW.outBitW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  val io = IO(new Bundle {
    val x    = Flipped(Valid(new DecomposedReal(spec)))
    // 1-bit pos/neg selector (= xsgn) prepended to the natural addrW bits
    val adr  = Valid(UInt((1 + adrW).W))
    val dx   = if (order != 0) Some(Valid(UInt(dxW.W))) else None
    val info = Valid(new SoftPlusPostProcInfo(spec))
  })

  val rangeMaxLog2 = SoftPlusTable.rangeMaxLog2(fracW)
  val rangeMaxEx   = exBias + rangeMaxLog2 - 1

  val xExLargeEnough = rangeMaxEx.U(exW.W) < io.x.bits.ex
  val xinf  = io.x.bits.inf && !io.x.bits.nan
  val znan  = io.x.bits.nan
  val zinf  = xinf && (io.x.bits.sgn === 0.U)
  val zzero = (xinf || xExLargeEnough) && (io.x.bits.sgn === 1.U)
  val xlarge = xExLargeEnough && (io.x.bits.sgn === 0.U) && !xinf && !znan

  val xScaleShift = (rangeMaxEx + 1).U(exW.W) - io.x.bits.ex
  val xScaled     = Cat(1.U(1.W), io.x.bits.man) >> xScaleShift

  // pos/neg selector goes in the MSB of the adr; TableCoeff uses it to mux.
  val adr = Cat(io.x.bits.sgn, xScaled(manW - 1, dxW))
  io.adr := Pipe(io.x.valid, adr, stage)

  if (order != 0) {
    val dx = Cat(~xScaled(dxW - 1), xScaled(dxW - 2, 0))
    io.dx.get := Pipe(io.x.valid, dx, stage)
  }

  val info = Wire(new SoftPlusPostProcInfo(spec))
  info.xsgn   := io.x.bits.sgn
  info.znan   := znan
  info.zinf   := zinf
  info.zzero  := zzero
  info.xlarge := xlarge
  info.xVal   := Cat(io.x.bits.sgn, io.x.bits.ex, io.x.bits.man)
  io.info := Pipe(io.x.valid, info, stage)
}

// -------------------------------------------------------------------------
//   table coefficient ROM (concat'd pos+neg, Mux'd by adr's MSB)
// -------------------------------------------------------------------------

class SoftPlusTableCoeff(
  val spec:  RealSpec,
  // concat [pos coeffs] ++ [neg coeffs]
  val coeff: PolynomialCoefficients,
  val stage: Int,
) extends Module {

  val adrW   = coeff.bitW.addrW
  val coefBW = coeff.bitW.coefBits.reduce(_+_)

  val io = IO(new Bundle {
    val adr = Flipped(Valid(UInt((1 + adrW).W)))
    val cs  = Valid(UInt(coefBW.W))
  })

  val all = PolynomialGenerator.packTable(coeff)
  require(all.length == 2 * (1 << adrW),
    s"expected 2 * 2^${adrW} softplus table entries, got ${all.length}")

  val n     = 1 << adrW
  val vecP  = VecInit(all.slice(0, n        ).map(c => c.U(coefBW.W)))
  val vecN  = VecInit(all.slice(n, 2 * n    ).map(c => c.U(coefBW.W)))

  val selNeg = io.adr.bits(adrW)       // MSB = xsgn
  val subAdr = io.adr.bits(adrW - 1, 0)

  val cs = Mux(selNeg, vecN(subAdr), vecP(subAdr))
  io.cs := Pipe(io.adr.valid, cs, stage)
}

// -------------------------------------------------------------------------
//   post-process
// -------------------------------------------------------------------------

class SoftPlusPostProcess(
  val spec:  RealSpec,
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val fracW     = bitW.outBitW
  val calcW     = fracW
  val extraBits = fracW - manW

  val rangeMaxLog2 = SoftPlusTable.rangeMaxLog2(fracW)

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info = new SoftPlusPostProcInfo(spec)
      val zres = UInt(fracW.W)
    }))
    val z = Valid(UInt(spec.W.W))
  })

  val info = io.in.bits.info
  val zres = io.in.bits.zres

  // ---- normalize: leftshift until top bit reaches position calcW ----
  // zShift = (calcW + 1) - binaryWidth(zres)
  //        = leading zeros count of zres (within fracW bits) + 1
  // Use +& so the +1 doesn't overflow when zres == 0 (the all-zero branch
  // returns `zres.getWidth.U` whose width is log2Up(fracW+1); for fracW that
  // is a power of two + 1, that width is just log2Up(fracW), and a plain
  // `+ 1.U` would wrap to 0).
  val zShift = Mux(zres === 0.U, zres.getWidth.U, PriorityEncoder(Reverse(zres))) +& 1.U
  val zmanW1 = zres << zShift

  val zmanW1Rounded = if (extraBits > 0) {
    (zmanW1 >> extraBits) + zmanW1(extraBits - 1)
  } else {
    zmanW1
  }

  val zmanRoundedMoreThan1 = zmanW1Rounded >= (BigInt(1) << (manW + 1)).U
  val zman = Mux(zmanRoundedMoreThan1, Fill(manW, 1.U(1.W)), zmanW1Rounded(manW - 1, 0))

  // ---- zEx: shifted bias, plus the pos branch's hidden (rangeMaxLog2+1) factor ----
  val zex = (exBias.U(exW.W) - zShift) +
            Mux(info.xsgn === 0.U, (rangeMaxLog2 + 1).U(exW.W), 0.U(exW.W))

  // ---- special-value overrides ----
  val zsgn = 0.U(1.W)
  val zNanOrInf = Cat(zsgn, Fill(exW, 1.U(1.W)), info.znan.asUInt, 0.U((manW - 1).W))
  val zPassthrough = Cat(zsgn, info.xVal(spec.W - 2, 0)) // strip x's sign (always positive)

  val z = Mux(info.znan || info.zinf, zNanOrInf,
          Mux(info.zzero,             0.U(spec.W.W),
          Mux(info.xlarge,            zPassthrough,
                                      Cat(zsgn, zex, zman))))

  io.z := Pipe(io.in.valid, z, stage)
}
