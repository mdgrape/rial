//% @file sqrt.scala
//
// square root function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import chisel3._
import chisel3.util._

import rial.approx._
import rial.arith.RealSpec
import rial.util.PipelineStageConfig

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
  val coeff: PolynomialCoefficients,
  val stage: PipelineStageConfig,
) extends Module {

  val nStage = stage.total
  val manW   = spec.manW

  val adrW   = coeff.bitW.addrW
  val fracW  = coeff.bitW.outBitW
  val dxW    = coeff.bitW.inBitW
  val order  = coeff.bitW.nOrder

  val io = IO(new Bundle {
    val x   = Flipped(Valid(new DecomposedReal(spec)))
    val adr = Output(UInt((1+adrW).W)) // we use LSB of x.ex
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  // use the LSB of x.ex to distinguish 1~2 and 2~4
  val adr = Cat(io.x.bits.ex(0), io.x.bits.man(manW-1, dxW))
  io.adr := Pipe(io.x.valid, adr, nStage)

  if(order != 0) {
    // convert UInt 0000~1111 -> SInt 1000,..1111,0000,..0111
    val dx = Cat(~io.x.bits.man(dxW-1), io.x.bits.man(dxW-2, 0))
    io.dx.get := Pipe(io.x.valid, dx, nStage)
  }
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
  val stage: PipelineStageConfig,
) extends Module {

  val nStage = stage.total

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
  io.cs := Pipe(io.adr.valid, packed(io.adr.bits), nStage)
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
  val coeff: PolynomialCoefficients,
  val stage: PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val x = Flipped(Valid(new DecomposedReal(spec)))
    val z = Valid(new RoundingPostProcInfo(spec))
  })

  val out = WireDefault(0.U.asTypeOf(new RoundingPostProcInfo(spec)))
  io.z := Pipe(io.x.valid, out, nStage)

  val x = io.x.bits

  val xneg  = if(spec.disableSign) {false.B} else {x.sgn.asBool}
  val znan  = x.nan
  val zinf  = x.inf
  val zzero = x.zero || xneg

  out.zsgn := 0.U(1.W)

  val xExNobias  = x.ex - exBias.U
  val zExShifted = Cat(xExNobias(exW-1), xExNobias(exW-1, 1)) // arith right shift
  out.zex  := Mux(zinf || znan, Fill(exW, 1.U(1.W)),
              Mux(zzero,        0.U(exW.W),
                                zExShifted + exBias.U))
  out.znan := znan
  out.zinf := zinf
  out.zIsNonTable := znan || zinf || zzero
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

//
// see RoundingPostProc.
//
