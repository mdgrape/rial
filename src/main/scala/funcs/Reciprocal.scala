//% @file Reciprocal.scala
//
// Reciprocal function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import chisel3._
import chisel3.util._

import rial.approx._
import rial.arith.RealSpec
import rial.util.PipelineStageConfig

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
  val stage: Int,
) extends Module {

  val manW   = spec.manW

  val adrW   = bitW.addrW
  val fracW  = bitW.outBitW
  val dxW    = bitW.inBitW
  val order  = bitW.nOrder

  val io = IO(new Bundle {
    val x   = Flipped(Valid(new DecomposedReal(spec)))
    val adr = Valid(UInt(adrW.W))
    val dx  = if(order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
  })

  val adr = io.x.bits.man(manW-1, dxW)
  io.adr := Pipe(io.x.valid, adr, stage)

  if(order != 0) {
    // convert UInt 0000~1111 -> SInt 1000,..1111,0000,..0111
    val dx = Cat(~io.x.bits.man(dxW-1), io.x.bits.man(dxW-2, 0))
    io.dx.get := Pipe(io.x.valid, dx, stage)
  }
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class ReciprocalTableCoeff(
  val spec:  RealSpec,
  val coeff: PolynomialCoefficients,
  val stage: Int,
) extends Module {

  val manW   = spec.manW
  val adrW   = coeff.bitW.addrW
  val fracW  = coeff.bitW.outBitW
  val order  = coeff.bitW.nOrder

  val io = IO(new Bundle {
    val adr = Flipped(Valid(UInt(adrW.W)))
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
class ReciprocalNonTable(
  val spec:  RealSpec, // Input / Output floating spec
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val x = Flipped(Valid(new DecomposedReal(spec)))
    val z = Valid(new RoundingPostProcInfo(spec))
  })

  val out = WireDefault(0.U.asTypeOf(new RoundingPostProcInfo(spec)))
  io.z := Pipe(io.x.valid, out, stage)

  val x = io.x.bits

  val znan  = x.nan
  val zinf  = x.zero
  val zzero = x.inf

  val xmanNonZero = x.man.orR
  val zex0        = (exBias<<1).U((exW+1).W) - x.ex - xmanNonZero.asUInt
  val zexMSB      = zex0(exW)

  // Check overflow/underflow
  // 1/x = 2^(-e-1) * 2/1.m
  val invExMax = -spec.exMin - 1
  val invExMin = -spec.exMax - 1
  val (zinfOvf, zzeroUnf) = if(spec.exMax < invExMax) { // check overflow
    (zexMSB, zex0 === 0.U)
  } else if (invExMin < spec.exMin) { // check underflow
    (false.B, zexMSB || zex0 === 0.U)
  } else { // nothing happens.
    (false.B, zex0 === 0.U)
  }

  val zman0 = (x.man === 0.U) && (!x.zero)
  // if x.man == 0, that means x = 2^N. then z = 2^-N, so zman = 0.

  val zIsNonTable = znan || zinf || zzero || zinfOvf || zzeroUnf || zman0

  out.zsgn := Mux(znan || zzero || zinf || zinfOvf || zzeroUnf, 0.U, x.sgn)
  out.zex  := Mux(znan || zinf || zinfOvf, Fill(exW, 1.U(1.W)),
              Mux(zzero || zzeroUnf,        0.U(exW.W),
                                            zex0(exW-1,0)))
  out.znan := znan
  out.zinf := zinf || zinfOvf
  out.zIsNonTable := zIsNonTable
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
