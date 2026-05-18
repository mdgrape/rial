//% @file InvSqrt.scala
//
// inverse square root function (1/sqrt(x))
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import chisel3._
import chisel3.util._

import rial.approx._
import rial.arith.RealSpec
import rial.util.PipelineStageConfig

// An implementation of 1/sqrt using (normally 2nd-order) polynomial interpolation
//
// Since FP is composed of exponent and mantissa, the result of 1/sqrt will be
//     1/sqrt(2^ex * 1.m) = 2^(-ex/2) * 1/sqrt(1.m) .
// The fractional part of exponent appears when ex is an odd number, so
//     1/sqrt(2^ex * 1.m) = 2^(-floor(ex/2)-1) * 1/sqrt(1.m)   if ex is even
//                          2^(-floor(ex/2)-1) * 1/sqrt(2*1.m) otherwise.
// It requires 2 different tables for 1/sqrt(1.m) and 1/sqrt(2*1.m).
// To distinguish them, we use LSB of ex to detect whether ex is even or odd.

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

// invsqrt(x): floating => floating
// - if x < 0, returns inf.
class InvSqrtPreProcess(
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
    val adr = Valid(UInt((1+adrW).W)) // we use LSB of x.ex
    val dx  = if(order != 0) { Some(Valid(UInt(dxW.W))) } else { None }
  })

  // use the LSB of x.ex to distinguish 1~2 and 2~4
  val adr = Cat(io.x.bits.ex(0), io.x.bits.man(manW-1, dxW))
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

class InvSqrtTableCoeff(
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
class InvSqrtNonTable(
  val spec:  RealSpec, // Input / Output floating spec
  val bitW:  PolynomialBitWidth,
  val stage: Int,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val x = Flipped(Valid(new DecomposedReal(spec)))
    val z = Valid(new InvSqrtPostProcInfo(spec))
  })

  val out = WireDefault(0.U.asTypeOf(new InvSqrtPostProcInfo(spec)))
  io.z := Pipe(io.x.valid, out, stage)

  val x = io.x.bits

  val xneg  = if(spec.disableSign) {false.B} else {x.sgn.asBool}
  val znan  = x.nan
  val zinf  = x.zero || xneg
  val zzero = x.inf

  val zman0 = (x.man === 0.U) && (!x.zero) && (x.ex(0) === (exBias % 2).U(1.W))
  // if x.man == 0 && x.ex == 2N, that means x = 2^2N. then z = 2^-N, so zman = 0.

  val xExNobias  = x.ex - exBias.U
  val xExHalf    = Cat(xExNobias(exW-1), xExNobias(exW-1, 1)) // arith right shift
  val zEx0       = ~xExHalf // -(xex>>1)-1 = ~(xex>>1)+1-1 = ~(xex>>1)

  out.zex  := Mux(zinf || znan, Fill(exW, 1.U(1.W)),
              Mux(zzero,        0.U(exW.W),
                                zEx0 + zman0.asUInt + exBias.U))
  out.znan  := znan
  out.zinf  := zinf
  out.zman0 := zzero || zman0
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class InvSqrtPostProcInfo(val spec: RealSpec) extends Bundle {
  val zex   = UInt(spec.exW.W)
  val znan  = Bool()
  val zinf  = Bool()
  val zman0 = Bool()
}

class InvSqrtPostProcess(
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

  val nStage = stage.total

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info = new InvSqrtPostProcInfo(spec)
      val zres = UInt(fracW.W)
    }))
    val z = Valid(UInt(spec.W.W))
  })

  val zsgn  = 0.U
  val zex   = io.in.bits.info.zex
  val znan  = io.in.bits.info.znan
  val zinf  = io.in.bits.info.zinf
  val zman0 = io.in.bits.info.zman0

  val zIsNonTable  = znan || zinf || zman0
  val zmanNonTable = Cat(znan, 0.U((manW-1).W))

  val zmanRounded = WireDefault(0.U(manW.W))
  if(extraBits == 0) {
    zmanRounded := io.in.bits.zres
  } else {
    val zman0 = (io.in.bits.zres >> extraBits) +& io.in.bits.zres(extraBits-1)
    val polynomialOvf = zman0(manW) === 1.U
    zmanRounded := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW-1,0))
  }

  val zman = Mux(zIsNonTable, zmanNonTable, zmanRounded)
  val z = Cat(zsgn, zex, zman)

  io.z := Pipe(io.in.valid, z, nStage)
}
