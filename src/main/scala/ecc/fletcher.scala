//% @file crc.scala
//
// Fletcher's checksum
//
// part of RIAL - RIKEN Arithmetic and Logic library
// Copyright (C) Toru Niina RIKEN BDR 2025
//
package rial.ecc

import chisel3._
import chisel3.util._

// Fletcher-16
//
// Ai = Ai-1 + di
// Bi = Bi-1 + Ai
// checksum = { (Bi mod 255), (Ai mod 255) }
//
// A0 = 0                  | B0 = 0
// A1 = d1                 | B1 = d1
// A2 = d1 + d2            | B2 = d1*2 + d2
// A3 = d1 + d2 + d3       | B3 = d1*3 + d2*2 + d3
// A4 = d1 + d2 + d3 + d4  | B4 = d1*4 + d2*3 + d3*2 + d4
// ...                     | ...
// An = sum({di})          | Bn = d1*n + d2*(n-1) + ... + dn
//
// max(An) = 255 * n
// max(Bn) = 255 * (n + n-1 + ... + 1) = 255 * (n * (n+1) / 2)
//
// (0 <= xi <= 255)
// x = x0 + x1 * 256     + x2 * 256^2     + ... xn * 256^n     (mod255)
//   = x0 + x1 * (1+255) + x2 * (1+255)^2 + ... xn * (1+255)^n (mod255)
//   = x0 + x1           + x2             + ... xn             (mod255)
//
// if 0 <= x <= 254, mod255(x) == x. note that mod255(0xff) == 0
//

class Fletcher(
  checksumW: Int, // if 16, use Fletcher-16.
  inputW:    Int, // if not a multiple of checksumW/2, it will be 0-padded.
  latency:   Int  // latency.
) extends Module {

  assert(checksumW % 2 == 0)
  val elemW = checksumW / 2
  val mod = (1 << elemW) - 1

  def splitIntoElems(x: UInt): Seq[UInt] = {
    // 0-pad to make it a multiple of elemW
    val y = x.pad((x.getWidth + elemW - 1) / elemW * elemW)
    Seq.tabulate(y.getWidth / elemW)(i => { y((i+1)*elemW-1, i*elemW) })
  }

  def modsum(lhs: UInt, rhs: UInt): UInt = {
    assert(lhs.getWidth == elemW)
    assert(rhs.getWidth == elemW)

    val sum = lhs +& rhs
    sum(elemW-1,0) + sum(elemW)
  }

  // -------------------------------------------------------------------------

  val io = IO(new Bundle {
    val in  = Input(new Bundle {
      val valid = Bool()
      val bits  = UInt(inputW.W)
    })
    val out = Output(new Bundle {
      val valid = Bool()
      val bits  = UInt(checksumW.W)
    })
  })

  val d = VecInit(splitIntoElems(io.in.bits))
  val n = d.length

  val maxBn = ((1 << elemW) - 1) * (n * (n+1)) / 2
  val bnW = log2Ceil(maxBn+1)
  val bn = VecInit(d.zipWithIndex.map{
    case(v, i) => { (v * (d.length-i).U).pad(bnW) }
  }).reduceTree(_+_)

  val amod0 = d.reduceTree(modsum(_,_))
  val bmod0 = VecInit(splitIntoElems(bn)).reduceTree(modsum(_,_))

  val amod = Mux(amod0 === mod.U, 0.U(elemW.W), amod0)
  val bmod = Mux(bmod0 === mod.U, 0.U(elemW.W), bmod0)
  val checksum = Cat(bmod, amod)

  assert(amod.getWidth == elemW)
  assert(bmod.getWidth == elemW)
  assert(checksum.getWidth == checksumW)

  io.out.valid := ShiftRegister(io.in.valid, latency, false.B, true.B)
  io.out.bits  := ShiftRegister(checksum,    latency, 0.U,     true.B)
}
