package rial.rng

// Squares: A Fast Counter-Based RNG
//   Bernard Widynski1, 2022
// https://arxiv.org/pdf/2004.06278.pdf

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._

class Squares32 extends Module {

  val io = IO(new Bundle {
    val key   = Input (UInt(64.W))
    val count = Input (UInt(64.W))
    val rand  = Output(UInt(32.W))
  })

  val x0 = WireDefault(0.U(64.W))
  val y0 = WireDefault(0.U(64.W))
  val z0 = WireDefault(0.U(64.W))
  x0 := io.key * io.count
  y0 := x0
  z0 := io.key + y0

  val stage0 = Module(new Squares32RNGStage)
  val stage1 = Module(new Squares32RNGStage)
  val stage2 = Module(new Squares32RNGStage)
  val stage3 = Module(new Squares32RNGStage)

  stage0.io.x := x0
  stage0.io.y := y0

  val x1 = stage0.io.out
  val z1 = ShiftRegister(1, z0, 0.U, true.B)
  stage1.io.x := x1
  stage1.io.y := z1

  val x2 = stage1.io.out
  val y2 = ShiftRegister(2, y0, 0.U, true.B)
  stage2.io.x := x2
  stage2.io.y := y2

  val x3 = stage2.io.out
  val z3 = ShiftRegister(2, z1, 0.U, true.B)
  stage3.io.x := x3
  stage3.io.y := z3

  io.rand := stage3.io.out(31, 0)
}

class Squares32RNGStage extends Module {

  val io = IO(new Bundle {
    val x   = Input(UInt(64.W))
    val y   = Input(UInt(64.W))
    val out = Output(UInt(64.W))
  })

  val xsq = WireDefault(0.U(64.W))
  xsq := io.x * io.x

  val z0 = WireDefault(0.U(64.W))
  z0 := xsq + io.y

  val z = RegInit(0.U(64.W))
  z := Cat(z(31,0), z(63, 32))

  io.out := z
}
