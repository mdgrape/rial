//% @file threefry.scala
//
// Threefry Random Number Generation
//   based on J. Salomon et al.
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.rng

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._

//
// take random number from Threefry using ready/valid interface
//
class Threefry4x32Generator(
  r: Int = 20,
  rotStage: Int = 0,
) extends Module {

  val io = IO(new Bundle {
    val initialized = Output(Bool())
    val init = new Bundle {
      val en    = Input(Bool())
      val key   = Input(Vec(4, UInt(32.W)))
      val count = Input(Vec(4, UInt(32.W)))
    }
    val rand = Decoupled(Vec(4, UInt(32.W)))
  })

  val rng = Module(new Threefry4_32(r, rotStage))

  def getParam = { rng.getParam }
  def nStages  = { rng.nStages + 1 }

  val initialized = RegInit(false.B)
  io.initialized := initialized

  val key   = RegInit(VecInit.fill(4)(0.U(32.W)))
  val count = RegInit(VecInit.fill(4)(0.U(32.W)))
  when(io.init.en) {
    key   := io.init.key
    count := io.init.count
    initialized := true.B
  }

  val outQ = Module(new Queue(Vec(4, UInt(32.W)), 1, pipe=true, flow=false))
  io.rand <> outQ.io.deq

  val step = initialized && outQ.io.enq.ready

  rng.io.en          := step
  rng.io.input.valid := initialized
  rng.io.input.key   := key
  rng.io.input.count := count

  // update counter only when we need to push pipeline forward
  val countNext = WireDefault(VecInit(count.map(c => c + 1.U)))
  when(step) {
    count(0) := countNext(0)
    when(countNext(0) === 0.U) {
      count(1) := countNext(1)
      when(countNext(1) === 0.U) {
        count(2) := countNext(2)
        when(countNext(2) === 0.U) {
          count(3) := countNext(3)
        }
      }
    }
  }
  outQ.io.enq.valid := rng.io.output.valid
  outQ.io.enq.bits  := rng.io.output.rand

  // printf("Threefry4x32Gen: init = %b, out.ready = %b, out.valid = %b, step = %b\n",
  //   initialized, io.rand.ready, io.rand.valid, step)
  // printf("Threefry4x32Gen: rng.en = %b, input.valid = %b, output.valid = %b\n",
  //   rng.io.en, rng.io.input.valid, rng.io.output.valid)
  // printf("Threefry4x32Gen: key = [%x,%x,%x,%x], cnt = [%x,%x,%x,%x]\n",
  //   key(0), key(1), key(2), key(3), count(0), count(1), count(2), count(3))
  // printf("Threefry4x32Gen: out = [%x,%x,%x,%x]\n",
  //   rng.io.output.rand(0), rng.io.output.rand(1), rng.io.output.rand(2), rng.io.output.rand(3))
}

// Fully pipelined implementation, no loop
// r : Number of rotations
// rotStage : Insert pipeline stage every rotStage rotation
class Threefry4_32(
  r: Int = 20,
  rotStage: Int = 0
) extends Module {

  def getParam = { (r, rotStage) }
  def rotL( x : UInt, n : Int ) = { Cat( x.tail(n), x.head(n) ) }

  def nStages = {
    if(rotStage == 0) {
      0
    } else {
      Seq.tabulate(r)(j => {
        if( ((j+1) % rotStage) == 0 ) { 1 } else { 0 }
      }).fold(0)(_+_)
    }
  }

  // -------------------------------------------------------------------------

  val io = IO(new Bundle {
    val en    = Input(Bool())
    val input = new Bundle {
      val valid = Input(Bool())
      val key   = Input(Vec(4, UInt(32.W)))
      val count = Input(Vec(4, UInt(32.W)))
    }
    val output = new Bundle {
      val valid = Output(Bool())
      val rand  = Output(Vec(4, UInt(32.W)))
    }
  })

  val keys  = WireDefault(VecInit.fill(5  )(0.U(32.W)))
  val x     = WireDefault(VecInit.fill(r+1)(VecInit.fill(4)(0.U(32.W))))
  val xNext = WireDefault(VecInit.fill(r  )(VecInit.fill(4)(0.U(32.W))))

  // -------------------------------------------------------------------------
  // init key

  //val Threefish_C240 = 0x1BD11BDAA9FC1A22
  val c240 = 0x1BD11BDA.U

  for(i <- 0 until 4) {
    keys(i) := io.input.key(i)
  }
  keys(4) := io.input.key.foldLeft(c240)( (z,k) => z ^ k )

  // -------------------------------------------------------------------------
  // set xNext

  // Rotation constant
  val r0 = Array[Int]( 10, 11, 13, 23, 6, 17, 25, 18 )
  val r1 = Array[Int]( 26, 21, 27, 5, 20, 11, 10, 20 )

  for(j <- 0 until r) {
    if (j % 2 == 0) {
      xNext(j)(0) :=       x(j)(0)              + x    (j)(1)
      xNext(j)(1) := rotL( x(j)(1), r0(j % 8) ) ^ xNext(j)(0)
      xNext(j)(2) :=       x(j)(2)              + x    (j)(3)
      xNext(j)(3) := rotL( x(j)(3), r1(j % 8) ) ^ xNext(j)(2)
    } else {
      xNext(j)(0) :=       x(j)(0)              + x    (j)(3)
      xNext(j)(1) := rotL( x(j)(1), r1(j % 8) ) ^ xNext(j)(2)
      xNext(j)(2) :=       x(j)(2)              + x    (j)(1)
      xNext(j)(3) := rotL( x(j)(3), r0(j % 8) ) ^ xNext(j)(0)
    }
  }

  // -------------------------------------------------------------------------
  // set x

  val v = WireDefault(VecInit.fill(r+1)(false.B))
  v(0) := io.input.valid

  for(i <- 0 until 4) {
    x(0)(i) := io.input.count(i) + io.input.key(i)
  }

  for(j <- 0 until r) {

    val xj = WireDefault(VecInit.fill(4)(0.U(32.W)))
    if (j % 4 == 3) {
      val keyR = (j/4) + 1
      xj(0) := xNext(j)(0) + keys((keyR+0)%5)
      xj(1) := xNext(j)(1) + keys((keyR+1)%5)
      xj(2) := xNext(j)(2) + keys((keyR+2)%5)
      xj(3) := xNext(j)(3) + keys((keyR+3)%5) + keyR.U
    } else {
      xj    := xNext(j)
    }

    val addStage = (rotStage != 0) && ( ((j+1) % rotStage) == 0 )
    if(addStage) {
      val vReg = RegInit(false.B)
      val xReg = RegInit(VecInit.fill(4)(0.U(32.W)))
      when(io.en) {
        xReg := xj
        vReg := v(j)
      }
      x(j+1) := xReg
      v(j+1) := vReg
    } else {
      x(j+1) := xj
      v(j+1) := v(j)
    }
  }

  //for(j <- 0 to r) {
  //printf("%d %x %x %x %x %x\n", j.U, X(j)(0), X(j)(1), X(j)(2), X(j)(3), rotL( X(j)(1), R0(j%8) ))
  //}
  io.output.valid := v(r)
  io.output.rand  := x(r)
}

object Threefry4_32Driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
      Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new Threefry4_32(20,2)))
    )
}
