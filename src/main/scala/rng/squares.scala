package rial.rng

// Squares: A Fast Counter-Based RNG
//   Bernard Widynski, 2022
// https://arxiv.org/pdf/2004.06278.pdf

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._

/** Implementation of another counter-based RNG: Squares
 *
 * This implementation is based on  B. Widynski (2022)
 *
 * It keeps key and count inside the module and sequentially generates random
 * number.
 *
 * {{{
 * class Squares32Generator(...) extends Moudle {
 *   val io = IO(new Bundle {
 *     val initialized = Output(Bool())
 *     val init = new Bundle {
 *       val en    = Input(Bool())
 *       val key   = Input(UInt(64.W))
 *       val count = Input(UInt(64.W))
 *     }
 *     val rand = Decoupled(UInt(32.W))
 *   })
 *   // ...
 * }
 * }}}
 *
 * @constructor create a chisel Module. No parameter.
 */
class Squares32Generator() extends Module {

  val io = IO(new Bundle {
    val initialized = Output(Bool())
    val init = new Bundle {
      val en    = Input(Bool())
      val key   = Input(UInt(64.W))
      val count = Input(UInt(64.W))
    }
    val rand = Decoupled(UInt(32.W))
  })

  val rng = Module(new Squares32)

  def nStage = { rng.nStage + 1 }

  val initialized = RegInit(false.B)
  io.initialized := initialized

  val key   = RegInit(0.U(64.W))
  val count = RegInit(0.U(64.W))
  when(io.init.en) {
    key   := io.init.key
    count := io.init.count
    initialized := true.B
  }

  val outQ = Module(new Queue(UInt(32.W), 1, pipe=true, flow=false))
  io.rand <> outQ.io.deq

  val step = initialized && outQ.io.enq.ready

  rng.io.en          := step
  rng.io.input.valid := initialized
  rng.io.input.key   := key
  rng.io.input.count := count

  // update counter only when we need to push pipeline forward
  val countNext = count + 1.U
  when(step) {
    count := countNext
  }

  outQ.io.enq.valid := rng.io.output.valid
  outQ.io.enq.bits  := rng.io.output.rand
}

/** Fully pipelined implementation of Squares32, no loop
 *
 * It generates random number from input key and count.
 * It does not have internal state.
 *
 * {{{
 * class Squares32(...) extends Moudle {
 *   val io = IO(new Bundle {
 *     val en = Input(Bool())
 *                                      
 *     val input = new Bundle {
 *       val valid = Input (Bool())
 *       val key   = Input (UInt(64.W))
 *       val count = Input (UInt(64.W))
 *     }
 *                                      
 *     val output = new Bundle {
 *       val valid = Output(Bool())
 *       val rand  = Output(UInt(32.W))
 *     }
 *   })
 *   // ...
 * }
 * }}}
 *
 * @constructor create a chisel Module.
 */
class Squares32 extends Module {

  def nStage = 4

  val io = IO(new Bundle {
    val en = Input(Bool())

    val input = new Bundle {
      val valid = Input (Bool())
      val key   = Input (UInt(64.W))
      val count = Input (UInt(64.W))
    }

    val output = new Bundle {
      val valid = Output(Bool())
      val rand  = Output(UInt(32.W))
    }
  })

  val x0 = WireDefault(0.U(64.W))
  val y0 = WireDefault(0.U(64.W))
  val z0 = WireDefault(0.U(64.W))
  x0 := io.input.key * io.input.count
  y0 := x0
  z0 := io.input.key + y0

  val stage0 = Module(new Squares32RNGStage)
  val stage1 = Module(new Squares32RNGStage)
  val stage2 = Module(new Squares32RNGStage)
  val stage3 = Module(new Squares32RNGStage)

  stage0.io.en := io.en
  stage0.io.input.valid := io.input.valid
  stage0.io.input.x     := x0
  stage0.io.input.y     := y0
  stage0.io.input.c     := z0

  stage1.io.en := io.en
  stage1.io.input.valid := stage0.io.output.valid
  stage1.io.input.x     := stage0.io.output.z
  stage1.io.input.y     := stage0.io.output.c
  stage1.io.input.c     := stage0.io.output.y

  stage2.io.en := io.en
  stage2.io.input.valid := stage1.io.output.valid
  stage2.io.input.x     := stage1.io.output.z
  stage2.io.input.y     := stage1.io.output.c
  stage2.io.input.c     := stage1.io.output.y

  stage3.io.en := io.en
  stage3.io.input.valid := stage2.io.output.valid
  stage3.io.input.x     := stage2.io.output.z
  stage3.io.input.y     := stage2.io.output.c
  stage3.io.input.c     := stage2.io.output.y

  io.output.valid := stage3.io.output.valid
  io.output.rand  := stage3.io.output.z(31, 0)
}

private[rial] class Squares32RNGStage extends Module {

  val io = IO(new Bundle {
    val en = Input(Bool())

    val input = new Bundle {
      val valid = Input(Bool())
      val x = Input(UInt(64.W))
      val y = Input(UInt(64.W))
      val c = Input(UInt(64.W))
    }
    val output = new Bundle {
      val valid = Output(Bool())
      val z     = Output(UInt(64.W))
      val y     = Output(UInt(64.W))
      val c     = Output(UInt(64.W))
    }
  })

  val xsq = WireDefault(0.U(64.W))
  val z0  = WireDefault(0.U(64.W))
  val z   = WireDefault(0.U(64.W))
  xsq := io.input.x * io.input.x
  z0  := xsq + io.input.y
  z   := Cat(z0(31,0), z0(63, 32))

  val vReg = RegInit(false.B)
  val zReg = RegInit(0.U(64.W))
  val yReg = RegInit(0.U(64.W))
  val cReg = RegInit(0.U(64.W))

  when(io.en) {
    vReg := io.input.valid
    zReg := z
    yReg := io.input.y
    cReg := io.input.c
  }

  io.output.valid := vReg
  io.output.z     := zReg
  io.output.y     := yReg
  io.output.c     := cReg
}
