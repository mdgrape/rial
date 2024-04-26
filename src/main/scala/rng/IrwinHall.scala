package rial.rng

import rial.arith.RealSpec

// generate approx-gaussian random number based on Irwin-Hall distribution.

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._

class IrwinHallGenerator(
  val spec: RealSpec = RealSpec.Float32Spec,
  val rndW: Int = 32
) extends Module {

  val io = IO(new Bundle {
    val input  = Flipped(Decoupled(Input(Vec(12, UInt(rndW.W)))))
    val output = Decoupled(UInt(spec.W.W))
  })

  val rng = Module(new IrwinHall(spec, rndW))

  val inQ = Module(new Queue(Vec(12, UInt(rndW.W)), 1, pipe=true, flow=false))
  inQ.io.enq <> io.input

  val outQ = Module(new Queue(UInt(spec.W.W), 1, pipe=true, flow=false))
  io.output <> outQ.io.deq

  // if outQ is ready, then we can keep output from rng. push rng pipe forward.
  val step = outQ.io.enq.ready
  inQ.io.deq.ready := step

  rng.io.en          := step
  rng.io.input.valid := inQ.io.deq.valid
  rng.io.input.rnds  := inQ.io.deq.bits

  outQ.io.enq.valid := rng.io.output.valid
  outQ.io.enq.bits  := rng.io.output.rnd
}

class IrwinHall(
  spec: RealSpec = RealSpec.Float32Spec,
  rndW: Int = 32
) extends Module {

  assert(rndW > spec.manW)

  val io = IO(new Bundle {
    val en    = Input(Bool())
    val input = new Bundle {
      val valid = Input(Bool())
      val rnds  = Input(Vec(12, UInt(rndW.W)))
    }
    val output = new Bundle {
      val valid = Output(Bool())
      val rnd   = Output(UInt(spec.W.W))
    }
  })

  def nStage = { 4 + 1 } // 4 for TreeSum, 1 for FPConv

  val node0 = Module(new IrwinHallTreeSumNode(12, rndW))
  val node1 = Module(new IrwinHallTreeSumNode( 6, rndW+1))
  val node2 = Module(new IrwinHallTreeSumNode( 3, rndW+2))
  val node3 = Module(new IrwinHallTreeSumNode( 2, rndW+3))

  node0.io.en := io.en
  node0.io.input.valid := io.input.valid
  node0.io.input.rnds  := io.input.rnds

  node1.io.en := io.en
  node1.io.input.valid := node0.io.output.valid
  node1.io.input.rnds  := node0.io.output.rnds

  node2.io.en := io.en
  node2.io.input.valid := node1.io.output.valid
  node2.io.input.rnds  := node1.io.output.rnds

  node3.io.en := io.en
  node3.io.input.valid := node2.io.output.valid
  node3.io.input.rnds  := node2.io.output.rnds

  val fpconv = Module(new IrwinHallFPConverter(spec, rndW))

  fpconv.io.en := io.en
  fpconv.io.input.valid := node3.io.output.valid
  fpconv.io.input.rnd   := node3.io.output.rnds(0)
  assert(node3.io.output.rnds.length == 1)

  io.output.valid := fpconv.io.output.valid
  io.output.rnd   := fpconv.io.output.rnd
}

class IrwinHallTreeSumNode(n: Int, w: Int) extends Module {

  val m = (n+1) / 2

  val io = IO(new Bundle {
    val en = Input(Bool())
    val input = new Bundle {
      val valid = Input(Bool())
      val rnds  = Input(Vec(n, UInt(w.W)))
    }
    val output = new Bundle {
      val valid = Output(Bool())
      val rnds  = Output(Vec(m, UInt((w+1).W)))
    }
  })

  val x = WireDefault(VecInit.fill(2 * m)(0.U(w.W)))
  for(i <- 0 until n) {
    x(i) := io.input.rnds(i)
  }
  val z = WireDefault(VecInit.fill(m)(0.U((w+1).W)))
  for(i <- 0 until m) {
    z(i) := x(i) +& x(i+m)
  }

  val vReg = RegInit(false.B)
  val zReg = RegInit(VecInit.fill(m)(0.U((w+1).W)))
  when(io.en) {
    vReg := io.input.valid
    zReg := z
  }

  io.output.valid := vReg
  io.output.rnds  := zReg
}

class IrwinHallFPConverter(
  spec: RealSpec, rndW: Int
) extends Module {

  val io = IO(new Bundle {
    val en = Input(Bool())
    val input = new Bundle {
      val valid = Input(Bool())
      val rnd   = Input(UInt((rndW+4).W))
    }
    val output = new Bundle {
      val valid = Output(Bool())
      val rnd = Output(UInt(spec.W.W))
    }
  })

  val sum    = io.input.rnd
  val six    = Cat(6.U(4.W), 0.U(rndW.W))
  val zZero  = sum === six
  val zSgn   = Mux(sum >= six, 0.U(1.W),  1.U(1.W))
  val zAbs   = Mux(sum >= six, sum - six, six - sum)

  val zShift = PriorityEncoder(Reverse(zAbs))

  val zShifted0 = (zAbs << zShift)(zAbs.getWidth-1, 0)
  assert(!io.input.valid || zShifted0(zAbs.getWidth-1) === 1.U || zZero)

  val zShifted = zShifted0.tail(1)
  val zRounded = zShifted.head(spec.manW) +& zShifted.head(spec.manW+1)(0)

  val zMoreThan2 = zRounded.head(1)
  val zMan       = zRounded.tail(1)

  val zEx = Mux(zZero, 0.U(spec.exW.W),
                       (spec.exBias+3).U(spec.exW.W) - zShift + zMoreThan2)

  // printf("IrwinHallFPConv: sum      = %b\n", sum)
  // printf("IrwinHallFPConv: 6        = %b, zSgn = %b\n", six, zSgn)
  // printf("IrwinHallFPConv: zAbs     = %b, shift = %d, zEx = %d\n", zAbs, zShift, zEx)
  // printf("IrwinHallFPConv: shifted0 = %b\n", zShifted0)
  // printf("IrwinHallFPConv: shifted  =  %b\n", zShifted)
  // printf("IrwinHallFPConv: rounded  = %b\n", zRounded)
  // printf("IrwinHallFPConv: zMan     =  %b, moreThan2? = %b\n", zMan, zMoreThan2)

  val z = Cat(zSgn, zEx, zMan)
  assert(z.getWidth == spec.W)

  val vReg = RegInit(false.B)
  val zReg = RegInit(0.U(spec.W.W))
  when(io.en) {
    vReg := io.input.valid
    zReg := z
  }

  io.output.valid := vReg
  io.output.rnd   := zReg
}
