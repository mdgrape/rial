package rial.rng

import rial.arith.RealSpec

// generate approx-gaussian random number based on Irwin-Hall distribution.

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._

class IrwinHall(
  spec: RealSpec = RealSpec.Float32Spec,
  rndW: Int = 32
) extends Module {

  assert(rndW > spec.W)

  val io = IO(new Bundle {
    val input  = Flipped(Decoupled(Vec(12, UInt(rndW.W))))
    val output = Decoupled(UInt(spec.W.W))
  })

  val node0 = Module(new IrwinHallTreeSumNode(12, rndW))
  val node1 = Module(new IrwinHallTreeSumNode( 6, rndW+1))
  val node2 = Module(new IrwinHallTreeSumNode( 3, rndW+2))
  val node3 = Module(new IrwinHallTreeSumNode( 2, rndW+3))

  io.input.ready := node0.io.input.ready
  node0.io.input.valid := io.input.ready
  node0.io.input.bits  := io.input.bits

  node1.io.input <> node0.io.output
  node2.io.input <> node1.io.output
  node3.io.input <> node2.io.output

  val fpconv = Module(new IrwinHallFPConverter(spec, rndW))

  fpconv.io.input <> node3.io.output
  io.output <> fpconv.io.output
}

class IrwinHallFPConverter(
  spec: RealSpec, rndW: Int
) extends Module {

  val io = IO(new Bundle {
    val input = Flipped(Decoupled(UInt((rndW+4).W)))
    val output = Decoupled(UInt(spec.W.W))
  })

  val sum    = io.input.bits
  val six    = Cat(6.U(4.W), 0.U(rndW.W))
  val zZero  = sum === six
  val zSgn   = Mux(sum >= six, 0.U(1.W),  1.U(1.W))
  val zAbs   = Mux(sum >= six, sum - six, six - sum)

  val zPE    = PriorityEncoder(zAbs)
  val zShift = (zAbs.getWidth-1).U - zPE

  val zShifted0 = (zAbs << zShift)(zAbs.getWidth-1, 0)
  assert(zShifted0(zAbs.getWidth-1) === 1.U || zZero)

  val zShifted = zShifted0.tail(1)
  val zRounded = zShifted.head(spec.manW) +& zShifted.head(spec.manW+1)(0)

  val zMoreThan2 = zRounded.head(1)
  val zMan       = zRounded.tail(1)

  val zEx = Mux(zZero, 0.U(spec.exW.W),
                       (spec.exBias+4).U(spec.exW.W) - zPE + zMoreThan2)

  val z = Cat(zSgn, zEx, zMan)
  assert(z.getWidth == spec.W)

  // 1 cycle to wait
  val zQ = Module(new Queue(UInt(spec.W.W), 1, pipe=true, flow=false))

  io.input.ready  := zQ.io.enq.ready
  zQ.io.enq.valid := io.input.valid
  zQ.io.enq.bits  := z

  io.output <> zQ.io.deq
}

class IrwinHallTreeSumNode(n: Int, w: Int) extends Module {

  val m = (n+1) / 2

  val io = IO(new Bundle {
    val input  = Flipped(Decoupled(Vec(n, UInt(w.W))))
    val output =         Decoupled(Vec(m, UInt((w+1).W)))
  })

  val x = WireDefault(VecInit.fill(2 * m)(0.U(w.W)))
  for(i <- 0 until n) {
    x(i) := io.input.bits(i)
  }
  val z = WireDefault(VecInit.fill(m)(0.U((w+1).W)))
  for(i <- 0 until m) {
    z(i) := x(i) +& x(i+m)
  }

  val zQ = Module(new Queue(Vec(m, UInt((w+1).W)), 1, pipe=true, flow=false))
  io.input.ready  := zQ.io.enq.ready
  zQ.io.enq.valid := io.input.valid
  zQ.io.enq.bits  := z

  io.output <> zQ.io.deq
}
