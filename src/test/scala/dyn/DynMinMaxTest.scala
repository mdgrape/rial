package rial.tests

import org.scalatest._
import org.scalatest.funspec.AnyFunSpec

import rial.dyn._

import chisel3._
import chisel3.util._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import scala.collection.mutable.Queue
import scala.util.Random
import scala.math._

class DynMinMaxTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))

  // split x into elements of elemW bits. [0] is the LSB side.
  def toElems(x: BigInt, width: Int, elemW: Int, signed: Boolean): Seq[BigInt] = {
    Seq.tabulate(width / elemW)(i => {
      val u = (x >> (elemW * i)) & ((BigInt(1) << elemW) - 1)
      if (signed && u.testBit(elemW-1)) { u - (BigInt(1) << elemW) } else { u }
    })
  }

  def fromElems(xs: Seq[BigInt], elemW: Int): BigInt = {
    val mask = (BigInt(1) << elemW) - 1
    xs.zipWithIndex.foldLeft(BigInt(0))( (acc, x) => acc | ((x._1 & mask) << (elemW * x._2)) )
  }

  def minmaxSoft(
    a: BigInt, b: BigInt, width: Int, elemW: Int, signed: Boolean
  ): (BigInt, BigInt) = {
    val ae = toElems(a, width, elemW, signed)
    val be = toElems(b, width, elemW, signed)
    val mn = ae.zip(be).map{ case(x, y) => x.min(y) }
    val mx = ae.zip(be).map{ case(x, y) => x.max(y) }
    (fromElems(mn, elemW), fromElems(mx, elemW))
  }

  def runtest(laneW: Int, nLanes: Int, signed: Boolean, latency: Int, seed: Int) = {
    it(f"DynMinMax(laneW=${laneW}, nLanes=${nLanes}, signed=${signed}, latency=${latency})") {
      val n = getOption[Int]("n").getOrElse(1000)
      val width = laneW * nLanes
      val elemWs = Seq.tabulate(log2Ceil(nLanes)+1)(m => laneW << m)

      simulate(new DynMinMax(laneW, nLanes, signed, latency)) {
        c => {
          val rng = new Random(seed)

          for(elemW <- elemWs) {
            val mode = DynLane.modeOf(elemW, laneW)
            val q = new Queue[(BigInt, BigInt, BigInt, BigInt)]

            for(i <- 0 until n + latency) {
              // make equal elements appear often enough
              val ain = BigInt(width, rng)
              val bin = if (rng.nextInt(4) == 0) { ain } else { BigInt(width, rng) }

              val (mnref, mxref) = minmaxSoft(ain, bin, width, elemW, signed)
              q += ((ain, bin, mnref, mxref))

              c.io.in.a   .poke(ain.U(width.W))
              c.io.in.b   .poke(bin.U(width.W))
              c.io.in.mode.poke(mode.U)

              val mnout = c.io.out.min.peek().litValue
              val mxout = c.io.out.max.peek().litValue

              c.clock.step(1)

              if(latency <= i) {
                val (x, y, mn, mx) = q.dequeue()
                assert(mnout == mn,
                  f"elemW=${elemW}: min(${x}%x, ${y}%x) = ${mnout}%x, but ref = ${mn}%x")
                assert(mxout == mx,
                  f"elemW=${elemW}: max(${x}%x, ${y}%x) = ${mxout}%x, but ref = ${mx}%x")
              }
            }
          }
        }
      }
    }
  }

  describe("DynMinMax compared to software impl") {
    // int8x4 / int16x2 / int32
    runtest(8, 4, true,  0, 123456789)
    runtest(8, 4, false, 0, 123456789)
    runtest(8, 4, true,  2, 123456789)

    // other configurations
    runtest( 4, 8, true,  0, 123456789)
    runtest( 4, 8, false, 0, 123456789)
    runtest(16, 4, true,  0, 123456789)
    runtest( 8, 1, true,  0, 123456789)
    runtest( 1, 4, true,  0, 123456789)
  }
}
