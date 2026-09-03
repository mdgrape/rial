package rial.tests

import org.scalatest._
import org.scalatest.funspec.AnyFunSpec

import rial.dyn._
import rial.testUtil.RialChiselSim

import chisel3._
import chisel3.util._
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import scala.collection.mutable.Queue
import scala.util.Random
import scala.math._

class DynCompareTest extends AnyFunSpec with RialChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))

  // split x into elements of elemW bits. [0] is the LSB side.
  def toElems(x: BigInt, width: Int, elemW: Int, signed: Boolean): Seq[BigInt] = {
    Seq.tabulate(width / elemW)(i => {
      val u = (x >> (elemW * i)) & ((BigInt(1) << elemW) - 1)
      if (signed && u.testBit(elemW-1)) { u - (BigInt(1) << elemW) } else { u }
    })
  }

  // expected per-lane flags. the result of an element is broadcast to all the
  // lanes it occupies.
  def compareSoft(
    a: BigInt, b: BigInt, width: Int, laneW: Int, elemW: Int, signed: Boolean
  ): Seq[(Boolean, Boolean, Boolean)] = {
    val ae = toElems(a, width, elemW, signed)
    val be = toElems(b, width, elemW, signed)
    Seq.tabulate(width / laneW)(i => {
      val e = i / (elemW / laneW)
      (ae(e) == be(e), ae(e) < be(e), ae(e) > be(e))
    })
  }

  def runtest(laneW: Int, nLanes: Int, signed: Boolean, latency: Int, seed: Int) = {
    it(f"DynCompare(laneW=${laneW}, nLanes=${nLanes}, signed=${signed}, latency=${latency})") {
      val n = getOption[Int]("n").getOrElse(1000)
      val width = laneW * nLanes
      val elemWs = Seq.tabulate(log2Ceil(nLanes)+1)(m => laneW << m)

      simulate(new DynCompare(laneW, nLanes, signed, latency)) {
        c => {
          val rng = new Random(seed)

          for(elemW <- elemWs) {
            val mode = DynLane.modeOf(elemW, laneW)
            val q = new Queue[(BigInt, BigInt, Seq[(Boolean, Boolean, Boolean)])]

            for(i <- 0 until n + latency) {
              // make equal elements appear often enough
              val ain = BigInt(width, rng)
              val bin = if (rng.nextInt(4) == 0) { ain } else { BigInt(width, rng) }

              q += ((ain, bin, compareSoft(ain, bin, width, laneW, elemW, signed)))

              c.io.in.a   .poke(ain.U(width.W))
              c.io.in.b   .poke(bin.U(width.W))
              c.io.in.mode.poke(mode.U)

              val out = Seq.tabulate(nLanes)(l => (
                c.io.out.eq(l).peek().litValue == 1,
                c.io.out.lt(l).peek().litValue == 1,
                c.io.out.gt(l).peek().litValue == 1
              ))

              c.clock.step(1)

              if(latency <= i) {
                val (x, y, ref) = q.dequeue()
                for(l <- 0 until nLanes) {
                  assert(out(l) == ref(l),
                    f"elemW=${elemW}, lane=${l}: cmp(${x}%x, ${y}%x) = ${out(l)}, but ref = ${ref(l)}")
                }
              }
            }
          }
        }
      }
    }
  }

  describe("DynCompare compared to software impl") {
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
