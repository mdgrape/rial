package rial.tests

import org.scalatest._
import org.scalatest.funspec.AnyFunSpec

import rial.ecc._

import chisel3._
import chisel3.util._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import scala.collection.mutable.{Queue, ArrayBuffer}
import scala.util.Random
import scala.math._

class FletchersChecksumTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))

  // Ai = Ai-1 + di
  // Bi = Bi-1 + Ai
  def fletcher16Soft(w: Int, data: BigInt): Int = {

    val ds = Seq.tabulate((w+7) / 8)(i => {
      ((data >> (8 * i)) & 0x00FF).toInt
    })

    var a = 0
    var b = 0
    for(d <- ds) {
      a = a + d
      b = b + a
    }
    val amod = a % 255
    val bmod = b % 255

    (bmod << 8) + amod
  }

  def runtest(dataW: Int, latency: Int, seed: Int) = {
    it(f"Fletcher16(dataW=${dataW}, latency=${latency})") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new Fletcher(16, dataW, latency)) {
        c => {
          val rng = new Random(seed)

          val q = new Queue[(BigInt, Int)]
          for(i <- 0 until n + latency) {

            val xin = BigInt(dataW, rng).abs
            val zref = fletcher16Soft(dataW, xin)
            q += ((xin, zref))

            c.io.in.valid.poke(true.B)
            c.io.in.bits .poke(xin.U(dataW.W))
            val zvalid = c.io.out.valid.peek().litValue == 1
            val zout   = c.io.out.bits.peek().litValue

            // println(f"xin=${xin}%x, zout=${zout}%x, zref=${zref}%x")

            c.clock.step(1)

            if(latency <= i) {
              val (x, z) = q.dequeue()
              assert(zvalid)
              assert(z == zout, f"checksum(${x}%x) = ${z}%x, but checksumSoft(${x}%x) = ${zout}%x")
            }
          }
        }
      }
    }
  }

  describe("Fletcher's checksum compared to software impl") {
    runtest( 32, 0, 123456789)
    runtest( 32, 1, 123456789)
    runtest( 32, 2, 123456789)

    runtest( 64, 0, 123456789)
    runtest( 64, 1, 123456789)
    runtest( 64, 2, 123456789)

    runtest(128, 0, 123456789)
    runtest(128, 1, 123456789)
    runtest(128, 2, 123456789)
  }
}
