package rial.tests

import org.scalatest._

import rial.ecc._

import chisel3._
import chiseltest._
import firrtl.AnnotationSeq

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import scala.collection.mutable.{Queue, ArrayBuffer}
import scala.util.Random
import scala.math._

class FletchersChecksumTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Fletcher's checksum compared to software impl"

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

  def runtest(dataW: Int, latency: Int, n: Int, seed: Int, anno: AnnotationSeq) = {
    it should f"Fletcher16(dataW=${dataW}, latency=${latency})" in {
      test(new Fletcher(16, dataW, latency)).withAnnotations(anno) {
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

  runtest( 32, 0, 1000, 123456789, Seq())
  runtest( 32, 1, 1000, 123456789, Seq())
  runtest( 32, 2, 1000, 123456789, Seq())

  runtest( 64, 0, 1000, 123456789, Seq())
  runtest( 64, 1, 1000, 123456789, Seq())
  runtest( 64, 2, 1000, 123456789, Seq())

  runtest(128, 0, 1000, 123456789, Seq())
  runtest(128, 1, 1000, 123456789, Seq())
  runtest(128, 2, 1000, 123456789, Seq())
}
