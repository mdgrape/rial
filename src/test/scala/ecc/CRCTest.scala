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

class CRC16Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "CRC-16(g=0x18005) result compared to software impl"

  // use the simplest approach because it is the reference
  def crc16Soft(data: Seq[Int]): Int = {
    val gen = 0x18005

    val input = data.foldLeft(BigInt(0))( (u, l) => {
      (u << 16) | (l & 0xFFFF)
    })

    var x = input << 16 // add 16x zeroes at the LSB side to save the remnant
    var g = BigInt(gen)

    // println(f"x = ${x.toString(2)}")
    for(i <- 16*(data.length+1) to 16 by -1) {
      val msb = x & (BigInt(1) << i)
      if(msb != 0) {
        x = x ^ (g << (i-16))
        // println(f"g = ${(g << (i-16)).toString(2)}")
        // println(f"x = ${x.toString(2)}")
      }
    }
    x.toInt
  }

  def runtest(dlen: Int, n: Int, seed: Int, anno: AnnotationSeq) = {
    it should f"CRC16(g=0x18005, w=16)" in {
      test(new CRCcore(0x18005, 16)).withAnnotations(anno) {
        c => {
          val rng = new Random(seed)

          for(testcase <- 0 until n) {
            val xin = Seq.tabulate(dlen)(i => rng.nextInt() & 0xFFFF)

            var rem = 0
            c.io.ri.poke(rem.U)

            for(i <- 0 until dlen) {

              c.io.x .poke(xin(i).U(16.W))
              c.io.ri.poke(rem.U(16.W))
              rem = c.io.ro.peek().litValue.toInt

              c.clock.step(1)
            }

            val ref = crc16Soft(xin)
            assert(rem == ref, f"rem=0x${rem}%x, ref=0x${ref}%x, g=0x18005, xin=${xin.foldLeft(BigInt(0))( (u, l) => {(u << 16) | (l & 0xFFFF)})}%x")
            println(f"rem=0x${rem}%x, ref=0x${ref}%x, g=0x18005, xin=${xin.foldLeft(BigInt(0))( (u, l) => {(u << 16) | (l & 0xFFFF)})}%x")

            c.io.x .poke(rem.U(16.W))
            c.io.ri.poke(rem.U(16.W))

            val zero = c.io.ro.peek().litValue.toInt
            assert(zero == 0)
          }
        }
      }
    }
  }
  runtest(4, 1000, 123456789, Seq())
}
