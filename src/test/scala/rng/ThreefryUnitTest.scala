package rial.tests

import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chisel3.util._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption
import org.scalatest.funspec.AnyFunSpec

import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig
import rial.rng._

import scala.util.Random
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

//
// Known-answer vectors anchoring [[ThreefrySim]] to the reference implementation.
//
class ThreefryKatTest extends AnyFunSpec {

  // (count, key, expected)
  val kat: Seq[(Seq[Long], Seq[Long], Seq[Long])] = Seq(
    (Seq(0x00000000L, 0x00000000L, 0x00000000L, 0x00000000L),
     Seq(0x00000000L, 0x00000000L, 0x00000000L, 0x00000000L),
     Seq(0x9c6ca96aL, 0xe17eae66L, 0xfc10ecd4L, 0x5256a7d8L)),
    (Seq(0xffffffffL, 0xffffffffL, 0xffffffffL, 0xffffffffL),
     Seq(0xffffffffL, 0xffffffffL, 0xffffffffL, 0xffffffffL),
     Seq(0x2a881696L, 0x57012287L, 0xf6c7446eL, 0xa16a6732L)),
    (Seq(0x243f6a88L, 0x85a308d3L, 0x13198a2eL, 0x03707344L),
     Seq(0xa4093822L, 0x299f31d0L, 0x082efa98L, 0xec4e6c89L),
     Seq(0x59cd1dbbL, 0xb8879579L, 0x86b5d00cL, 0xac8b6d84L)),
    (Seq(0x00000000L, 0x00000000L, 0x00000000L, 0x00000000L),
     Seq(0x11111111L, 0x22222222L, 0x00000000L, 0x00000000L),
     Seq(0x81831df8L, 0xd92ae0d7L, 0x30ef273eL, 0x63d4e43bL)),
    (Seq(0x00000001L, 0x00000000L, 0x00000000L, 0x00000000L),
     Seq(0x11111111L, 0x22222222L, 0x00000000L, 0x00000000L),
     Seq(0xf5a9e3deL, 0x3a23e829L, 0x0dffafa2L, 0x68bc1ac8L)),
    (Seq(0xffffffffL, 0x00000000L, 0x00000000L, 0x00000000L),
     Seq(0x11111111L, 0x22222222L, 0x00000000L, 0x00000000L),
     Seq(0xae6b3e91L, 0xa7a9034eL, 0x36bac8a8L, 0x1d364609L)),
    (Seq(0x00000000L, 0x00000000L, 0x00000000L, 0x00000001L),
     Seq(0x0badc0deL, 0xdeadbeefL, 0x00000001L, 0x80000000L),
     Seq(0x939eba18L, 0x3dea0587L, 0x550e4e3aL, 0x94773692L)),
    (Seq(0x89abcdefL, 0x01234567L, 0xfedcba98L, 0x76543210L),
     Seq(0x00000001L, 0x00000000L, 0x00000000L, 0x00000000L),
     Seq(0x55902d76L, 0xf4968568L, 0x8c3da80bL, 0xa00f6461L))
  )

  it("ThreefrySim should reproduce the Random123 known-answer vectors") {
    for ((count, key, expected) <- kat) {
      val got = new ThreefrySim(20, key).gen(count)
      assert(got == expected,
        f"count=${count.map(c => f"$c%08x").mkString(",")} " +
        f"key=${key.map(k => f"$k%08x").mkString(",")}: " +
        f"got=${got.map(g => f"$g%08x").mkString(",")} " +
        f"expected=${expected.map(e => f"$e%08x").mkString(",")}")
    }
  }

  it("BitRngSim.threefry4x32 should agree with ThreefrySim (jint word interface)") {
    for ((count, key, expected) <- kat) {
      val got = BitRngSim.threefry4x32(
        count.map(_.toInt).toArray, key.map(_.toInt).toArray, 20)
      assert(got.toSeq.map(w => w.toLong & 0xFFFFFFFFL) == expected)
    }
  }
}

//
// Threefry4_32 vs ThreefrySim.
//
class ThreefryTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000)

  val nKey = 3

  def runTest(r: Int, rotStage: Int) = {
    simulate(new Threefry4_32(r, rotStage)) {
      c => {
        val (r, rotStage) = c.getParam
        println(f"Threefry4_32 parameters rotation=$r%d stage between rotation=$rotStage%d, nStage = ${c.nStage}")

        val threefry = c
        val nStage   = c.nStage
        val rng      = new Random(19660809 + rotStage)

        // counts: edge cases first, then random over the full 128-bit range.
        def count(i: Int): Seq[Long] = i match {
          case 0 => Seq(0x00000000L, 0x00000000L, 0x00000000L, 0x00000000L)
          case 1 => Seq(0xffffffffL, 0xffffffffL, 0xffffffffL, 0xffffffffL)
          case 2 => Seq(0x80000000L, 0x00000000L, 0x00000000L, 0x80000000L)
          case _ => Seq.fill(4)(rng.nextInt().toLong & 0xffffffffL)
        }

        val keys = Seq(Seq(0x11111111L, 0x22222222L, 0L, 0L)) ++
          Seq.tabulate(nKey - 1)(_ => Seq.fill(4)(rng.nextInt().toLong & 0xffffffffL))

        threefry.io.en.poke(true.B)
        threefry.io.input.valid.poke(true.B)

        for ((key, ki) <- keys.zipWithIndex) {
          for (j <- 0 to 3) threefry.io.input.key(j).poke(key(j).U(32.W))

          val sim = new ThreefrySim(r, key)
          val q   = new Queue[Seq[Long]]

          for (i <- 0 to n + nStage - 1) {
            val cnt = count(i)
            q += cnt
            for (j <- 0 to 3) { threefry.io.input.count(j).poke(cnt(j).U(32.W)) }

            // valid only ramps up on the first key; afterwards the pipeline
            // stays full across the key change.
            if (ki == 0) {
              threefry.io.output.valid.expect((i >= nStage).B)
            }

            if (i >= nStage) {
              val expected = sim.gen(q.dequeue())
              for (j <- 0 to 3) {
                val y = expected(j)
                val x = threefry.io.output.rand(j).peek().litValue.toLong
                if (x != y) {
                  println(f"actual=${x}%08x, expected=${y}%08x")
                }
                threefry.io.output.rand(j).expect(y.U)
              }
            }
            c.clock.step(1)
          }
          q.clear()
        }
      }
    }
  }

  it("Threefry should be equal to software impl") {
    runTest(20, 0)
    runTest(20, 20)
    runTest(20, 10)
    runTest(20, 5)
    runTest(20, 4)
    runTest(20, 2)
  }
}

// test of wrapper
class ThreefryGeneratorTest extends AnyFunSpec
  with ChiselSim with EmitVcd {


  addOption(CliOption.int("n", "number of test cases (default 1000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000)

  /** 128-bit little-endian counter increment, matching Threefry4x32Generator. */
  def increment(count: Seq[Long]): Seq[Long] = {
    val next = count.toArray
    var j    = 0
    var carry = true
    while (carry && j < 4) {
      next(j) = (next(j) + 1) & 0xffffffffL
      carry   = next(j) == 0
      j      += 1
    }
    next.toSeq
  }

  def runTest(r: Int, rotStage: Int) = {
    simulate(new Threefry4x32Generator(r, rotStage)) {
      c => {
        val (r, rotStage) = c.getParam
        println(f"Threefry4x32 parameters rotation=${r}%d stage between rotation=${rotStage}%d, nStage = ${c.nStage}")

        val threefry = c

        val key = Seq( 0x11111111L, 0x22222222L, 0L, 0L )
        val ctr = Seq( 0L, 0L, 0L, 0L )

        threefry.io.initialized.expect(false.B)

        threefry.io.init.en.poke(true.B)
        for (j <- 0 to 3) { threefry.io.init.key  (j).poke(key(j).U(32.W)) }
        for (j <- 0 to 3) { threefry.io.init.count(j).poke(ctr(j).U(32.W)) }

        threefry.clock.step(1) // --------------------------------------------

        threefry.io.initialized.expect(true.B)

        threefry.io.init.en.poke(false.B)
        for (j <- 0 to 3) { threefry.io.init.key  (j).poke(0.U) }
        for (j <- 0 to 3) { threefry.io.init.count(j).poke(0.U) }

        val sim   = new ThreefrySim(r, key)
        var count = ctr

        threefry.io.rand.ready.poke(true.B)
        for (i <- 0 to n-1) {

          threefry.io.rand.valid.expect( (i >= (threefry.nStage)).B )
          if(threefry.io.rand.valid.peek().litValue == 1) {

            val expected = sim.gen(count)
            count = increment(count)
            for (j <- 0 to 3) {
              val y = expected(j)
              val x = threefry.io.rand.bits(j).peek().litValue.toLong
              if(x != y) {
                println(f"actual=${x}%08x, expected=${y}%08x")
              }
              assert(x == y)
              // threefry.io.rand.bits(j).expect(y.U)
            }
          }
          c.clock.step(1)
        }
      }
    }
  }

  it("Threefry should be equal to software impl") {
    runTest(20, 0)  // latency  0 + 1 (queue)
    runTest(20, 20) // latency  1 + 1
    runTest(20, 10) // latency  2 + 1
    runTest(20, 5)  // latency  4 + 1
    runTest(20, 4)  // latency  5 + 1
    runTest(20, 2)  // latency 10 + 1
  }
}
