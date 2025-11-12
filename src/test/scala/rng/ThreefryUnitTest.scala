package rial.tests

import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chisel3.util._
import chiseltest._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig
import rial.util.DebugControlMaster
import rial.rng._

import java.io.File

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

import crial._

class ThreefryTest extends AnyFlatSpec
  with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Threefry"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  def runTest(r: Int, rotStage: Int) = {
    test(new Threefry4_32(r, rotStage)) {
      c => {
        val (r, rotStage) = c.getParam
        println(f"Threefry4_32 parameters rotation=$r%d stage between rotation=$rotStage%d, nStage = ${c.nStage}")

        val threefry = c
        threefry.io.en.poke(true.B)

        val key = Array( 0x11111111, 0x22222222, 0, 0 )
        val ctr = Array( 0, 0, 0, 0 )
        threefry.io.input.valid.poke(true.B)
        threefry.io.input.key(0).poke(0x11111111.U)
        threefry.io.input.key(1).poke(0x22222222.U)
        threefry.io.input.key(2).poke(0.U)
        threefry.io.input.key(3).poke(0.U)

        val rng = crial.threefry4x32_init( ctr, key )

        for (i <- 0 to n-1) {
          //val key = "".U
          //val count = Wire(UInt(128.W)); count = i.U

          threefry.io.input.count(0).poke(i.U)
          threefry.io.input.count(1).poke(0.U)
          threefry.io.input.count(2).poke(0.U)
          threefry.io.input.count(3).poke(0.U)

          threefry.io.output.valid.expect((i >= threefry.nStage).B)

          if(i >= threefry.nStage) {
            val r = crial.threefry4x32( rng )
            for (j <- 0 to 3) {
              val y : Long = if (r(j)<0) {
                r(j) + 0x100000000L
              } else {
                r(j)
              }
              val x = threefry.io.output.rand(j).peek().litValue.toLong
              if(x != y) {
                println(f"actual=${x}%08x, expected=${y}%08x")
              }
              threefry.io.output.rand(j).expect(y.U)
            }
          }
          c.clock.step(1)
        }
      }
    }
  }

  it should "Threefry should be equal to crial C impl" in {
    runTest(20, 0)
    runTest(20, 20)
    runTest(20, 10)
    runTest(20, 5)
    runTest(20, 4)
    runTest(20, 2)
  }
}

class ThreefrySimTest extends AnyFlatSpec
  with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Threefry software impl"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  def runTest(r: Int, rotStage: Int) = {
    test(new Threefry4_32(r, rotStage)) {
      c => {
        val (r, rotStage) = c.getParam
        println(f"Threefry4_32 parameters rotation=$r%d stage between rotation=$rotStage%d, nStage = ${c.nStage}")

        val threefry = c
        threefry.io.en.poke(true.B)

        val key = Array( 0x11111111, 0x22222222, 0, 0 )
        val ctr = Array( 0, 0, 0, 0 )
        threefry.io.input.valid.poke(true.B)
        threefry.io.input.key(0).poke(0x11111111.U)
        threefry.io.input.key(1).poke(0x22222222.U)
        threefry.io.input.key(2).poke(0.U)
        threefry.io.input.key(3).poke(0.U)

        val sim = new ThreefrySim(r, Seq(0x11111111L, 0x22222222L, 0L, 0L))

        for (i <- 0 to n-1) {
          //val key = "".U
          //val count = Wire(UInt(128.W)); count = i.U

          threefry.io.input.count(0).poke(i.U)
          threefry.io.input.count(1).poke(0.U)
          threefry.io.input.count(2).poke(0.U)
          threefry.io.input.count(3).poke(0.U)

          threefry.io.output.valid.expect((i >= threefry.nStage).B)

          if(i >= threefry.nStage) {
            val r = sim.gen(Seq(i-threefry.nStage, 0, 0, 0))
            for (j <- 0 to 3) {
              val y = r(j)
              val x = threefry.io.output.rand(j).peek().litValue.toLong
              if(x != y) {
                println(f"actual=${x}%08x, expected=${y}%08x")
                c.clock.step(1)
              }
              threefry.io.output.rand(j).expect(y.U)
            }
          }
          c.clock.step(1)
        }
      }
    }
  }

  it should "Threefry should be equal to software impl" in {
    runTest(20, 0)
    // runTest(20, 20)
    // runTest(20, 10)
    // runTest(20, 5)
    // runTest(20, 4)
    // runTest(20, 2)
  }
}

// test of wrapper
class ThreefryGeneratorTest extends AnyFlatSpec
  with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test ThreefryGenerator"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  def runTest(r: Int, rotStage: Int) = {
    test(new Threefry4x32Generator(r, rotStage)) {
      c => {
        val (r, rotStage) = c.getParam
        println(f"Threefry4x32 parameters rotation=${r}%d stage between rotation=${rotStage}%d, nStage = ${c.nStage}")

        val threefry = c

        val key = Array( 0x11111111, 0x22222222, 0, 0 )
        val ctr = Array( 0, 0, 0, 0 )

        threefry.io.initialized.expect(false.B)

        threefry.io.init.en.poke(true.B)
        threefry.io.init.key  (0).poke(key(0).U)
        threefry.io.init.key  (1).poke(key(1).U)
        threefry.io.init.key  (2).poke(key(2).U)
        threefry.io.init.key  (3).poke(key(3).U)
        threefry.io.init.count(0).poke(ctr(0).U)
        threefry.io.init.count(1).poke(ctr(1).U)
        threefry.io.init.count(2).poke(ctr(2).U)
        threefry.io.init.count(3).poke(ctr(3).U)

        threefry.clock.step(1) // --------------------------------------------

        threefry.io.initialized.expect(true.B)

        threefry.io.init.en.poke(false.B)
        threefry.io.init.key(0).poke(0.U)
        threefry.io.init.key(1).poke(0.U)
        threefry.io.init.key(2).poke(0.U)
        threefry.io.init.key(3).poke(0.U)
        threefry.io.init.count(0).poke(0.U)
        threefry.io.init.count(1).poke(0.U)
        threefry.io.init.count(2).poke(0.U)
        threefry.io.init.count(3).poke(0.U)

        val rng = crial.threefry4x32_init( ctr, key )

        threefry.io.rand.ready.poke(true.B)
        for (i <- 0 to n-1) {

          threefry.io.rand.valid.expect( (i >= (threefry.nStage)).B )
          if(threefry.io.rand.valid.peek().litValue == 1) {

            val r = crial.threefry4x32( rng )
            for (j <- 0 to 3) {
              val y : Long = if (r(j)<0) {
                r(j) + 0x100000000L
              } else {
                r(j)
              }
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

  it should "Threefry should be equal to crial C impl" in {
    runTest(20, 0)  // latency  0 + 1 (queue)
    runTest(20, 20) // latency  1 + 1
    runTest(20, 10) // latency  2 + 1
    runTest(20, 5)  // latency  4 + 1
    runTest(20, 4)  // latency  5 + 1
    runTest(20, 2)  // latency 10 + 1
  }
}


// /**
//   * This is a trivial example of how to run this Specification
//   * From within sbt use:
//   * {{{
//   * testOnly pipe.ThreefryTester
//   * }}}
//   *
//   * From a terminal shell use:
//   * {{{
//   * sbt 'testOnly gcd.ThreefryTester'
//   * }}}
//   */
// class ThreefryQTester extends ChiselFlatSpec {
//   private val r = 20
//   private val rotStage = 0 
//   // Disable this until we fix isCommandAvailable to swallow stderr along with stdout
//   private val backendNames = if (firrtl.FileUtils.isCommandAvailable(Seq("verilator", "--version"))) {
//     Array("firrtl", "verilator")
//   }
//   else {
//     Array("firrtl")
//   }
//   for (backendName <- backendNames) {
//     "Threefry" should s"calculate distance between atoms (with $backendName)" in {
//       Driver(() => new Threefry4_32(r, rotStage), backendName) {
//         c => new ThreefryUnitTester(c)
//       } should be(true)
//     }
//   }
//
//   "Basic test using Driver.execute" should "be used as an alternative way to run specification" in {
//     iotesters.Driver.execute(Array(), () => new Threefry4_32(r, rotStage)) {
//       c => new ThreefryUnitTester(c)
//     } should be(true)
//   }
//
//   if (backendNames.contains("verilator")) {
//     "using --backend-name verilator" should "be an alternative way to run using verilator" in {
//       iotesters.Driver.execute(Array("--backend-name", "verilator"), () => new Threefry4_32(r, rotStage)) {
//         c => new ThreefryUnitTester(c)
//       } should be(true)
//     }
//   }
//
//   "running with --is-verbose" should "show more about what's going on in your tester" in {
//     iotesters.Driver.execute(Array("--is-verbose"), () => new Threefry4_32(r, rotStage)) {
//       c => new ThreefryUnitTester(c)
//     } should be(true)
//   }
//
//   /**
//     *
//     * By default verilator backend produces vcd file, and firrtl and treadle backends do not.
//     * Following examples show you how to turn on vcd for firrtl and treadle and how to turn it off for verilator
//     */
//
//   "running with --generate-vcd-output on" should "create a vcd file from your test" in {
//     iotesters.Driver.execute(
//       Array("--generate-vcd-output", "on", "--target-dir", "test_run_dir/make_a_vcd", "--top-name", "make_a_vcd"),
//       () => new Threefry4_32(r, rotStage)
//     ) {
//
//       c => new ThreefryUnitTester(c)
//     } should be(true)
//
//     new File("test_run_dir/make_a_vcd/make_a_vcd.vcd").exists should be(true)
//   }
//
//   "running with --generate-vcd-output off" should "not create a vcd file from your test" in {
//     iotesters.Driver.execute(
//       Array("--generate-vcd-output", "off", "--target-dir", "test_run_dir/make_no_vcd", "--top-name", "make_no_vcd",
//         "--backend-name", "verilator"),
//       () => new Threefry4_32(r, rotStage)
//     ) {
//
//       c => new ThreefryUnitTester(c)
//     } should be(true)
//
//     new File("test_run_dir/make_no_vcd/make_a_vcd.vcd").exists should be(false)
//
//   }
//
// }
