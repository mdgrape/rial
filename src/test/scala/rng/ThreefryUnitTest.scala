//

package rial.rng

import java.io.File
import scala.math._
import scala.util.Random
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

import chisel3.util._
import chisel3.iotesters
import chisel3.iotesters.{ChiselFlatSpec, Driver, PeekPokeTester}

import crial._

class ThreefryUnitTester(c: Threefry4_32, ncycle: Int = 1000) extends PeekPokeTester(c) {

  val (r, rotStage) = c.getParam()
  println(f"Threefry4_32 parameters rotation=$r%d stage between rotation=$rotStage%d")

  private val threefry = c

  val key = Array( 0x11111111, 0x22222222, 0, 0 )
  val ctr = Array( 0, 0, 0, 0 )
  poke(threefry.io.key(0), 0x11111111)
  poke(threefry.io.key(1), 0x22222222)
  poke(threefry.io.key(2), 0)
  poke(threefry.io.key(3), 0)


  val rng = crial.threefry4x32_init( ctr, key )

  for (i <- 0 to ncycle-1) {
    //val key = "".U
    //val count = Wire(UInt(128.W)); count = i.U

    poke(threefry.io.count(0), i)
    poke(threefry.io.count(1), 0)
    poke(threefry.io.count(2), 0)
    poke(threefry.io.count(3), 0)

    val r = crial.threefry4x32( rng )

    for (j <- 0 to 3) {
      val x = peek(threefry.io.rand(j))
      val y : Long = if (r(j)<0) {
        r(j) + 0x100000000L
      } else {
        r(j)
      }
      //print(f" $x%08x $y%08x")
      expect(threefry.io.rand(j),y)
    }
    //println()
    step(1)
  }
}

/**
  * This is a trivial example of how to run this Specification
  * From within sbt use:
  * {{{
  * testOnly pipe.ThreefryTester
  * }}}
  *
  * From a terminal shell use:
  * {{{
  * sbt 'testOnly gcd.ThreefryTester'
  * }}}
  */
class ThreefryQTester extends ChiselFlatSpec {
  private val r = 20
  private val rotStage = 0 
  // Disable this until we fix isCommandAvailable to swallow stderr along with stdout
  private val backendNames = if (firrtl.FileUtils.isCommandAvailable(Seq("verilator", "--version"))) {
    Array("firrtl", "verilator")
  }
  else {
    Array("firrtl")
  }
  for (backendName <- backendNames) {
    "Threefry" should s"calculate distance between atoms (with $backendName)" in {
      Driver(() => new Threefry4_32(r, rotStage), backendName) {
        c => new ThreefryUnitTester(c)
      } should be(true)
    }
  }

  "Basic test using Driver.execute" should "be used as an alternative way to run specification" in {
    iotesters.Driver.execute(Array(), () => new Threefry4_32(r, rotStage)) {
      c => new ThreefryUnitTester(c)
    } should be(true)
  }

  if (backendNames.contains("verilator")) {
    "using --backend-name verilator" should "be an alternative way to run using verilator" in {
      iotesters.Driver.execute(Array("--backend-name", "verilator"), () => new Threefry4_32(r, rotStage)) {
        c => new ThreefryUnitTester(c)
      } should be(true)
    }
  }

  "running with --is-verbose" should "show more about what's going on in your tester" in {
    iotesters.Driver.execute(Array("--is-verbose"), () => new Threefry4_32(r, rotStage)) {
      c => new ThreefryUnitTester(c)
    } should be(true)
  }

  /**
    *
    * By default verilator backend produces vcd file, and firrtl and treadle backends do not.
    * Following examples show you how to turn on vcd for firrtl and treadle and how to turn it off for verilator
    */

  "running with --generate-vcd-output on" should "create a vcd file from your test" in {
    iotesters.Driver.execute(
      Array("--generate-vcd-output", "on", "--target-dir", "test_run_dir/make_a_vcd", "--top-name", "make_a_vcd"),
      () => new Threefry4_32(r, rotStage)
    ) {

      c => new ThreefryUnitTester(c)
    } should be(true)

    new File("test_run_dir/make_a_vcd/make_a_vcd.vcd").exists should be(true)
  }

  "running with --generate-vcd-output off" should "not create a vcd file from your test" in {
    iotesters.Driver.execute(
      Array("--generate-vcd-output", "off", "--target-dir", "test_run_dir/make_no_vcd", "--top-name", "make_no_vcd",
        "--backend-name", "verilator"),
      () => new Threefry4_32(r, rotStage)
    ) {

      c => new ThreefryUnitTester(c)
    } should be(true)

    new File("test_run_dir/make_no_vcd/make_a_vcd.vcd").exists should be(false)

  }

}
