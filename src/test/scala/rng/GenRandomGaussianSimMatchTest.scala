package rial.tests

import org.scalatest.funspec.AnyFunSpec
import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import rial.rng._
import rial.arith.{MulImpl, RealSpec}
import rial.util.PipelineStageConfig

import scala.util.Random

// Bit-exact match between the GenRandomGaussian hardware and
// GenRandomGaussianSim. A uniform random sweep never reaches the deep tail
// (P(xex0 == 0) = 2^-39), so instead we draw a few samples from EVERY interval
// window -- each main octave (xex), each deep octave (xman leading-zero count),
// the lumped deepest interval, and the central page-0 region -- by forcing the
// bit fields that select that window. This exercises every table entry.
class GenRandomGaussianSimMatchTest extends AnyFunSpec with ChiselSim with EmitVcd {
  addOption(CliOption.int("n", "samples per interval window (default 16)"))

  val spec   = RealSpec.Float32Spec
  val rndW   = 64
  val nOrder = 3
  val lump   = 2
  val manW   = spec.manW
  val mainOctaves = rndW - 2 - manW       // 39 (xex0 width)
  val deepIndiv   = manW - lump           // 21

  private val rng = new Random(20260528)
  private val manMask = (BigInt(1) << manW) - 1
  private def bits(n: Int): BigInt = BigInt(n, rng)
  private def signPage(page: Int): BigInt =
    (BigInt(rng.nextInt(2)) << (rndW - 1)) | (BigInt(page) << (rndW - 2))

  private def pack(sign_page: BigInt, xex0: BigInt, xman: BigInt): BigInt =
    sign_page | (xex0 << manW) | (xman & manMask)

  // page 0 (central): xex0 and xman free.
  private def centralInput(): BigInt = pack(signPage(0), bits(mainOctaves), bits(manW))

  // page 1, main octave k: xex0's lowest set bit at k (bits above free), xman free.
  private def mainInput(k: Int): BigInt = {
    val highW = mainOctaves - k - 1
    val high  = if (highW > 0) bits(highW) else BigInt(0)
    pack(signPage(1), (BigInt(1) << k) | (high << (k + 1)), bits(manW))
  }

  // page 1, deep octave c (xex0 == 0): xman's leading 1 at position manW-1-c,
  // the (manW-1-c) bits below it free -> leading-zero count clz = c.
  private def deepInput(c: Int): BigInt = {
    val lowW = manW - 1 - c
    val low  = if (lowW > 0) bits(lowW) else BigInt(0)
    pack(signPage(1), BigInt(0), (BigInt(1) << (manW - 1 - c)) | low)
  }

  // page 1, lumped deepest interval (xex0 == 0, clz >= deepIndiv): xman < 2^lump.
  private def lumpInput(): BigInt =
    pack(signPage(1), BigInt(0), BigInt(rng.nextInt(1 << lump)))

  def matchTest(polyMulImpl: MulImpl,
    scaleMulImpl: MulImpl = MulImpl.plain,
    splitW: Int = 3, guardBits: Int = 3): Unit = {

    val n = getOption[Int]("n").getOrElse(16)
    simulate(new GenRandomGaussian(rndW, spec, nOrder, PipelineStageConfig.none,
      lump, splitW, guardBits, polyMulImpl, scaleMulImpl)) { c =>
      c.io.rnd.valid.poke(true.B)
      var checked = 0
      val fails = scala.collection.mutable.ArrayBuffer.empty[String]

      def check(window: String, rnd: BigInt): Unit = {
        c.io.rnd.bits.poke(rnd.U(rndW.W))
        c.clock.step(1 + c.getStage) // hold the input and flush the pipeline
        val hw  = c.io.z.bits.peek().litValue
        val sim = GenRandomGaussianSim.genGaussian(rnd, rndW, spec, nOrder,
          lump, splitW, guardBits).value.toBigInt
        checked += 1
        if (hw != sim && fails.length < 20)
          fails += f"$window: rnd=0x$rnd%016x hw=0x$hw%08x sim=0x$sim%08x"
        else if (hw != sim) fails += "..."
      }

      for (_ <- 0 until n) check("central", centralInput())
      for (k <- 0 until mainOctaves; _ <- 0 until n) check(f"main xex=$k", mainInput(k))
      for (cc <- 0 until deepIndiv;  _ <- 0 until n) check(f"deep clz=$cc", deepInput(cc))
      for (_ <- 0 until n) check("lump", lumpInput())

      info(s"checked $checked inputs across ${mainOctaves + deepIndiv + 2} windows")
      assert(fails.isEmpty, "\n  " + fails.mkString("\n  "))
    }
  }

  it("GenRandomGaussian HW == GenRandomGaussianSim across every interval window") {
    matchTest(MulImpl.plain)
  }

  it("... with the products split at 17 too") {
    matchTest(MulImpl.split(17), MulImpl.split(17))
  }

  it("... and with the polynomial's multiply-adds at 17") {
    // the product of the mantissas has no addend to fold, so it stays plain
    matchTest(MulImpl.macc(17))
  }

  it("... and with a multiply-add polynomial over a split mantissa product") {
    matchTest(MulImpl.macc(17), MulImpl.split(17))
  }

  it("... and with the small forms of both") {
    matchTest(MulImpl.maccSmall(17), MulImpl.splitSmall(17))
  }

  it("... and with one polynomial per interval and no guard bits") {
    matchTest(MulImpl.plain, MulImpl.plain, splitW = 0, guardBits = 0)
  }
}
