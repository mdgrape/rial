package rial.tests

import org.scalatest.funspec.AnyFunSpec
import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import rial.rng._
import rial.arith.MulImpl
import rial.fixed.FixedSpec
import rial.util.PipelineStageConfig

import scala.collection.mutable.ArrayBuffer

// Bit-exact match between the GenRandomGaussianFixed hardware and
// GenRandomGaussianFixedSim. Random words rarely reach the deep octaves, so
// every octave is also forced once per sign by setting its scan bit and
// clearing the ones below, next to the words with no scan bit set and the
// all-zero and all-one words.
class GenRandomGaussianFixedSimMatchTest extends AnyFunSpec with ChiselSim with EmitVcd {
  addOption(CliOption.int("n", "random words per configuration (default 3000)"))

  private val sr = new java.util.SplittableRandom(20261027L)
  private def bits(n: Int): BigInt = {
    (0 until n by 32).foldLeft(BigInt(0))((acc, i) => {
      acc | (BigInt(sr.nextInt() & 0xffffffffL) << i)
    }) & ((BigInt(1) << n) - 1)
  }

  private def words(rndW: Int, zSpec: FixedSpec): Seq[BigInt] = {
    val positionW = GenRandomGaussianFixed.positionWidth(zSpec)
    val scanW     = rndW - 2 - positionW
    val nOctaves  = GenRandomGaussianFixed.octaves(rndW, zSpec)
    val n         = getOption[Int]("n").getOrElse(3000)
    val all       = (BigInt(1) << rndW) - 1

    def word(sgn: Int, page: Int, scan: BigInt, pos: BigInt): BigInt = {
      (BigInt(sgn) << (rndW - 1)) | (BigInt(page) << (rndW - 2)) |
        (scan << positionW) | pos
    }
    val edges = ArrayBuffer[BigInt](BigInt(0), all)
    for (sgn <- 0 to 1) {
      for (pos <- Seq(BigInt(0), (BigInt(1) << positionW) - 1, bits(positionW))) {
        edges += word(sgn, 0, bits(scanW), pos)
        edges += word(sgn, 1, BigInt(0), pos)
      }
      // octave i: bit i of the scan set, the bits below it clear, the ones
      // above, the ignored ones included, random
      for (i <- 0 until nOctaves) {
        val above = bits(scanW - i - 1) << (i + 1)
        edges += word(sgn, 1, above | (BigInt(1) << i), bits(positionW))
      }
      // the ignored bits of the scan field do not reach the last row
      if (scanW > nOctaves) {
        edges += word(sgn, 1, bits(scanW - nOctaves) << nOctaves, bits(positionW))
      }
    }
    edges.toSeq ++ Seq.fill(n)(bits(rndW))
  }

  def matchTest(rndW: Int, zSpec: FixedSpec, stage: PipelineStageConfig,
    polyMulImpl: MulImpl): Unit = {

    it(f"Q${zSpec.intW}.${zSpec.fracW} rndW=$rndW stage=(${stage.getString}) " +
       f"poly=$polyMulImpl matches the model") {
      val ws = words(rndW, zSpec)
      simulate(new GenRandomGaussianFixed(rndW, zSpec, stage = stage,
          polyMulImpl = polyMulImpl)) { c =>
        val lat = c.getStage

        // the latency: a single valid word comes out after getStage cycles
        c.io.rnd.valid.poke(false.B)
        c.clock.step(lat + 1)
        c.io.rnd.valid.poke(true.B)
        c.io.rnd.bits.poke(0.U)
        var seen = -1
        for (k <- 0 to lat + 1) {
          if (seen < 0 && c.io.z.valid.peek().litValue == 1) { seen = k }
          c.clock.step(1)
          c.io.rnd.valid.poke(false.B)
        }
        assert(seen == lat, s"z.valid after $seen cycles, getStage = $lat")

        // one word per cycle; the output of word i appears at cycle i + lat
        val fails = ArrayBuffer.empty[String]
        c.io.rnd.valid.poke(true.B)
        for (k <- 0 until ws.length + lat) {
          if (k < ws.length) { c.io.rnd.bits.poke(ws(k).U(rndW.W)) }
          if (k == ws.length) { c.io.rnd.valid.poke(false.B) }
          if (k >= lat) {
            val w   = ws(k - lat)
            assert(c.io.z.valid.peek().litValue == 1)
            val hw  = c.io.z.bits.peek().litValue
            val sim = GenRandomGaussianFixedSim.genGaussian(w, rndW, zSpec)
            if (hw != sim && fails.length < 20) {
              fails += f"rnd=0x$w%x hw=0x$hw%x sim=0x$sim%x"
            }
          }
          c.clock.step(1)
        }
        info(s"checked ${ws.length} words, latency $lat")
        assert(fails.isEmpty, "\n  " + fails.mkString("\n  "))
      }
    }
  }

  matchTest(64, FixedSpec(16, 12), PipelineStageConfig.none, MulImpl.plain)
  matchTest(48, FixedSpec(24, 16), PipelineStageConfig.default(4), MulImpl.macc(17))
  matchTest(32, FixedSpec(16, 13), PipelineStageConfig.specify(
    "in" -> 1, "decode" -> 1, "table" -> 1, "poly.prod1" -> 1, "poly.step2" -> 1,
    "poly.out" -> 1, "round" -> 1, "out" -> 2), MulImpl.plain)
}
