package rial.tests

import org.scalatest.funspec.AnyFunSpec

import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import circt.stage.ChiselStage

import rial.math.MathFuncConfig
import rial.table.{PolynomialBitWidth, PolynomialEval}
import rial.util.{PipelineStageConfig, PipelineStages}

import scala.collection.mutable.ArrayBuffer
import scala.util.Random

/** Tests of the pipelined polynomial evaluator.
 *
 * The pipeline registers must not change what the module computes, so every
 * pipelined configuration is compared, bit for bit, against the combinational
 * one, with the results shifted by the latency the configuration resolves to.
 */
class PolynomialEvalTest extends AnyFunSpec with ChiselSim {

  /** The polynomial of the standard FP32 math function set (order 2). */
  val bitW2: PolynomialBitWidth = MathFuncConfig.standardFP32.bitW

  /** A third order variant of the same shape. */
  val bitW3: PolynomialBitWidth = PolynomialBitWidth(
    nOrder   = 3,
    addrW    = bitW2.addrW,
    inBitW   = bitW2.inBitW,
    outBitW  = bitW2.outBitW,
    coefBits = Seq(27, 19, 12, 8),
    calcBits = Seq(27, 19, 12, 8))

  /** number of random vectors per configuration */
  val nVec = 200

  /** Random coefficients and positions.
   *
   * Each coefficient is drawn two bits short of its field so that no Horner
   * intermediate overflows its `calcBits`, which the module asserts.
   */
  def randVectors(bitW: PolynomialBitWidth, n: Int, r: Random): Seq[(Seq[BigInt], BigInt)] = {
    Seq.fill(n) {
      val cs = bitW.coefBits.map { w =>
        BigInt(w - 1, r) - (BigInt(1) << (w - 2))
      }
      val x = BigInt(bitW.inBitW, r) - (BigInt(1) << (bitW.inBitW - 1))
      (cs, x)
    }
  }

  /** Runs the vectors through one configuration and returns the results.
   *
   * The `i`-th vector is applied at cycle `i`; the result read back `latency`
   * cycles later is the one that belongs to it.
   */
  def run(bitW: PolynomialBitWidth, stage: PipelineStageConfig, clamp: Boolean,
    vectors: Seq[(Seq[BigInt], BigInt)]): Seq[BigInt] = {

    val res = ArrayBuffer.empty[BigInt]
    simulate(new PolynomialEval(bitW, stage, clamp)) { c =>
      val latency = c.getStage
      assert(latency == stage.resolve(PolynomialEval.stagePoints(bitW.nOrder)).total)
      c.io.in.valid.poke(true.B)
      for (i <- 0 until vectors.length + latency) {
        val (cs, x) = vectors(scala.math.min(i, vectors.length - 1))
        for (k <- cs.indices) {
          c.io.in.bits.coeffs(k).poke(cs(k).S(bitW.coefBits(k).W))
        }
        if (bitW.nOrder != 0) {
          c.io.in.bits.x.get.poke(x.S(bitW.inBitW.W))
        }
        if (i >= latency) {
          assert(c.io.out.valid.peek().litToBoolean,
            s"output is not valid at cycle $i (latency $latency)")
          res += c.io.out.bits.peek().litValue
        }
        c.clock.step(1)
      }
    }
    res.toSeq
  }

  /** Compares every configuration against the combinational module. */
  def compare(name: String, bitW: PolynomialBitWidth, clamp: Boolean,
    stages: Seq[PipelineStageConfig]): Unit = {

    val r = new Random(19660809)
    val vectors = randVectors(bitW, nVec, r)
    val ref = run(bitW, PipelineStageConfig.none, clamp, vectors)
    assert(ref.length == nVec)
    for (stage <- stages) {
      val got = run(bitW, stage, clamp, vectors)
      assert(got.length == nVec)
      for (i <- 0 until nVec) {
        assert(got(i) == ref(i),
          f"$name [${stage.getString}] vector $i: ${got(i)} != ${ref(i)} " +
          f"(coeffs ${vectors(i)._1.mkString(",")}, x ${vectors(i)._2})")
      }
    }
  }

  describe("PolynomialEval.stagePoints") {

    it("declares one boundary per multiply-add, then the output") {
      assert(PolynomialEval.stagePoints(0).map(_.name) == Seq("out"))
      assert(PolynomialEval.stagePoints(1).map(_.name) == Seq("step1", "out"))
      assert(PolynomialEval.stagePoints(3).map(_.name) ==
        Seq("step1", "step2", "step3", "out"))
    }

    it("balances the multiply-adds when the module places the registers") {
      val s = PipelineStageConfig.default(2).resolve(PolynomialEval.stagePoints(2))
      assert(s.at("step1") == 1)
      assert(s.at("step2") == 1)
      assert(s.at("out") == 0)
      assert(s.total == 2)
    }
  }

  describe("PolynomialEval: a plain Int") {

    it("puts every register at the output") {
      var seen: PipelineStages = null
      ChiselStage.emitCHIRRTL {
        val m = new PolynomialEval(bitW2, 2)
        seen = m.stages
        m
      }
      assert(seen.at("out") == 2)
      assert(seen.at("in") == 0)
      assert(seen.at("step1") == 0)
      assert(seen.at("step2") == 0)
      assert(seen.total == 2)
    }

    it("puts every register at the output when clamping too") {
      var seen: PipelineStages = null
      ChiselStage.emitCHIRRTL {
        val m = new PolynomialEval(bitW2, 3, clampNegToZero = true)
        seen = m.stages
        m
      }
      assert(seen.at("out") == 3)
      assert(seen.total == 3)
    }
  }

  describe("PolynomialEval: the pipeline registers do not change the result") {

    it(f"order ${bitW2.nOrder} matches the combinational module") {
      compare("order2", bitW2, false, Seq(
        PipelineStageConfig.default(2),
        PipelineStageConfig.specify("step1" -> 1, "step2" -> 1, "out" -> 1),
        PipelineStageConfig.atOut(2)))
    }

    it(f"order ${bitW3.nOrder} matches the combinational module, with clamping") {
      compare("order3", bitW3, true, Seq(
        PipelineStageConfig.default(2),
        PipelineStageConfig.specify(
          "step1" -> 1, "step2" -> 1, "step3" -> 1, "out" -> 1)))
    }
  }
}
