package rial.tests

import org.scalatest.funspec.AnyFunSpec

import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import circt.stage.ChiselStage

import rial.arith.MulImpl
import rial.math.MathFuncConfig
import rial.table.{PolynomialBitWidth, PolynomialCoefficients, PolynomialEval}
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
    vectors: Seq[(Seq[BigInt], BigInt)],
    mulImpl: MulImpl = MulImpl.plain): Seq[BigInt] = {

    val res = ArrayBuffer.empty[BigInt]
    simulate(new PolynomialEval(bitW, stage, clamp, mulImpl)) { c =>
      val latency = c.getStage
      assert(latency ==
        stage.resolve(PolynomialEval.stagePoints(bitW.nOrder, mulImpl)).total)
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
    stages: Seq[PipelineStageConfig], mulImpl: MulImpl = MulImpl.plain): Unit = {

    val r = new Random(19660809)
    val vectors = randVectors(bitW, nVec, r)
    val ref = run(bitW, PipelineStageConfig.none, clamp, vectors)
    assert(ref.length == nVec)
    for (stage <- stages) {
      val got = run(bitW, stage, clamp, vectors, mulImpl)
      assert(got.length == nVec)
      for (i <- 0 until nVec) {
        assert(got(i) == ref(i),
          f"$name [${stage.getString}] vector $i: ${got(i)} != ${ref(i)} " +
          f"(coeffs ${vectors(i)._1.mkString(",")}, x ${vectors(i)._2})")
      }
    }
  }

  /** The result the module must produce, from the software model of the
   * polynomial in [[PolynomialEval.eval]] and the clamp of its output stage.
   */
  def swEval(bitW: PolynomialBitWidth, clamp: Boolean,
    cs: Seq[BigInt], x: BigInt): BigInt = {

    val xu = if (bitW.nOrder == 0) {
      BigInt(0)
    } else {
      x + (BigInt(1) << bitW.inPrec)
    }
    val coeffs = PolynomialCoefficients(bitW, Seq(cs.map(_.toLong)))
    PolynomialEval.clamp(PolynomialEval.eval(0, xu, coeffs), bitW, clamp)
  }

  describe("PolynomialEval.stagePoints") {

    it("declares the multiplier and the end of every multiply-add, then the output") {
      assert(PolynomialEval.stagePoints(0).map(_.name) == Seq("out"))
      assert(PolynomialEval.stagePoints(1).map(_.name) ==
        Seq("prod1", "step1", "out"))
      assert(PolynomialEval.stagePoints(3).map(_.name) ==
        Seq("prod1", "step1", "prod2", "step2", "prod3", "step3", "out"))
    }

    it("cuts the multipliers first when the module places the registers") {
      val s = PipelineStageConfig.default(2).resolve(PolynomialEval.stagePoints(2))
      assert(s.at("prod1") == 1)
      assert(s.at("prod2") == 1)
      assert(s.at("step1") == 0)
      assert(s.at("step2") == 0)
      assert(s.at("out") == 0)
      assert(s.total == 2)
    }
  }

  describe("PolynomialEval: a plain Int") {

    it("puts every register at the output") {
      var seen: PipelineStages = null
      ChiselStage.emitCHIRRTL {
        val m = new PolynomialEval(bitW2, PipelineStageConfig.atOut(2))
        seen = m.stages
        m
      }
      assert(seen.at("out") == 2)
      assert(seen.at("in") == 0)
      assert(seen.at("prod1") == 0)
      assert(seen.at("step1") == 0)
      assert(seen.at("prod2") == 0)
      assert(seen.at("step2") == 0)
      assert(seen.total == 2)
    }

    it("puts every register at the output when clamping too") {
      var seen: PipelineStages = null
      ChiselStage.emitCHIRRTL {
        val m = new PolynomialEval(bitW2, PipelineStageConfig.atOut(3), clampNegToZero = true)
        seen = m.stages
        m
      }
      assert(seen.at("out") == 3)
      assert(seen.total == 3)
    }
  }

  describe("PolynomialEval: the module matches the software model") {

    for ((nm, bw, cl) <- Seq(("order2", bitW2, false), ("order3", bitW3, true));
         (mn, mi) <- Seq(("", MulImpl.plain), (", split at 17", MulImpl.split(17)),
           (", a multiply-add at 17", MulImpl.macc(17)),
           (", split small at 17", MulImpl.splitSmall(17)),
           (", a small multiply-add at 17", MulImpl.maccSmall(17)))) {
      it(f"$nm reproduces PolynomialEval.eval bit for bit$mn") {
        val r = new Random(19660809)
        val vectors = randVectors(bw, nVec, r)
        val got = run(bw, PipelineStageConfig.none, cl, vectors, mi)
        assert(got.length == nVec)
        for (i <- 0 until nVec) {
          val (cs, x) = vectors(i)
          val want = swEval(bw, cl, cs, x)
          assert(got(i) == want,
            f"$nm$mn vector $i: ${got(i)} != $want " +
            f"(coeffs ${cs.mkString(",")}, x $x)")
        }
      }
    }
  }

  describe("PolynomialEval: a result outside [0, 1)") {

    /** The same shape with room above and below `[0.0, 1.0)` in the last
     * accumulator, so that a coefficient alone can leave the range. */
    val bitWWide: PolynomialBitWidth = PolynomialBitWidth(
      nOrder   = 2,
      addrW    = bitW2.addrW,
      inBitW   = bitW2.inBitW,
      outBitW  = 16,
      coefBits = Seq(20, 14, 10),
      calcBits = Seq(20, 14, 10))

    val resMax = (BigInt(1) << bitWWide.outBitW) - 1
    val xMax   = (BigInt(1) << bitWWide.inPrec) - 1
    val xMin   = -(BigInt(1) << bitWWide.inPrec)
    val c1Max  = (BigInt(1) << 13) - 1

    /** coefficients and position whose result lands below zero */
    val below = Seq(
      (Seq(BigInt(-1), BigInt(0), BigInt(0)), BigInt(0)),
      (Seq(-(BigInt(1) << 17), BigInt(0), BigInt(0)), BigInt(0)),
      (Seq(BigInt(-8000), c1Max, BigInt(0)), xMin)) // the last step takes it down

    /** ... and above the range */
    val above = Seq(
      (Seq(resMax + 1, BigInt(0), BigInt(0)), BigInt(0)),
      (Seq((BigInt(1) << 19) - 1, BigInt(0), BigInt(0)), BigInt(0)),
      (Seq(resMax - 100, c1Max, BigInt(0)), xMax)) // the last step takes it up

    /** ... and inside it, where nothing is clamped */
    val inside = Seq(
      (Seq(BigInt(0), BigInt(0), BigInt(0)), BigInt(0)),
      (Seq(resMax - 1, BigInt(0), BigInt(0)), BigInt(0)),
      (Seq(resMax, BigInt(0), BigInt(0)), BigInt(0)))

    for (cl <- Seq(false, true)) {
      it(f"saturates where PolynomialEval.clamp does (clampNegToZero = $cl)") {
        val vectors = below ++ above ++ inside
        val got = run(bitWWide, PipelineStageConfig.none, cl, vectors)
        for (i <- vectors.indices) {
          val (cs, x) = vectors(i)
          assert(got(i) == swEval(bitWWide, cl, cs, x),
            f"vector $i: ${got(i)} != ${swEval(bitWWide, cl, cs, x)} " +
            f"(coeffs ${cs.mkString(",")}, x $x)")
        }
        // a negative result reaches the output as its two's complement bits,
        // which are above the range unless the module clamps them to zero
        val want = below.map(_ => if (cl) { BigInt(0) } else { resMax }) ++
                   above.map(_ => resMax) ++
                   inside.map { case (cs, _) => cs(0) }
        assert(got == want, f"$got != $want")
      }
    }
  }

  describe("PolynomialEval: the pipeline registers do not change the result") {

    it(f"order ${bitW2.nOrder} matches the combinational module") {
      compare("order2", bitW2, false, Seq(
        PipelineStageConfig.default(2),
        PipelineStageConfig.specify(
          "prod1" -> 1, "step1" -> 1, "prod2" -> 1, "step2" -> 1, "out" -> 1),
        PipelineStageConfig.atOut(2)))
    }

    it(f"order ${bitW3.nOrder} matches the combinational module, with clamping") {
      compare("order3", bitW3, true, Seq(
        PipelineStageConfig.default(2),
        PipelineStageConfig.specify(
          "prod1" -> 1, "step1" -> 1, "prod2" -> 1, "step2" -> 1,
          "prod3" -> 1, "step3" -> 1, "out" -> 1)))
    }

    it(f"order ${bitW3.nOrder} matches it with a split product too") {
      compare("order3 split", bitW3, true, Seq(
        PipelineStageConfig.default(2),
        PipelineStageConfig.specify(
          "prod1" -> 1, "step1" -> 1, "prod2" -> 1, "step2" -> 1,
          "prod3" -> 1, "step3" -> 1, "out" -> 1)), MulImpl.split(17))
    }

    it(f"order ${bitW3.nOrder} matches it with a multiply-add too") {
      compare("order3 macc", bitW3, true, Seq(
        PipelineStageConfig.default(2),
        PipelineStageConfig.specify(
          "prod1.pp" -> 1, "step1" -> 1, "prod2.pp" -> 1, "step2" -> 1,
          "prod3.pp" -> 1, "step3" -> 1, "out" -> 1),
        PipelineStageConfig.specify(
          "prod1" -> 1, "step1" -> 1, "prod2" -> 1, "step2" -> 1,
          "prod3" -> 1, "step3" -> 1, "out" -> 1)), MulImpl.macc(17))
    }
  }
}
