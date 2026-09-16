package rial.tests

import org.scalatest._

import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import org.scalatest.funspec.AnyFunSpec

import spire.math.SafeLong
import spire.implicits._

import rial.math._
import rial.math.FuncKind._
import rial.arith._
import rial.util.PipelineStageConfig
import rial.util.ScalaUtil._

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue

class ACosPhase2Test extends AnyFunSpec
    with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of random test cases per case (default 1000)"))
  val r = new Random(123456789)

  def generateRealWithin(from: Double, to: Double, spec: RealSpec, r: Random) = {
    val rD: Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  var counter = 0
  val specialValues = Seq(
       0.0, -0.0,  1.0, -1.0,
       1.1, -1.1,
       Double.PositiveInfinity, Double.NegativeInfinity,
    )
  def generateSpecialValues(spec: RealSpec, r: Random) = {
    val idx = counter
    counter += 1
    if (counter >= specialValues.length) counter = 0
    new RealGeneric(spec, specialValues(idx))
  }

  // Full chain reference: takes original x, applies phase 1 (SW), then phase 2.
  // The HW under test receives phase-1's output as its input.
  private def fullRef(cfg: MathFuncConfig): RealGeneric => (RealGeneric, RealGeneric) = {
    val tSqrt = cfg.coeffs(Sqrt)
    val tACos = cfg.coeffs(ACosPhase2)
    (x: RealGeneric) => {
      val phase1 = ACosPhase1SimGeneric(tSqrt, x)
      val phase2 = ACosPhase2SimGeneric(tACos, phase1)
      (phase1, phase2)
    }
  }

  private def runtest(
    spec:  RealSpec,
    cfg:   MathFuncConfig,
    stage: PipelineStageConfig,
    generatorStr: String, generator: (RealSpec, Random) => RealGeneric,
    ref: RealGeneric => (RealGeneric, RealGeneric),
  ) = {
    val total      = stage.resolve(MathFunctions.stagePoints(cfg)).total
    val pipeconfig = MathTestUtil.label(cfg, stage)

    it(f"acosPhase2(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new MathFunctions(cfg, stage)) { c =>
        counter = 0
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          // Generate original x, compute phase-1 in SW, feed phase-1's output to HW.
          val x      = generator(spec, r)
          val (p1, p2) = ref(x)
          q += ((p1.value.toBigInt, p2.value.toBigInt))

          c.io.in.valid.poke(true.B)
          c.io.in.bits.sel.poke(cfg.signal(ACosPhase2))
          c.io.in.bits.x.poke(p1.value.toBigInt.U(spec.W.W))

          val zi = c.io.out.bits.z.peek().litValue.toBigInt
          if (i > total) {
            assert(c.io.out.valid.peek().litToBoolean,
              f"io.out.valid is not set at cycle $i%d (latency $total%d)")
            val (xid, z0d) = q.dequeue()

            val xidsgn = bit(spec.W - 1, xid).toInt
            val xidexp = slice(spec.manW, spec.exW, xid)
            val xidman = xid & maskSL(spec.manW)

            val zisgn = bit(spec.W - 1, zi).toInt
            val ziexp = slice(spec.manW, spec.exW, zi)
            val ziman = zi & maskSL(spec.manW)

            val z0dsgn = bit(spec.W - 1, z0d).toInt
            val z0dexp = slice(spec.manW, spec.exW, z0d)
            val z0dman = z0d & maskSL(spec.manW)

            assert(zi == z0d, f"phase1 in = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
                              f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != " +
                              f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
          }
          c.clock.step(1)
        }
      }
    }
  }

  // --------------------------------------------------------------------------
  // FP32

  val fp32Config = MathFuncConfig.standardFP32
  val refFP32 = fullRef(fp32Config)

  val noPipeStage     = PipelineStageConfig.none
  val simplePipeStage = PipelineStageConfig.default(4)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test 0 < x <= 2^-4",   generateRealWithin(0.0, pow(2.0, -4), _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test 2^-4 < x < 1-2^-4", generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4), _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test 1-2^-4 < x <= 1", generateRealWithin(1.0 - pow(2.0, -4), 1.0, _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test -1 <= x < 0",     generateRealWithin(-1.0, 0.0, _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test special values",  generateSpecialValues(_, _), refFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test 2^-4 < x < 1-2^-4", generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4), _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test -1 <= x < 0",     generateRealWithin(-1.0, 0.0, _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test special values",  generateSpecialValues(_, _), refFP32)

  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test 2^-4 < x < 1-2^-4", generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4), _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test -1 <= x < 0",     generateRealWithin(-1.0, 0.0, _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test special values",  generateSpecialValues(_, _), refFP32)

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = MathFuncConfig.standardBF16
  val refBF16 = fullRef(bf16Config)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test 2^-4 < x < 1-2^-4", generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4), _, _), refBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test -1 <= x < 0",     generateRealWithin(-1.0, 0.0, _, _), refBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test special values",  generateSpecialValues(_, _), refBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test 2^-4 < x < 1-2^-4", generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4), _, _), refBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test special values",  generateSpecialValues(_, _), refBF16)
}
