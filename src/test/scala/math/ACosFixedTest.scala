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

class ACosFixedTest extends AnyFunSpec
    with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of random test cases per case (default 1000)"))
  val r = new Random(123456789)

  /** a random x in [-1, 1). */
  def generateWithin(spec: RealSpec, r: Random) = {
    new RealGeneric(spec, 2.0 * (r.nextDouble() - 0.5))
  }

  var edgeCounter = 0
  /** x next to +1 or -1: +-(1 - 2^-k) with the bits below it random. */
  def generateNearEdge(spec: RealSpec, r: Random) = {
    val k   = 1 + (edgeCounter % 24)
    val sgn = if ((edgeCounter / 24) % 2 == 0) { 1.0 } else { -1.0 }
    edgeCounter += 1
    val x = new RealGeneric(spec, sgn * (1.0 - pow(2.0, -k)))
    new RealGeneric(spec, x.value + SafeLong(r.nextInt(4)))
  }

  /** a tiny x, where acos(x) is pi/2 to the working precision. */
  def generateTiny(spec: RealSpec, r: Random) = {
    val ex = 1 + r.nextInt(spec.exBias - 20)
    new RealGeneric(spec, (SafeLong(r.nextInt(2)) << (spec.W - 1)) +
                          (SafeLong(ex) << spec.manW) +
                          SafeLong(BigInt(spec.manW, r)))
  }

  var counter = 0
  val specialValues = Seq(
      +0.0,       -0.0,
      1.0,        -1.0,
      0.5,        -0.5,
      1.5,        -1.5,
      Double.NaN, Double.PositiveInfinity, Double.NegativeInfinity,
    )
  def generateSpecialValues(spec: RealSpec, r: Random) = {
    val idx = counter
    counter += 1
    if (counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }

  private def runtest(
    spec:  RealSpec,
    cfg:   MathFuncConfig,
    stage: PipelineStageConfig,
    generatorStr: String, generator: (RealSpec, Random) => RealGeneric,
    reference: RealGeneric => RealGeneric,
  ) = {
    val total      = stage.resolve(MathFunctions.stagePoints(cfg)).total
    val pipeconfig = MathTestUtil.label(cfg, stage)

    it(f"acosFixed(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new MathFunctions(cfg, spec, stage)) { c =>
        counter     = 0
        edgeCounter = 0
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val xi  = generator(spec, r)
          val z0r = reference(xi)
          q += ((xi.value.toBigInt, z0r.value.toBigInt))

          c.io.in.valid.poke(true.B)
          c.io.in.bits.sel.poke(cfg.signal(ACosFixed))
          c.io.in.bits.x.poke(xi.value.toBigInt.U(spec.W.W))

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

            assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
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

  val fp32Config = MathFuncConfig.standardFixedFP32
  // Use the same merged coefficients the hardware sees, so HW == SW bit-exact.
  val acosFixedRefFP32 =
    ACosFixedSimGeneric(fp32Config.coeffs(ACosFixed), _: RealGeneric)

  val noPipeStage     = PipelineStageConfig.none
  val simplePipeStage = PipelineStageConfig.default(4)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Within [-1, 1)", generateWithin(_, _), acosFixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test x next to +-1", generateNearEdge(_, _), acosFixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test a tiny x", generateTiny(_, _), acosFixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), acosFixedRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test Within [-1, 1)", generateWithin(_, _), acosFixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test x next to +-1", generateNearEdge(_, _), acosFixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test special values", generateSpecialValues(_, _), acosFixedRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test Within [-1, 1)", generateWithin(_, _), acosFixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test x next to +-1", generateNearEdge(_, _), acosFixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test special values", generateSpecialValues(_, _), acosFixedRefFP32)

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = MathFuncConfig.standardFixedBF16
  val acosFixedRefBF16 =
    ACosFixedSimGeneric(bf16Config.coeffs(ACosFixed), _: RealGeneric)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Within [-1, 1)", generateWithin(_, _), acosFixedRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test x next to +-1", generateNearEdge(_, _), acosFixedRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test a tiny x", generateTiny(_, _), acosFixedRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), acosFixedRefBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test Within [-1, 1)", generateWithin(_, _), acosFixedRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test special values", generateSpecialValues(_, _), acosFixedRefBF16)
}
