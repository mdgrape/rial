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

class Log2FixedTest extends AnyFunSpec
    with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of random test cases per case (default 1000)"))
  val r = new Random(123456789)

  /** any bit pattern of the format, specials and negatives included. */
  def generateFullRange(spec: RealSpec, r: Random) = {
    new RealGeneric(spec, SafeLong(BigInt(spec.W, r)))
  }

  def generateRealWithin(from: Double, to: Double, spec: RealSpec, r: Random) = {
    new RealGeneric(spec, from + r.nextDouble() * (to - from))
  }

  var nearOneCounter = 0
  /** x = 1 +- 2^-k for k = 1 .. 23, and random arguments around them. */
  def generateNearOne(spec: RealSpec, r: Random) = {
    val k = 1 + (nearOneCounter % 23)
    nearOneCounter += 1
    val d = if ((nearOneCounter % 2) == 0) {
      pow(2.0, -k)
    } else {
      pow(2.0, -k) * (1.0 + r.nextDouble())
    }
    new RealGeneric(spec, if (r.nextBoolean()) { 1.0 + d } else { 1.0 - d })
  }

  var counter = 0
  val specialValues = Seq(
      +0.0,
      -0.0,
      1.0,
      2.0,
      0.5,
      4.0,
      0.25,
      Double.PositiveInfinity,
      Double.NegativeInfinity,
      Double.NaN,
      -1.0,
      -128.0,
    )
  def generateSpecialValues(spec: RealSpec, r: Random) = {
    val idx = counter
    counter += 1
    if (counter >= specialValues.length) counter = 0
    new RealGeneric(spec, specialValues(idx))
  }

  var smallestCounter = 0
  /** x in the smallest normal bands of the format. */
  def generateSmallestNormals(spec: RealSpec, r: Random) = {
    val ex = 1 + (smallestCounter % 4)
    smallestCounter += 1
    new RealGeneric(spec, (SafeLong(ex) << spec.manW) +
                          SafeLong(BigInt(spec.manW, r)))
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

    it(f"log2Fixed(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new MathFunctions(cfg, stage)) { c =>
        counter = 0
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val xi  = generator(spec, r)
          val z0r = reference(xi)
          q += ((xi.value.toBigInt, z0r.value.toBigInt))

          c.io.in.valid.poke(true.B)
          c.io.in.bits.sel.poke(cfg.signal(Log2Fixed))
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
  val log2FixedRefFP32 =
    Log2FixedSimGeneric(fp32Config.coeffs(Log2Fixed), _: RealGeneric)

  val noPipeStage     = PipelineStageConfig.none
  val simplePipeStage = PipelineStageConfig.default(4)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test full range", generateFullRange(_, _), log2FixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Within (0, 128)", generateRealWithin(0.0, 128.0, _, _), log2FixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test around 1", generateNearOne(_, _), log2FixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), log2FixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test the smallest normals", generateSmallestNormals(_, _), log2FixedRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test full range", generateFullRange(_, _), log2FixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test around 1", generateNearOne(_, _), log2FixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test special values", generateSpecialValues(_, _), log2FixedRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test full range", generateFullRange(_, _), log2FixedRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test around 1", generateNearOne(_, _), log2FixedRefFP32)

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = MathFuncConfig.standardFixedBF16
  val log2FixedRefBF16 =
    Log2FixedSimGeneric(bf16Config.coeffs(Log2Fixed), _: RealGeneric)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test full range", generateFullRange(_, _), log2FixedRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Within (0, 128)", generateRealWithin(0.0, 128.0, _, _), log2FixedRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), log2FixedRefBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test full range", generateFullRange(_, _), log2FixedRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test around 1", generateNearOne(_, _), log2FixedRefBF16)
}
