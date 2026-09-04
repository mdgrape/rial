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

class LogTest extends AnyFunSpec
    with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of random test cases per case (default 1000)"))
  val r = new Random(123456789)

  def generateRealWithin(from: Double, to: Double, spec: RealSpec, r: Random) = {
    val rD: Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  var counter = 0
  val specialValues = Seq(
      1.0,
      math.E / 2.0,
      math.E,
      math.E * 2.0,
      pow(math.E, 2),
      pow(math.E, 3),
      pow(math.E, 4),
    )
  def generateSpecialValues(spec: RealSpec, r: Random) = {
    val idx = counter
    counter += 1
    if (counter >= specialValues.length) counter = 0
    new RealGeneric(spec, specialValues(idx))
  }

  private def logRef(cfg: MathFuncConfig): RealGeneric => RealGeneric = {
    val lt = cfg.tables(Log).asInstanceOf[LogTable]
    val cN = lt.genNormalCoeffs  (cfg.bitW)
    val cP = lt.genSmallPosCoeffs(cfg.bitW)
    val cG = lt.genSmallNegCoeffs(cfg.bitW)
    (x: RealGeneric) => LogSimGeneric(cN, cP, cG, x)
  }

  private def runtest(
    spec:  RealSpec,
    cfg:   MathFuncConfig,
    stage: PipelineStageConfig,
    generatorStr: String, generator: (RealSpec, Random) => RealGeneric,
    reference: RealGeneric => RealGeneric,
  ) = {
    val total      = stage.resolve(MathFunctions.stagePoints(cfg)).total
    val pipeconfig = stage.getString

    it(f"log(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new MathFunctions(cfg, spec, stage)) { c =>
        counter = 0
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val xi  = generator(spec, r)
          val z0r = reference(xi)
          q += ((xi.value.toBigInt, z0r.value.toBigInt))

          c.io.in.valid.poke(true.B)
          c.io.in.bits.sel.poke(cfg.signal(Log))
          c.io.in.bits.x.poke(xi.value.toBigInt.U(spec.W.W))

          val zi = c.io.out.bits.peek().litValue.toBigInt
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

  val fp32Config = MathFuncConfig.standardFP32
  val logRefFP32 = logRef(fp32Config)

  val noPipeStage     = PipelineStageConfig.none
  val simplePipeStage = PipelineStageConfig.default(4)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test [2, 128]", generateRealWithin(2.0, 128.0, _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test [1, 2)",   generateRealWithin(1.0, 2.0 - pow(2.0, -23), _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test [0.5, 1)", generateRealWithin(0.5, 1.0 - pow(2.0, -24), _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test [0, 0.5)", generateRealWithin(0.0, 0.5 - pow(2.0, -25), _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test negative", generateRealWithin(-128.0, 0.0, _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), logRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test [2, 128]", generateRealWithin(2.0, 128.0, _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test [1, 2)",   generateRealWithin(1.0, 2.0 - pow(2.0, -23), _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test [0.5, 1)", generateRealWithin(0.5, 1.0 - pow(2.0, -24), _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test [0, 0.5)", generateRealWithin(0.0, 0.5 - pow(2.0, -25), _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test special values", generateSpecialValues(_, _), logRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test [2, 128]", generateRealWithin(2.0, 128.0, _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test [1, 2)",   generateRealWithin(1.0, 2.0 - pow(2.0, -23), _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test [0.5, 1)", generateRealWithin(0.5, 1.0 - pow(2.0, -24), _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test [0, 0.5)", generateRealWithin(0.0, 0.5 - pow(2.0, -25), _, _), logRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test special values", generateSpecialValues(_, _), logRefFP32)

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = MathFuncConfig.standardBF16
  val logRefBF16 = logRef(bf16Config)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test [2, 128]", generateRealWithin(2.0, 128.0, _, _), logRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test [1, 2)",   generateRealWithin(1.0, 2.0 - pow(2.0, -7), _, _), logRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test [0.5, 1)", generateRealWithin(0.5, 1.0 - pow(2.0, -8), _, _), logRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test [0, 0.5)", generateRealWithin(0.0, 0.5 - pow(2.0, -9), _, _), logRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), logRefBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test [2, 128]", generateRealWithin(2.0, 128.0, _, _), logRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test [1, 2)",   generateRealWithin(1.0, 2.0 - pow(2.0, -7), _, _), logRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test special values", generateSpecialValues(_, _), logRefBF16)
}
