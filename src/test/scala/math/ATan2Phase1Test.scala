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
import rial.util.ScalaUtil._

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue

class ATan2Phase1Test extends AnyFunSpec
    with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of random test cases per case (default 1000)"))
  val r = new Random(123456789)

  def generateRealWithin(from: Double, to: Double, spec: RealSpec, r: Random) = {
    val rD: Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  private def runtest(
    spec:  RealSpec,
    cfg:   MathFuncConfig,
    stage: MathFuncPipelineConfig,
    generatorStr: String,
    generatorX:      (RealSpec, Random) => RealGeneric,
    generatorYoverX: (RealSpec, Random) => RealGeneric,
    reference: (RealGeneric, RealGeneric) => RealGeneric,
  ) = {
    val total      = stage.total
    val pipeconfig = stage.getString

    it(f"atan2Phase1(y, x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new MathFunctions(cfg, spec, stage)) { c =>
        val q = new Queue[(BigInt, BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val yOverX = generatorYoverX(spec, r)
          val x      = generatorX(spec, r)
          val x0     = x.toDouble
          val y      = new RealGeneric(spec, x0 * yOverX.toDouble)

          val z0r = reference(y, x)
          q += ((x.value.toBigInt, y.value.toBigInt, z0r.value.toBigInt))

          c.io.sel.poke(cfg.signal(ATan2Phase1))
          c.io.x.poke(x.value.toBigInt.U(spec.W.W))
          c.io.y.get.poke(y.value.toBigInt.U(spec.W.W))

          val zi = c.io.z.peek().litValue.toBigInt
          if (i > total) {
            val (xid, yid, z0d) = q.dequeue()

            val zisgn = bit(spec.W - 1, zi).toInt
            val ziexp = slice(spec.manW, spec.exW, zi)
            val ziman = zi & maskSL(spec.manW)
            val z0dsgn = bit(spec.W - 1, z0d).toInt
            val z0dexp = slice(spec.manW, spec.exW, z0d)
            val z0dman = z0d & maskSL(spec.manW)

            assert(zi == z0d, f"x = ${xid.toString(16)}, y = ${yid.toString(16)}, " +
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
  val phase1RefFP32 = (y: RealGeneric, x: RealGeneric) =>
    ATan2Phase1SimGeneric(fp32Config.coeffs(Reciprocal), y, x)

  val noPipeStage     = new MathFuncPipelineConfig(0, 0, 0, 0, 0, 0)
  val simplePipeStage = new MathFuncPipelineConfig(1, 1, 1, 1, 0, 0)
  val mulPipeStage    = new MathFuncPipelineConfig(2, 1, 1, 2, 0, 1)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test y/x >  2^24",         generateRealWithin(-1.0, 1.0, _, _),         generateRealWithin( pow(2.0,  24),  pow(2.0, 128), _, _), phase1RefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test 1 < y/x < 2^12",      generateRealWithin(-1.0, 1.0, _, _),         generateRealWithin( 1.0,             pow(2.0,   8), _, _), phase1RefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test 2^-12 < y/x < 1",     generateRealWithin(-1.0, 1.0, _, _),         generateRealWithin( pow(2.0, -12),  1.0,           _, _), phase1RefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test y/x < 2^-12",         generateRealWithin(-1.0, 1.0, _, _),         generateRealWithin( 0.0,             pow(2.0, -12), _, _), phase1RefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test y/x ~ 1",             generateRealWithin(-1.0, 1.0, _, _),         generateRealWithin( 0.999998,        1.000002,      _, _), phase1RefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test 2^-12 < y/x < 1",     generateRealWithin(-1.0, 1.0, _, _),         generateRealWithin( pow(2.0, -12),  1.0,           _, _), phase1RefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test y/x ~ 1",             generateRealWithin(-1.0, 1.0, _, _),         generateRealWithin( 0.999998,        1.000002,      _, _), phase1RefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, mulPipeStage,
    "Test 2^-12 < y/x < 1",     generateRealWithin(-1.0, 1.0, _, _),         generateRealWithin( pow(2.0, -12),  1.0,           _, _), phase1RefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, mulPipeStage,
    "Test y/x ~ 1",             generateRealWithin(-1.0, 1.0, _, _),         generateRealWithin( 0.999998,        1.000002,      _, _), phase1RefFP32)

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = MathFuncConfig.standardBF16
  val phase1RefBF16 = (y: RealGeneric, x: RealGeneric) =>
    ATan2Phase1SimGeneric(bf16Config.coeffs(Reciprocal), y, x)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test y/x > 2^24",      generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, 24),  pow(2.0, 128), _, _), phase1RefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test 1 < y/x < 2^12",  generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(1.0,           pow(2.0,   8), _, _), phase1RefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, -12), 1.0,           _, _), phase1RefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test y/x < 2^-12",     generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(0.0,           pow(2.0, -12), _, _), phase1RefBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, -12), 1.0,           _, _), phase1RefBF16)
}
