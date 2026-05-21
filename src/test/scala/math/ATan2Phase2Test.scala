package rial.tests

import org.scalatest._

import chisel3._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.implicits._

import rial.math._
import rial.funcs.FuncKind._
import rial.arith._
import rial.util.ScalaUtil._

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue

class ATan2Phase2Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test atan2Phase2 with funcs/MathFunctions"

  val n = 1000
  val r = new Random(123456789)

  def generateRealWithin(from: Double, to: Double, spec: RealSpec, r: Random) = {
    val rD: Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  // Full chain reference: takes original (y, x), applies phase 1 (SW), then phase 2.
  // The HW under test only runs phase 2 and consumes phase-1's output as its `x`.
  private def fullRef(cfg: MathFuncConfig): (RealGeneric, RealGeneric) => (RealGeneric, RealGeneric) = {
    val tRec  = cfg.coeffs(Reciprocal)
    val tATan = cfg.coeffs(ATan2Phase2)
    (y: RealGeneric, x: RealGeneric) => {
      val phase1 = ATan2Phase1SimGeneric(tRec,  y, x)
      val phase2 = ATan2Phase2SimGeneric(tATan, phase1)
      (phase1, phase2)
    }
  }

  private def runtest(
    spec:  RealSpec,
    cfg:   MathFuncConfig,
    stage: MathFuncPipelineConfig,
    n: Int, generatorStr: String,
    generatorX: (RealSpec, Random) => RealGeneric,
    generatorY: (RealSpec, Random) => RealGeneric,
    ref: (RealGeneric, RealGeneric) => (RealGeneric, RealGeneric),
  ) = {
    val total      = stage.total
    val pipeconfig = stage.getString

    it should f"atan2Phase2(y, x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}" in {
      test(new MathFunctions(cfg, spec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val x = generatorX(spec, r)
          val y = generatorY(spec, r)
          val (p1, p2) = ref(y, x)
          q += ((p1.value.toBigInt, p2.value.toBigInt))

          c.io.sel.poke(cfg.signal(ATan2Phase2))
          c.io.x.poke(p1.value.toBigInt.U(spec.W.W))
          // y is unused by Phase2 but the port must still be driven.
          if (c.io.y.isDefined) c.io.y.get.poke(0.U(spec.W.W))

          val zi = c.io.z.peek().litValue.toBigInt
          if (i > total) {
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

  val noPipeStage     = new MathFuncPipelineConfig(0, 0, 0, 0, 0, 0)
  val simplePipeStage = new MathFuncPipelineConfig(1, 1, 1, 1, 0, 0)
  val mulPipeStage    = new MathFuncPipelineConfig(1, 1, 1, 2, 0, 1)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, "Test 2^12 < |y/x|",
    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, 12), pow(2.0, 24), _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, "Test 1 < |y/x| < 2^12",
    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(1.0,          pow(2.0,  8), _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, "Test 2^-12 < |y/x| < 1",
    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, -12), 1.0,         _, _), refFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, "Test |y/x| < 2^-12",
    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(0.0,          pow(2.0, -12), _, _), refFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    n, "Test 2^-12 < |y/x| < 1",
    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, -12), 1.0,         _, _), refFP32)

  runtest(RealSpec.Float32Spec, fp32Config, mulPipeStage,
    n, "Test 2^-12 < |y/x| < 1",
    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, -12), 1.0,         _, _), refFP32)

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = MathFuncConfig.standardBF16
  val refBF16 = fullRef(bf16Config)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    n, "Test 1 < |y/x| < 2^12",
    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(1.0,          pow(2.0,  8), _, _), refBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    n, "Test 2^-12 < |y/x| < 1",
    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, -12), 1.0,         _, _), refBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    n, "Test 2^-12 < |y/x| < 1",
    generateRealWithin(-1.0, 1.0, _, _), generateRealWithin(pow(2.0, -12), 1.0,         _, _), refBF16)
}
