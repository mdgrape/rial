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
import rial.math.FuncKind._
import rial.arith._
import rial.util.ScalaUtil._

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue

class CosTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test cos with funcs/MathFunctions"

  val n = 1000
  val r = new Random(123456789)

  def generateRealWithin(p: Double, spec: RealSpec, r: Random) = {
    val rD: Double = (r.nextDouble() - 0.5) * p
    val x = new RealGeneric(spec, rD)
    new RealGeneric(spec, (x.value & (maskSL(spec.exW + 1) << spec.manW)) +
                          SafeLong(BigInt(spec.manW, r)))
  }

  var counter = 0
  val specialValues = Seq(
      +0.0,                -0.0,
      math.Pi / 4.0,       -math.Pi / 4.0,
      math.Pi / 3.0,       -math.Pi / 3.0,
      math.Pi / 2.0,       -math.Pi / 2.0,
      math.Pi * 2.0 / 3.0, -math.Pi * 2.0 / 3.0,
      math.Pi * 3.0 / 4.0, -math.Pi * 3.0 / 4.0,
      math.Pi,             -math.Pi,
      math.Pi * 5.0 / 4.0, -math.Pi * 5.0 / 4.0,
      math.Pi * 4.0 / 3.0, -math.Pi * 4.0 / 3.0,
      math.Pi * 3.0 / 2.0, -math.Pi * 3.0 / 2.0,
      math.Pi * 5.0 / 3.0, -math.Pi * 5.0 / 3.0,
      math.Pi * 7.0 / 4.0, -math.Pi * 7.0 / 4.0,
      math.Pi * 2.0,       -math.Pi * 2.0
    )
  def generateSpecialValues(spec: RealSpec, r: Random) = {
    val idx = counter
    counter += 1
    if (counter >= specialValues.length) counter = 0
    new RealGeneric(spec, specialValues(idx))
  }

  private def runtest(
    spec:  RealSpec,
    cfg:   MathFuncConfig,
    stage: MathFuncPipelineConfig,
    n: Int, generatorStr: String, generator: (RealSpec, Random) => RealGeneric,
    reference: RealGeneric => RealGeneric,
  ) = {
    val total      = stage.total
    val pipeconfig = stage.getString

    it should f"cos(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}" in {
      test(new MathFunctions(cfg, spec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        counter = 0
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val xi  = generator(spec, r)
          val z0r = reference(xi)
          q += ((xi.value.toBigInt, z0r.value.toBigInt))

          c.io.sel.poke(cfg.signal(Cos))
          c.io.x.poke(xi.value.toBigInt.U(spec.W.W))

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
  // Use the same merged coefficients the hardware sees, so HW == SW bit-exact.
  val cosRefFP32 = SinCosSimGeneric(false, fp32Config.coeffs(Cos), _: RealGeneric)

  val noPipeStage     = new MathFuncPipelineConfig(0, 0, 0, 0, 0, 0)
  val simplePipeStage = new MathFuncPipelineConfig(1, 1, 1, 1, 0, 0)
  val mulPipeStage    = new MathFuncPipelineConfig(2, 1, 1, 2, 1, 1)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, "Test special values", generateSpecialValues(_, _), cosRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    n, "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    n, "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    n, "Test special values", generateSpecialValues(_, _), cosRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, mulPipeStage,
    n, "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, mulPipeStage,
    n, "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, mulPipeStage,
    n, "Test special values", generateSpecialValues(_, _), cosRefFP32)

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = MathFuncConfig.standardBF16
  val cosRefBF16 = SinCosSimGeneric(false, bf16Config.coeffs(Cos), _: RealGeneric)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    n, "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    n, "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    n, "Test special values", generateSpecialValues(_, _), cosRefBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    n, "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    n, "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    n, "Test special values", generateSpecialValues(_, _), cosRefBF16)
}
