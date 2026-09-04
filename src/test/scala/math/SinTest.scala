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

class SinTest extends AnyFunSpec
    with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of random test cases per case (default 1000)"))
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
    stage: PipelineStageConfig,
    generatorStr: String, generator: (RealSpec, Random) => RealGeneric,
    reference: RealGeneric => RealGeneric,
  ) = {
    val total      = stage.resolve(MathFunctions.stagePoints(cfg)).total
    val pipeconfig = MathTestUtil.label(cfg, stage)

    it(f"sin(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new MathFunctions(cfg, spec, stage)) { c =>
        counter = 0
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val xi  = generator(spec, r)
          val z0r = reference(xi)
          q += ((xi.value.toBigInt, z0r.value.toBigInt))

          c.io.in.valid.poke(true.B)
          c.io.in.bits.sel.poke(cfg.signal(Sin))
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
  // Use the same merged coefficients the hardware sees, so HW == SW bit-exact.
  val sinRefFP32 = SinCosSimGeneric(true, fp32Config.coeffs(Sin), _: RealGeneric)

  // pipeline configurations -- exercise pre/post multipliers at 0 and 1 cycles
  val noPipeStage     = PipelineStageConfig.none
  val simplePipeStage = PipelineStageConfig.default(4)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), sinRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), sinRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), sinRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), sinRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), sinRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test special values", generateSpecialValues(_, _), sinRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), sinRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), sinRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test special values", generateSpecialValues(_, _), sinRefFP32)

  // one register inside the sin/cos argument reduction, one at every path of
  // the result assembly
  val sinCosPipeStage =
    PipelineStageConfig.specify("pre.sinCos.fold" -> 1, "post" -> 1)

  runtest(RealSpec.Float32Spec, fp32Config, sinCosPipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), sinRefFP32)

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = MathFuncConfig.standardBF16
  val sinRefBF16 = SinCosSimGeneric(true, bf16Config.coeffs(Sin), _: RealGeneric)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), sinRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), sinRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), sinRefBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), sinRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), sinRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test special values", generateSpecialValues(_, _), sinRefBF16)
}
