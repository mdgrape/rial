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

class CosTest extends AnyFunSpec
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

  var underflowCounter = 0
  /** |x| in the bands where |x|/pi underflows, and just above them. */
  def generateUnderflow(spec: RealSpec, r: Random) = {
    val ex = 1 + (underflowCounter % 4)
    underflowCounter += 1
    new RealGeneric(spec, (SafeLong(r.nextInt(2)) << (spec.W - 1)) +
                          (SafeLong(ex) << spec.manW) +
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

    it(f"cos(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new MathFunctions(cfg, stage)) { c =>
        counter = 0
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val xi  = generator(spec, r)
          val z0r = reference(xi)
          q += ((xi.value.toBigInt, z0r.value.toBigInt))

          c.io.in.valid.poke(true.B)
          c.io.in.bits.sel.poke(cfg.signal(Cos))
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

  val fp32Config = MathFuncConfig.standardFP32
  // Use the same merged coefficients the hardware sees, so HW == SW bit-exact.
  val cosRefFP32 = SinCosSimGeneric(false, fp32Config.coeffs(Cos), _: RealGeneric)

  val noPipeStage     = PipelineStageConfig.none
  val simplePipeStage = PipelineStageConfig.default(4)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), cosRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test special values", generateSpecialValues(_, _), cosRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test special values", generateSpecialValues(_, _), cosRefFP32)

  // one register inside the sin/cos argument reduction, one at every path of
  // the result assembly
  val sinCosPipeStage =
    PipelineStageConfig.specify("pre.sinCos.fold" -> 1, "post" -> 1)

  runtest(RealSpec.Float32Spec, fp32Config, sinCosPipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefFP32)

  // --------------------------------------------------------------------------
  // the smallest arguments, where |x|/pi underflows

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test the smallest exponents", generateUnderflow(_, _), cosRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test the smallest exponents", generateUnderflow(_, _), cosRefFP32)

  // --------------------------------------------------------------------------
  // the same arguments, with the reduction restricted to their range

  val fp32TwoPi = fp32Config.restrictSinCos(
    Some(SinCosDomain.twoPi(RealSpec.Float32Spec)))
  val fp32Pi = fp32Config.restrictSinCos(
    Some(SinCosDomain.pi(RealSpec.Float32Spec)))
  val cosRefTwoPi = SinCosSimGeneric(false, fp32TwoPi.coeffs(Cos), _: RealGeneric)
  val cosRefPi    = SinCosSimGeneric(false, fp32Pi.coeffs(Cos),    _: RealGeneric)

  runtest(RealSpec.Float32Spec, fp32TwoPi, noPipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefTwoPi)
  runtest(RealSpec.Float32Spec, fp32TwoPi, noPipeStage,
    "Test special values", generateSpecialValues(_, _), cosRefTwoPi)
  runtest(RealSpec.Float32Spec, fp32TwoPi, noPipeStage,
    "Test the smallest exponents", generateUnderflow(_, _), cosRefTwoPi)
  runtest(RealSpec.Float32Spec, fp32TwoPi, MathTestUtil.allPoints(fp32TwoPi),
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefTwoPi)

  runtest(RealSpec.Float32Spec, fp32Pi, noPipeStage,
    "Test Within (-pi, pi)", generateRealWithin(2 * Pi, _, _), cosRefPi)
  runtest(RealSpec.Float32Spec, fp32Pi, noPipeStage,
    "Test the smallest exponents", generateUnderflow(_, _), cosRefPi)
  runtest(RealSpec.Float32Spec, fp32Pi, MathTestUtil.allPoints(fp32Pi),
    "Test Within (-pi, pi)", generateRealWithin(2 * Pi, _, _), cosRefPi)

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = MathFuncConfig.standardBF16
  val cosRefBF16 = SinCosSimGeneric(false, bf16Config.coeffs(Cos), _: RealGeneric)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test special values", generateSpecialValues(_, _), cosRefBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test Within (-2pi, 2pi)", generateRealWithin(4 * Pi, _, _), cosRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test Within (-128, 128)", generateRealWithin(128.0, _, _), cosRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test special values", generateSpecialValues(_, _), cosRefBF16)
}
