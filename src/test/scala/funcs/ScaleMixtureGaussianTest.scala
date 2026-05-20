package rial.tests.funcs

import org.scalatest._

import chisel3._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.implicits._

import rial.funcs._
import rial.funcs.FuncKind._
import rial.arith._
import rial.util.ScalaUtil._

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue

class ScaleMixtureGaussianTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test scaleMixtureGaussian with funcs/MathFunctions"

  val n = 1000
  val r = new Random(123456789)

  def generateRealWithin(p: Double, spec: RealSpec, r: Random) = {
    val rD: Double = (r.nextDouble() - 0.5) * 2.0 * p
    val x = new RealGeneric(spec, rD)
    new RealGeneric(spec, (x.value & (maskSL(spec.exW + 1) << spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull(spec: RealSpec, r: Random) = {
    new RealGeneric(spec, SafeLong(BigInt(spec.W, r)))
  }

  var counter = 0
  val specialValues = Seq(0.0, -0.0)
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
    n: Int, generatorStr: String,
    generator: (RealSpec, Random) => RealGeneric,
  ) = {
    val total      = stage.total
    val pipeconfig = stage.getString
    val (sgmA, sgmB) = cfg.scaleMixtureGaussianSigma.get

    val reference = ScaleMixtureGaussianSimGeneric(
      cfg.coeffs(ScaleMixtureGaussian), _: RealGeneric, sgmA, sgmB)

    it should f"smg(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}" in {
      test(new MathFunctions(cfg, spec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        counter = 0
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val xi  = generator(spec, r)
          val z0r = reference(xi)
          q += ((xi.value.toBigInt, z0r.value.toBigInt))

          c.io.sel.poke(cfg.signal(ScaleMixtureGaussian))
          c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
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

            assert(zi == z0d,
              f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
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

  val sgmA16 = exp(-1.0)
  val sgmB16 = exp(-6.0)

  val fp32Config = new MathFuncConfig(
    funcs   = Seq(ScaleMixtureGaussian),
    nOrder  = 2,
    addrW   = 8,
    inBitW  = RealSpec.Float32Spec.manW - 8,
    outBitW = RealSpec.Float32Spec.manW + 8,
    scaleMixtureGaussianSigma = Some((sgmA16, sgmB16)),
  )

  val transitionPointF32 = pow(2.0,
    ScaleMixtureGaussianTable.tableDomainDigit(fp32Config.outBitW, sgmA16, sgmB16))

  val noPipeStage     = new MathFuncPipelineConfig(0, 0, 0, 0, 0, 0)
  val simplePipeStage = new MathFuncPipelineConfig(1, 1, 1, 1, 0, 0)
  val complexPipeline = new MathFuncPipelineConfig(2, 1, 1, 2, 1, 1)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, f"Test Within (-${transitionPointF32},${transitionPointF32}) where table is used",
    generateRealWithin(transitionPointF32, _, _))
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, "Test All range", generateRealFull(_, _))
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    n, "Test special values", generateSpecialValues(_, _))

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    n, "Test All range", generateRealFull(_, _))
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    n, "Test special values", generateSpecialValues(_, _))

  runtest(RealSpec.Float32Spec, fp32Config, complexPipeline,
    n, "Test All range", generateRealFull(_, _))
  runtest(RealSpec.Float32Spec, fp32Config, complexPipeline,
    n, "Test special values", generateSpecialValues(_, _))

  // --------------------------------------------------------------------------
  // BFloat16

  val bf16Config = new MathFuncConfig(
    funcs   = Seq(ScaleMixtureGaussian),
    nOrder  = 0,
    addrW   = RealSpec.BFloat16Spec.manW,
    inBitW  = 0,
    outBitW = RealSpec.BFloat16Spec.manW + 10,
    scaleMixtureGaussianSigma = Some((sgmA16, sgmB16)),
  )

  val transitionPointBF16 = pow(2.0,
    ScaleMixtureGaussianTable.tableDomainDigit(bf16Config.outBitW, sgmA16, sgmB16))

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    n, f"Test Within (-${transitionPointBF16},${transitionPointBF16}) where table is used",
    generateRealWithin(transitionPointBF16, _, _))
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    n, "Test All range", generateRealFull(_, _))

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    n, "Test All range", generateRealFull(_, _))
}
