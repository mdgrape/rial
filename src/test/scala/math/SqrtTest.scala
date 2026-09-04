package rial.tests

import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import org.scalatest.funspec.AnyFunSpec

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math._
import rial.math.FuncKind._
import rial.arith._
import rial.util.PipelineStageConfig
import rial.util.ScalaUtil._

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

//
// Testing Sqrt using ChiselTest with new funcs/MathFunctions
//

class SqrtTest extends AnyFunSpec
    with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of random test cases per case (default 1000)"))
  val r = new Random(123456789)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }
  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  var counter = 0
  val specialValues = Seq(
      0.0,
     -0.0,
      1.0,
      2.0,
      4.0,
  )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }

  private def runtest (
    spec:  RealSpec,
    cfg:   MathFuncConfig,
    stage: PipelineStageConfig,
    generatorStr: String, generator: (RealSpec, Random) => RealGeneric,
    reference: RealGeneric => RealGeneric
  ) = {
    val total      = stage.resolve(MathFunctions.stagePoints(cfg)).total
    val pipeconfig = MathTestUtil.label(cfg, stage)

    it(f"sqrt(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new MathFunctions(cfg, spec, stage)) { c =>
        {
          counter = 0

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+total) {
            val xi  = generator(spec,r)
            val z0r = reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))

            c.io.in.valid.poke(true.B)
            c.io.in.bits.sel.poke(cfg.signal(Sqrt))
            c.io.in.bits.x.poke(xi.value.toBigInt.U(spec.W.W))

            val zi = c.io.out.bits.peek().litValue.toBigInt
            if (i > total) {
              assert(c.io.out.valid.peek().litToBoolean,
                f"io.out.valid is not set at cycle $i%d (latency $total%d)")
              val (xid,z0d) = q.dequeue()

              val xidsgn = bit(spec.W-1, xid).toInt
              val xidexp = slice(spec.manW, spec.exW, xid)
              val xidman = xid & maskSL(spec.manW)

              val zisgn = bit(spec.W-1, zi).toInt
              val ziexp = slice(spec.manW, spec.exW, zi)
              val ziman = zi & maskSL(spec.manW)

              val z0dsgn = bit(spec.W-1, z0d).toInt
              val z0dexp = slice(spec.manW, spec.exW, z0d)
              val z0dman = z0d & maskSL(spec.manW)

              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), "+
                                f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != " +
                                f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  // --------------------------------------------------------------------------
  // Test Float32
  // --------------------------------------------------------------------------

  val fp32Config = MathFuncConfig.standardFP32
  val sqrtRefFP32 = SqrtSimGeneric(fp32Config.coeffs(Sqrt), _)

  val noPipeStage     = PipelineStageConfig.none
  val simplePipeStage = PipelineStageConfig.default(4)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), sqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test All range", generateRealFull(_,_), sqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Special Values", generateSpecialValues(_,_), sqrtRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), sqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test All range", generateRealFull(_,_), sqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test Special Values", generateSpecialValues(_,_), sqrtRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), sqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test All range", generateRealFull(_,_), sqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, MathTestUtil.allPoints(fp32Config),
    "Test Special Values", generateSpecialValues(_,_), sqrtRefFP32)

  // --------------------------------------------------------------------------
  // Test BFloat16
  // --------------------------------------------------------------------------

  val bf16Config = MathFuncConfig.standardBF16
  val sqrtRefBF16 = SqrtSimGeneric(bf16Config.coeffs(Sqrt), _)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), sqrtRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test All range", generateRealFull(_,_), sqrtRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Special Values", generateSpecialValues(_,_), sqrtRefBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), sqrtRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test All range", generateRealFull(_,_), sqrtRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test Special Values", generateSpecialValues(_,_), sqrtRefBF16)

  // --------------------------------------------------------------------------
  // simple config (Sqrt / InvSqrt / Reciprocal only -- no multiplier).
  // Re-run on top of the standard sweeps to cover the no-multiplier build.
  // --------------------------------------------------------------------------

  val fp32ConfigSimple  = MathFuncConfig.simpleFP32
  val sqrtRefFP32Simple = SqrtSimGeneric(fp32ConfigSimple.coeffs(Sqrt), _)

  runtest(RealSpec.Float32Spec, fp32ConfigSimple, noPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), sqrtRefFP32Simple)
  runtest(RealSpec.Float32Spec, fp32ConfigSimple, noPipeStage,
    "Test All range", generateRealFull(_,_), sqrtRefFP32Simple)
  runtest(RealSpec.Float32Spec, fp32ConfigSimple, noPipeStage,
    "Test Special Values", generateSpecialValues(_,_), sqrtRefFP32Simple)

  runtest(RealSpec.Float32Spec, fp32ConfigSimple, MathTestUtil.allPoints(fp32ConfigSimple),
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), sqrtRefFP32Simple)
  runtest(RealSpec.Float32Spec, fp32ConfigSimple, MathTestUtil.allPoints(fp32ConfigSimple),
    "Test All range", generateRealFull(_,_), sqrtRefFP32Simple)
  runtest(RealSpec.Float32Spec, fp32ConfigSimple, MathTestUtil.allPoints(fp32ConfigSimple),
    "Test Special Values", generateSpecialValues(_,_), sqrtRefFP32Simple)

  val bf16ConfigSimple  = MathFuncConfig.simpleBF16
  val sqrtRefBF16Simple = SqrtSimGeneric(bf16ConfigSimple.coeffs(Sqrt), _)

  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, noPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), sqrtRefBF16Simple)
  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, noPipeStage,
    "Test All range", generateRealFull(_,_), sqrtRefBF16Simple)
  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, noPipeStage,
    "Test Special Values", generateSpecialValues(_,_), sqrtRefBF16Simple)

  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, simplePipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), sqrtRefBF16Simple)
  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, simplePipeStage,
    "Test All range", generateRealFull(_,_), sqrtRefBF16Simple)
  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, simplePipeStage,
    "Test Special Values", generateSpecialValues(_,_), sqrtRefBF16Simple)
}
