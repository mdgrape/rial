package rial.tests

import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import org.scalatest.funspec.AnyFunSpec

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math._
import rial.math.FuncKind._
import rial.arith._
import rial.util.ScalaUtil._
import rial.testUtil.RialChiselSim

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

//
// Testing InvSqrt using ChiselTest with new funcs/MathFunctions
//

class InvSqrtTest extends AnyFunSpec
    with RialChiselSim with EmitVcd {

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
    stage: MathFuncPipelineConfig,
    generatorStr: String, generator: (RealSpec, Random) => RealGeneric,
    reference: RealGeneric => RealGeneric
  ) = {
    val total      = stage.total
    val pipeconfig = stage.getString

    it(f"invsqrt(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${cfg.getString}") {
      val n = getOption[Int]("n").getOrElse(1000)
      simulate(new MathFunctions(cfg, spec, stage)) { c =>
        {
          counter = 0

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+total) {
            val xi  = generator(spec,r)
            val z0r = reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))

            c.io.sel.poke(cfg.signal(InvSqrt))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))

            val zi = c.io.z.peek().litValue.toBigInt
            if (i > total) {
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
  val invsqrtRefFP32 = InvSqrtSimGeneric(fp32Config.coeffs(InvSqrt), _)

  val noPipeStage = new MathFuncPipelineConfig(0,0,0,0)
  val simplePipeStage = new MathFuncPipelineConfig(1,1,1,1)
  val complexPipeStage = new MathFuncPipelineConfig(1,1,3,1)

  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), invsqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test All range", generateRealFull(_,_), invsqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, noPipeStage,
    "Test Special Values", generateSpecialValues(_,_), invsqrtRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), invsqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test All range", generateRealFull(_,_), invsqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, simplePipeStage,
    "Test Special Values", generateSpecialValues(_,_), invsqrtRefFP32)

  runtest(RealSpec.Float32Spec, fp32Config, complexPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), invsqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, complexPipeStage,
    "Test All range", generateRealFull(_,_), invsqrtRefFP32)
  runtest(RealSpec.Float32Spec, fp32Config, complexPipeStage,
    "Test Special Values", generateSpecialValues(_,_), invsqrtRefFP32)

  // --------------------------------------------------------------------------
  // Test BFloat16
  // --------------------------------------------------------------------------

  val bf16Config = MathFuncConfig.standardBF16
  val invsqrtRefBF16 = InvSqrtSimGeneric(bf16Config.coeffs(InvSqrt), _)

  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), invsqrtRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test All range", generateRealFull(_,_), invsqrtRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, noPipeStage,
    "Test Special Values", generateSpecialValues(_,_), invsqrtRefBF16)

  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), invsqrtRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test All range", generateRealFull(_,_), invsqrtRefBF16)
  runtest(RealSpec.BFloat16Spec, bf16Config, simplePipeStage,
    "Test Special Values", generateSpecialValues(_,_), invsqrtRefBF16)

  // --------------------------------------------------------------------------
  // simple config (Sqrt / InvSqrt / Reciprocal only -- no multiplier).
  // --------------------------------------------------------------------------

  val fp32ConfigSimple     = MathFuncConfig.simpleFP32
  val invsqrtRefFP32Simple = InvSqrtSimGeneric(fp32ConfigSimple.coeffs(InvSqrt), _)

  runtest(RealSpec.Float32Spec, fp32ConfigSimple, noPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), invsqrtRefFP32Simple)
  runtest(RealSpec.Float32Spec, fp32ConfigSimple, noPipeStage,
    "Test All range", generateRealFull(_,_), invsqrtRefFP32Simple)
  runtest(RealSpec.Float32Spec, fp32ConfigSimple, noPipeStage,
    "Test Special Values", generateSpecialValues(_,_), invsqrtRefFP32Simple)

  runtest(RealSpec.Float32Spec, fp32ConfigSimple, complexPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), invsqrtRefFP32Simple)
  runtest(RealSpec.Float32Spec, fp32ConfigSimple, complexPipeStage,
    "Test All range", generateRealFull(_,_), invsqrtRefFP32Simple)
  runtest(RealSpec.Float32Spec, fp32ConfigSimple, complexPipeStage,
    "Test Special Values", generateSpecialValues(_,_), invsqrtRefFP32Simple)

  val bf16ConfigSimple     = MathFuncConfig.simpleBF16
  val invsqrtRefBF16Simple = InvSqrtSimGeneric(bf16ConfigSimple.coeffs(InvSqrt), _)

  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, noPipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), invsqrtRefBF16Simple)
  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, noPipeStage,
    "Test All range", generateRealFull(_,_), invsqrtRefBF16Simple)
  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, noPipeStage,
    "Test Special Values", generateSpecialValues(_,_), invsqrtRefBF16Simple)

  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, simplePipeStage,
    "Test Within (-128,128)", generateRealWithin(128.0,_,_), invsqrtRefBF16Simple)
  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, simplePipeStage,
    "Test All range", generateRealFull(_,_), invsqrtRefBF16Simple)
  runtest(RealSpec.BFloat16Spec, bf16ConfigSimple, simplePipeStage,
    "Test Special Values", generateSpecialValues(_,_), invsqrtRefBF16Simple)
}
