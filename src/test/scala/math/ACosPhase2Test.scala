import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math._
import rial.arith._
import rial.table._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

//
// Testing ACos using ChiselTest
//

class ACosPhase2Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test acos stage2"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRealWithin( from : Double, to : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = from + r.nextDouble() * (to - from)
    new RealGeneric(spec, rD)
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  var counter = 0
  val specialValues = Seq(
      0.0,                -0.0,
      1.0,                -1.0,
      1.0 / 2.0,          -1.0 / 2.0,
      1.1,                -1.1, // should be acos(+/-1.0)
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  private def runtest ( spec : RealSpec,
      nOrder : Int, adrW : Int, extraBits : Int, stage: MathFuncPipelineConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (RealSpec, Random) => RealGeneric),
      fncfg: MathFuncConfig = MathFuncConfig.all
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"acos(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(fncfg, spec, nOrder, adrW, extraBits, stage, None, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          import FuncKind._

          counter = 0
          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage
          val sqrtF32TableI = SqrtSim.sqrtTableGeneration(nOrder, adrW,
            spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))
          val acosF32TableI = ACosPhase2Sim.acosTableGeneration(nOrder, adrW,
            spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))

          val refacos1  = ACosPhase1Sim.acosPhase1SimGeneric(sqrtF32TableI, _ )
          val refacos2  = ACosPhase2Sim.acosPhase2SimGeneric(acosF32TableI, _ )

          for(i <- 1 to n) {
            val xi  = generator(spec,r)
            val z0r = refacos1(xi)
            val z1r = refacos2(z0r)

            c.io.sel.poke(fncfg.signal(ACosPhase1))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            if(fncfg.has(ATan2Phase1)) {
              c.io.y.get.poke(0.U(spec.W.W))
            }

            if(stage.total > 0) {
              c.clock.step(1)
              c.io.sel.poke(fncfg.signalNone())
              c.io.x.poke(0.U(spec.W.W))
              for(j <- 1 until max(1, stage.total)) {
                c.clock.step(1)
              }
            }

            val z0i = new RealGeneric(spec, c.io.z.peek().litValue.toBigInt)
            assert(z0i.value == z0r.value,
              f"x = (${xi.sgn}|${xi.ex}(${xi.ex - spec.exBias})|${xi.man.toLong.toBinaryString})(${xi.toDouble}), " +
              f"test(${z0i.sgn}|${z0i.ex}(${z0i.ex - spec.exBias})|${z0i.man.toLong.toBinaryString})(${z0i.toDouble}) != " +
              f"ref (${z0r.sgn}|${z0r.ex}(${z0r.ex - spec.exBias})|${z0r.man.toLong.toBinaryString})(${z0r.toDouble}) should be " +
              f"sqrt(1-|${xi.toDouble}|) = ${sqrt(1.0 - abs(xi.toDouble))}")

            if(stage.total == 0) {
              c.clock.step(1)
            }

            c.io.sel.poke(fncfg.signal(ACosPhase2))
            c.io.x.poke(z0i.value.toBigInt.U(spec.W.W))

            if(stage.total > 0) {
              c.clock.step(1)
              c.io.sel.poke(fncfg.signalNone())
              c.io.x.poke(0.U(spec.W.W))
              for(j <- 1 until max(1, stage.total)) {
                c.clock.step(1)
              }
            }

            val z1i = new RealGeneric(spec, c.io.z.peek().litValue.toBigInt)
            c.clock.step(1)

            assert(z1i.value == z1r.value,
              f"\nx     = (${xi.sgn}|${xi.ex}(${xi.ex - spec.exBias})|${xi.man.toLong.toBinaryString})(${xi.toDouble}), " +
              f"\ntest1 = (${z0i.sgn}|${z0i.ex}(${z0i.ex - spec.exBias})|${z0i.man.toLong.toBinaryString})(${z0i.toDouble}), " +
              f"\nref1  = (${z0r.sgn}|${z0r.ex}(${z0r.ex - spec.exBias})|${z0r.man.toLong.toBinaryString})(${z0r.toDouble}), " +
              f"\ntest2 = (${z1i.sgn}|${z1i.ex}(${z1i.ex - spec.exBias})|${z1i.man.toLong.toBinaryString})(${z1i.toDouble}) != " +
              f"\nref2  = (${z1r.sgn}|${z1r.ex}(${z1r.ex - spec.exBias})|${z1r.man.toLong.toBinaryString})(${z1r.toDouble}) should be " +
              f"\nacos(${xi.toDouble}) = ${acos(xi.toDouble)}")
          }
        }
      }
    }
  }
  val simplePipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      true, true, true)
  val complexPipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(2),
      PipelineStageConfig.atOut(2),
      PipelineStageConfig.atOut(1),
      true, true, true)

  val nOrderFP32 = 2
  val adrWFP32 = 8
  val extraBitsFP32 = 3
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within (-1+2^-4,  -0.5)",    generateRealWithin(-1.0+pow(2.0, -4), -0.5,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within (-0.5,  -2^-4)",    generateRealWithin(-0.5, -pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 2^-4,  0.5)",     generateRealWithin(pow(2.0, -4), 0.5,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 0.5,   1-2^-4)",     generateRealWithin(0.5, 1.0 - pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test special value", generateSpecialValues(_,_))


  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test special value", generateSpecialValues(_,_))


  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test special value", generateSpecialValues(_,_))

  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 1

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-1+2^-4,  -0.5)",    generateRealWithin(-1.0+pow(2.0, -4), -0.5,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-0.5,  -2^-4)",    generateRealWithin(-0.5, -pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 2^-4,  0.5)",     generateRealWithin(pow(2.0, -4), 0.5,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 0.5,   1-2^-4)",     generateRealWithin(0.5, 1.0 - pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test special value", generateSpecialValues(_,_))


  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test special value", generateSpecialValues(_,_))


  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test special value", generateSpecialValues(_,_))
}
