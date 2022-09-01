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
import rial.math.FuncKind._
import rial.arith._
import rial.table._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

//
// Testing Exp using ChiselTest
//

class ExpTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test exp"

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
      0.0,
      0.0 + (1.0                    ) * pow(2.0, -1022), // for Double
      0.0 + (1.0 + 1 * pow(2.0, -52)) * pow(2.0, -1022), // for Double
      0.0 + (1.0 + 2 * pow(2.0, -52)) * pow(2.0, -1022), // for Double
      0.0 + (1.0 + 3 * pow(2.0, -52)) * pow(2.0, -1022), // for Double
      0.0 + (1.1                    ) * pow(2.0, -1022), // for Double
      0.5,
      1.0,
      2.0,
      log(2.0),
      log(3.0),
      log(4.0),
    )
  def generateSpecialValues( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx))
  }

  private def runtest ( spec : RealSpec,
      nOrder : Int, adrW : Int, extraBits : Int, stage: MathFuncPipelineConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (RealSpec, Random) => RealGeneric),
      disableTimeout : Boolean = false,
      fncfg: MathFuncConfig = MathFuncConfig.all
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"exp(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${fncfg.getString} " in {
      test( new MathFunctions(fncfg, spec, nOrder, adrW, extraBits, stage, None, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          // To avoid timeoutException while testing z == inf/zero.
          // Detailed explanation follows.
          if(disableTimeout) { c.clock.setTimeout(0) }
          // Detailed explanation:
          //     Here, we disable timeout checking in chiseltest. By default,
          // chiseltest fails after 1000 idle cycles. This is a reasonable
          // behavior in most cases. However, unfortunately, the definition of
          // "an idle cycle" is like: a cycle where input or output does not
          // change. Since we are testing exp function, if we take a value
          // larger than 128 or something, the result becomes inf. So, the test
          // case [127, inf] fails after 1000 cycles. The behavior, returning
          // inf, is the correct, expected behavior. But it is indistinguishable
          // from doing nothing without any calculation. To avoid this false-
          // positive, we manually disable timeout setting if we test large and
          // small values that always return inf and zero, respectively. Note
          // that, assertion and `io.expect` still works so there is no problem
          // in the testing functionality.

          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage
          val reftable   = ExpSim.expTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit) )
          val reference  = ExpSim.expSimGeneric( reftable, _ )

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(fncfg.signal(Exp))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            if(fncfg.has(ATan2Phase1)) {
              c.io.y.get.poke(0.U(spec.W.W))
            }
            val zi = c.io.z.peek().litValue.toBigInt
            if (i > nstage) {

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

              val expx = new RealGeneric(spec, exp(new RealGeneric(spec, xid).toDouble))

              c.io.z.expect(z0d.U(spec.W.W))
              assert(zi == z0d, f"x = ${new RealGeneric(spec, xid).toDouble}(${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
                                f"exp(x) = ${expx.toDouble}(${expx.sgn}|${expx.ex}(${expx.ex - spec.exBias})|${expx.man.toLong.toBinaryString}) " +
                                f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != " +
                                f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
            }
            c.clock.step(1)
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
      PipelineStageConfig.atOut(3),
      PipelineStageConfig.atOut(2),
      PipelineStageConfig.atOut(1),
      true, true, true)


  val nOrderFP32 = 2
  val adrWFP32 = 8
  val extraBitsFP32 = 3

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), /*disableTimeout = */ true)
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), true)

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    "Test Special Values", generateSpecialValues(_,_))


  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), /*disableTimeout = */ true)
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), true)

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    "Test Special Values", generateSpecialValues(_,_))


  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), /*disableTimeout = */ true)
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), true)

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Special Values", generateSpecialValues(_,_))

  val expOnly = new MathFuncConfig(Seq(Exp))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_), false, expOnly)
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_), false, expOnly)

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), /*disableTimeout = */ true, expOnly)
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), true, expOnly)

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_), false, expOnly)
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_), false, expOnly)

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_), false, expOnly)
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_), false, expOnly)

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    "Test Special Values", generateSpecialValues(_,_), false, expOnly)

  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 1

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), /*disableTimeout = */ true)
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), true)

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    "Test Special Values", generateSpecialValues(_,_))


  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), /*disableTimeout = */ true)
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), true)

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    "Test Special Values", generateSpecialValues(_,_))



  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    "Test Safe Positive [1, 127]", generateRealWithin(1.0, 127.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    "Test Safe Negative [-126, -1]", generateRealWithin(-126.0, -1.0,_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    "Test Large Positive [127, inf]", generateRealWithin(127.0, Double.PositiveInfinity,_,_), /*disableTimeout = */ true)
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    "Test Large Negative [-inf, -126]", generateRealWithin(Double.NegativeInfinity, -126.0,_,_), true)

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    "Test Small Positive [2^-7, 1]", generateRealWithin(pow(2.0, -7),1.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    "Test Small Negative [-1, -2^-7]", generateRealWithin(-1.0, -pow(2.0, -7),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    "Test Tiny Positive [0, 2^-7]", generateRealWithin(0.0, pow(2.0, -7),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    "Test Tiny Negative [-2^-7, 0]", generateRealWithin(-pow(2.0, -7), 0.0,_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    "Test Special Values", generateSpecialValues(_,_))


}
