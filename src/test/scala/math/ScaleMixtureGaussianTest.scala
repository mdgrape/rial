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
// Testing ScaleMixtureGaussian using ChiselTest
//

class ScaleMixtureGaussianTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test scaleMixtureGaussian"

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
      1.0,
      math.E / 2.0,
      math.E,
      math.E * 2.0,
      pow(math.E, 2),
      pow(math.E, 3),
      pow(math.E, 4),
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

  private def runtest (sgmA: Double, sgmB: Double,
      spec : RealSpec,
      nOrder : Int, adrW : Int, extraBits : Int, stage: MathFuncPipelineConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (RealSpec, Random) => RealGeneric),
      disableTimeout : Boolean = false,
      fncfg: MathFuncConfig = MathFuncConfig.all
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"SMG ($pipeconfig, ${spec.toStringShort}) $generatorStr with ${fncfg.getString} " in {
      test( new MathFunctions(fncfg, spec, nOrder, adrW, extraBits, stage, None, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          import FuncKind._

          counter = 0
          // since table result depends on these values, it is unavoidable to
          // construct tables multiple times.
          val maxCbits   = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW

          val scaleMixtureGaussianTable = ScaleMixtureGaussianSim.tableGeneration(
           nOrder, adrW, spec.manW, spec.manW+extraBits, sgmA, sgmB, Some(maxCalcW), Some(maxCbits))

          val reference = ScaleMixtureGaussianSim.scaleMixtureGaussianSimGeneric(
            scaleMixtureGaussianTable, _, sgmA, sgmB, false)

          // To avoid timeoutException while testing z == neg.
          // Detailed explanation follows.
          if(disableTimeout) { c.clock.setTimeout(0) }
          // Detailed explanation:
          //     Here, we disable timeout checking in chiseltest. By default,
          // chiseltest fails after 1000 idle cycles. This is a reasonable
          // behavior in most cases. However, unfortunately, the definition of
          // "an idle cycle" is like: a cycle where input or output does not
          // change. Since we are testing scaleMixtureGaussian function, if we take a negative
          // value, the result becomes nan. So, the test case [-inf, 0] fails
          // after 1000 cycles. The behavior, returning nan, is the correct,
          // expected behavior. But it is indistinguishable from doing nothing
          // without any calculation. To avoid this false positive, we manually
          // disable timeout setting if we test negative values that always
          // return nan. Note that, assertion and `io.expect` still works so
          // there is no problem in the testing functionality.

          val nstage = c.getStage
          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
//             println("=========================================================")
            val xi = generator(spec,r)
            val z0r= reference(xi)

            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(fncfg.signal(ScaleMixtureGaussian))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            if(fncfg.has(ATan2Phase1)) {
              c.io.y.get.poke(0.U(spec.W.W))
            }

            val zi = c.io.z.peek().litValue.toBigInt
            c.clock.step(1)
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

              assert(zi == z0d, f"x = ${new RealGeneric(spec, xidsgn, xidexp.toInt, xidman).toDouble}(${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
                                f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != " +
                                f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
            }
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

  val sgmA = MathFuncConfig.all.scaleMixtureGaussianSigma.get._1 // exp(-1.0)
  val sgmB = MathFuncConfig.all.scaleMixtureGaussianSigma.get._2 // exp(-6.0)

  val nOrderFP32 = 2
  val adrWFP32 = 8
  val extraBitsFP32 = 8

  val scaleMixtureGaussianOnly = new MathFuncConfig(Seq(ScaleMixtureGaussian), Some((sgmA, sgmB)))

  val transitionPointF32 = pow(2.0, ScaleMixtureGaussianSim.tableDomainDigit(23, sgmA, sgmB))

  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    f"Test x Larger Than Pt [${transitionPointF32}, inf)",
    generateRealWithin(transitionPointF32, pow(2.0, 127.0),_,_))
  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    f"Test x Smaller Than Pt (-${transitionPointF32}, ${transitionPointF32})",
    generateRealWithin(-transitionPointF32, transitionPointF32,_,_))
  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r,
    f"Test -x Larger Than Pt [${transitionPointF32}, inf)",
    generateRealWithin(-pow(2.0, 127.0), -transitionPointF32,_,_))

  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    f"Test x Larger Than Pt [${transitionPointF32}, inf)",
    generateRealWithin(transitionPointF32, pow(2.0, 127.0),_,_))
  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    f"Test x Smaller Than Pt [0, ${transitionPointF32})",
    generateRealWithin(-transitionPointF32, transitionPointF32,_,_))
  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r,
    f"Test -x Larger Than Pt [${transitionPointF32}, inf)",
    generateRealWithin(-pow(2.0, 127.0), -transitionPointF32,_,_))

  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    f"Test x Larger Than Pt [${transitionPointF32}, inf)",
    generateRealWithin(transitionPointF32, pow(2.0, 127.0),_,_))
  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    f"Test x Smaller Than Pt [0, ${transitionPointF32})",
    generateRealWithin(-transitionPointF32, transitionPointF32,_,_))
  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    f"Test -x Larger Than Pt [${transitionPointF32}, inf)",
    generateRealWithin(-pow(2.0, 127.0), -transitionPointF32,_,_))

  // no other math func

  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    f"Test x Larger Than Pt [${transitionPointF32}, inf)",
    generateRealWithin(transitionPointF32, pow(2.0, 127.0),_,_), false, scaleMixtureGaussianOnly)
  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    f"Test x Smaller Than Pt [0, ${transitionPointF32})",
    generateRealWithin(-transitionPointF32, transitionPointF32,_,_), false, scaleMixtureGaussianOnly)
  runtest(sgmA, sgmB, RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r,
    f"Test -x Larger Than Pt [${transitionPointF32}, inf)",
    generateRealWithin(-pow(2.0, 127.0), -transitionPointF32,_,_), false, scaleMixtureGaussianOnly)

  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 10

  val transitionPointBF16 = pow(2.0, ScaleMixtureGaussianSim.tableDomainDigit(7, sgmA, sgmB))

  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    f"Test x Larger Than Pt [${transitionPointBF16}, inf)",
    generateRealWithin(transitionPointBF16, pow(2.0, 127.0),_,_))
  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    f"Test x Smaller Than Pt (-${transitionPointBF16}, ${transitionPointBF16})",
    generateRealWithin(-transitionPointBF16, transitionPointBF16,_,_))
  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r,
    f"Test -x Larger Than Pt [${transitionPointBF16}, inf)",
    generateRealWithin(-pow(2.0, 127.0), -transitionPointBF16,_,_))

  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    f"Test x Larger Than Pt [${transitionPointBF16}, inf)",
    generateRealWithin(transitionPointBF16, pow(2.0, 127.0),_,_))
  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    f"Test x Smaller Than Pt [0, ${transitionPointBF16})",
    generateRealWithin(-transitionPointBF16, transitionPointBF16,_,_))
  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r,
    f"Test -x Larger Than Pt [${transitionPointBF16}, inf)",
    generateRealWithin(-pow(2.0, 127.0), -transitionPointBF16,_,_))

  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    f"Test x Larger Than Pt [${transitionPointBF16}, inf)",
    generateRealWithin(transitionPointBF16, pow(2.0, 127.0),_,_))
  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    f"Test x Smaller Than Pt [0, ${transitionPointBF16})",
    generateRealWithin(-transitionPointBF16, transitionPointBF16,_,_))
  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    f"Test -x Larger Than Pt [${transitionPointBF16}, inf)",
    generateRealWithin(-pow(2.0, 127.0), -transitionPointBF16,_,_))

  // no other math func

  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    f"Test x Larger Than Pt [${transitionPointBF16}, inf)",
    generateRealWithin(transitionPointBF16, pow(2.0, 127.0),_,_), false, scaleMixtureGaussianOnly)
  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    f"Test x Smaller Than Pt [0, ${transitionPointBF16})",
    generateRealWithin(-transitionPointBF16, transitionPointBF16,_,_), false, scaleMixtureGaussianOnly)
  runtest(sgmA, sgmB, RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r,
    f"Test -x Larger Than Pt [${transitionPointBF16}, inf)",
    generateRealWithin(-pow(2.0, 127.0), -transitionPointBF16,_,_), false, scaleMixtureGaussianOnly)
}
