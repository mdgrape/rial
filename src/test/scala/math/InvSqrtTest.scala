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
// Testing InvSqrt using ChiselTest
//

class InvSqrtTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test invsqrt"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

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
      8.0,
      16.0,
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
    it should f"invsqrt(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr with ${fncfg.getString} " in {
      test( new MathFunctions(fncfg, spec, nOrder, adrW, extraBits, stage, None, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          counter = 0
          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage
          val invsqrtTable = InvSqrtSim.invsqrtTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit) )
          val reference  = InvSqrtSim.invsqrtSimGeneric( invsqrtTable, _ )

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(fncfg.signal(InvSqrt))
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

              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
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

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test All range",generateRealFull(_,_) )
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test SpecialValue",generateSpecialValues(_,_) )

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test All range",generateRealFull(_,_) )
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test SpecialValue",generateSpecialValues(_,_) )

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test All range",generateRealFull(_,_) )
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test SpecialValue",generateSpecialValues(_,_) )

  val invsqrtOnly = new MathFuncConfig(Seq(InvSqrt))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_), invsqrtOnly)
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test All range",generateRealFull(_,_), invsqrtOnly )
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test SpecialValue",generateSpecialValues(_,_), invsqrtOnly )



  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 1

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test All range",generateRealFull(_,_) )
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test SpecialValue",generateSpecialValues(_,_) )

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test All range",generateRealFull(_,_) )
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test SpecialValue",generateSpecialValues(_,_) )

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test All range",generateRealFull(_,_) )
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test SpecialValue",generateSpecialValues(_,_) )
}
