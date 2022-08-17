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

class ACosStage1Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test acos stage1"

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
      generator : ( (RealSpec, Random) => RealGeneric)
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"acos(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(spec, nOrder, adrW, extraBits, stage, None, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          counter = 0
          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage
          val sqrtF32TableI = SqrtSim.sqrtTableGeneration(nOrder, adrW,
            spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))

          val reference  = ACosStage1Sim.acosStage1SimGeneric(sqrtF32TableI, _ )
          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi  = generator(spec,r)
            val (z0r, xneg, special) = reference(xi)

            if(xi.toDouble < 0.0) {assert(xneg)} else {assert(!xneg)}
            if     (xi.toDouble ==  0.0) {assert(special == 0)}
            else if(xi.toDouble ==  1.0) {assert(special == 1)}
            else if(xi.toDouble == -1.0) {assert(special == 1)}

            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(SelectFunc.ACos1)
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            c.io.y.poke(0.U(spec.W.W))
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

              val rxid = new RealGeneric(spec, xid)
              val rz0d = new RealGeneric(spec, z0d)
              val rzi  = new RealGeneric(spec, zi)

              if (zi != z0d) {
                println(f"x = ${rxid.sgn}|${rxid.ex}(${rxid.ex - spec.exBias})|${rxid.man}(${rxid.toDouble})")
                println(f"z = ${rz0d.sgn}|${rz0d.ex}(${rz0d.ex - spec.exBias})|${rz0d.man}(${rz0d.toDouble}), should be ${sqrt(1.0 - rxid.toDouble)}")
                println(f"c = ${rzi .sgn}|${rzi .ex}(${rzi .ex - spec.exBias})|${rzi .man}(${rzi .toDouble})")
              }
              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
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
      true, true, true)
  val complexPipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(2),
      PipelineStageConfig.atOut(2),
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


//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
//     n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
//     n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
//     n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
//     n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
//     n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
//     n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
//     n, r, "Test special value", generateSpecialValues(_,_))
// 
// 
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
//     n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
//     n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
//     n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
//     n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
//     n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
//     n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
//     n, r, "Test special value", generateSpecialValues(_,_))
// 
//   val nOrderBF16 = 0
//   val adrWBF16 = 7
//   val extraBitsBF16 = 1
// 
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
//     n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
//     n, r, "Test Within (-1+2^-4,  -0.5)",    generateRealWithin(-1.0+pow(2.0, -4), -0.5,_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
//     n, r, "Test Within (-0.5,  -2^-4)",    generateRealWithin(-0.5, -pow(2.0, -4),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
//     n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
//     n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
//     n, r, "Test Within ( 2^-4,  0.5)",     generateRealWithin(pow(2.0, -4), 0.5,_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
//     n, r, "Test Within ( 0.5,   1-2^-4)",     generateRealWithin(0.5, 1.0 - pow(2.0, -4),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
//     n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
//     n, r, "Test special value", generateSpecialValues(_,_))
// 
// 
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
//     n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
//     n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
//     n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
//     n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
//     n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
//     n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
//     n, r, "Test special value", generateSpecialValues(_,_))
// 
// 
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
//     n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
//     n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
//     n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
//     n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
//     n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
//     n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
//   runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
//     n, r, "Test special value", generateSpecialValues(_,_))
}
