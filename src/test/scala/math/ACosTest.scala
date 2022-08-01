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

class ACosTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test acos"

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
//       sqrt(2.0) / 2.0,    -sqrt(2.0) / 2.0,
//       sqrt(3.0) / 2.0,    -sqrt(3.0) / 2.0,
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
          val acosF32TableI = ACosSim.acosTableGeneration(spec, nOrder, adrW,
            spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))
          val sqrtF32TableI = ACosSim.sqrtTableGeneration(nOrder, adrW,
            spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))
          val acosExTableI  = if(nOrder != 0) {None} else {
            Some(ACosSim.acosExTableGeneration(spec, nOrder, adrW,
              spec.manW, spec.manW+extraBits))
          }

          val reference  = ACosSim.acosSimGeneric(acosF32TableI, sqrtF32TableI, _, acosExTableI )
          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(SelectFunc.ACos)
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            c.io.y.poke(0.U(spec.W.W))
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

              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
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


// =============================================================================

class ACosOnlyTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test acos only"

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
      test( new ACosGeneric(spec, nOrder, adrW, extraBits, stage, None, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
//           println("----------------------------")
          val maxCbit    = c.getCbit
          val maxCalcW   = c.getCalcW
          val nstage     = c.getStage
          val acosTableI = ACosSim.acosTableGeneration(spec, nOrder, adrW,
            spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))
          val sqrtTableI = ACosSim.sqrtTableGeneration(nOrder, adrW,
            spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))

          val acosExTableI = if(nOrder != 0) {None} else {
            Some(ACosSim.acosExTableGeneration(spec, nOrder, adrW,
              spec.manW, spec.manW+extraBits))
          }

          val reference  = ACosSim.acosSimGeneric(acosTableI, sqrtTableI, _, acosExTableI )
          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.en.poke(true.B)
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
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

              val x = new RealGeneric(spec, xidsgn, xidexp.toInt, xidman)
              val z = new RealGeneric(spec, acos(x.toDouble))
              val zsgn = z.sgn
              val zexp = z.ex
              val zman = z.man

              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString})(${x.toDouble}), " +
                                f"test(${zisgn }|${ziexp }(${ziexp  - spec.exBias})|${ ziman.toLong.toBinaryString})(${new RealGeneric(spec,  zisgn,  ziexp.toInt,  ziman).toDouble}) != " +
                                f"ref(${ z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble}), " +
                                f"libc(${  zsgn}|${  zexp}(${  zexp - spec.exBias})|${  zman.toLong.toBinaryString})(${z.toDouble}), "
                                )
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
      true, true, true)
  val complexPipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(3),
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



  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 1
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-1+2^-4,  -0.5)",    generateRealWithin(-1.0+pow(2.0, -4), -0.5,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-0.5,  -2^-4)",    generateRealWithin(-0.5, -pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-0, 0)",   generateRealWithin(-0.0, 0.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 2^-4,  0.5)",     generateRealWithin(pow(2.0, -4), 0.5,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 0.5,   1-2^-4)",     generateRealWithin(0.5, 1.0 - pow(2.0, -4),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4), 1.0,_,_))

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


  val float48Spec = new RealSpec(12, 2047, 35)

  val nOrderFP48 = 3
  val adrWFP48 = 10
  val extraBitsFP48 = 4
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, MathFuncPipelineConfig.none,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, MathFuncPipelineConfig.none,
    n, r, "Test Within (-1+2^-4,  -0.5)",    generateRealWithin(-1.0+pow(2.0, -4), -0.5,_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, MathFuncPipelineConfig.none,
    n, r, "Test Within (-0.5,  -2^-4)",    generateRealWithin(-0.5, -pow(2.0, -4),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, MathFuncPipelineConfig.none,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 2^-4,  0.5)",     generateRealWithin(pow(2.0, -4), 0.5,_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, MathFuncPipelineConfig.none,
    n, r, "Test Within ( 0.5,   1-2^-4)",     generateRealWithin(0.5, 1.0 - pow(2.0, -4),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, MathFuncPipelineConfig.none,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))


  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, simplePipeline,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, simplePipeline,
    n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, simplePipeline,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, simplePipeline,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, simplePipeline,
    n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, simplePipeline,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))


  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, complexPipeline,
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, complexPipeline,
    n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, complexPipeline,
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, complexPipeline,
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, complexPipeline,
    n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
  runtest(float48Spec, nOrderFP48, adrWFP48, extraBitsFP48, complexPipeline,
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
}
