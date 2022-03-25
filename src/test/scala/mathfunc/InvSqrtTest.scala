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

import rial.math.InvSqrtSim
import rial.mathfunc._
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

class MathFuncInvSqrtTest extends AnyFlatSpec
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
    it should f"invsqrt(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(spec, nOrder, adrW, extraBits, stage, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
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
            c.io.sel.poke(SelectFunc.InvSqrt)
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

  val nOrder = 2
  val adrW = 8
  val extraBits = 3
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none,
    n, r, "Test All range",generateRealFull(_,_) )

  val simplePipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      true, true)

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline,
    n, r, "Test All range",generateRealFull(_,_) )

  val complexPipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(3),
      PipelineStageConfig.atOut(2),
      true, true)

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline,
    n, r, "Test All range",generateRealFull(_,_) )
}


class InvSqrtOnlyTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test invsqrt-only module"

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
    it should f"invsqrt(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new InvSqrtGeneric(spec, nOrder, adrW, extraBits, stage, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val maxCbit    = c.getCbit
          val maxCalcW   = c.getCalcW
          val nstage     = c.getStage
          val invsqrtTable = InvSqrtSim.invsqrtTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit) )
          val reference  = InvSqrtSim.invsqrtSimGeneric( invsqrtTable, _ )

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
              val xidd   = new RealGeneric(spec, xidsgn, xidexp.toInt, xidman).toDouble

              val zisgn = bit(spec.W-1, zi).toInt
              val ziexp = slice(spec.manW, spec.exW, zi)
              val ziman = zi & maskSL(spec.manW)
              val zidd   = new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble

              val z0dsgn = bit(spec.W-1, z0d).toInt
              val z0dexp = slice(spec.manW, spec.exW, z0d)
              val z0dman = z0d & maskSL(spec.manW)
              val z0dd   = new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble

              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString})(1/sqrt(${xidd}) = ${1.0 / sqrt(xidd)}), " +
                                f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${zidd}) != " +
                                f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${z0dd})")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  val nOrderFP32 = 2
  val adrWFP32 = 8
  val extraBitsFP32 = 3
  val complexPipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(3),
      PipelineStageConfig.atOut(2),
      true, true)
  val simplePipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      true, true)

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none,
    n, r, "Test All range",generateRealFull(_,_) )

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline,
    n, r, "Test All range",generateRealFull(_,_) )

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline,
    n, r, "Test All range",generateRealFull(_,_) )

  val nOrderBF16 = 0
  val adrWBF16 = 7
  val extraBitsBF16 = 0

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test Within (-128,128) BF16, ex=0",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none,
    n, r, "Test All range BF16, ex=0",generateRealFull(_,_) )

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test Within (-128,128) BF16, ex=0",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline,
    n, r, "Test All range BF16, ex=0",generateRealFull(_,_) )

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test Within (-128,128) BF16, ex=0",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline,
    n, r, "Test All range BF16, ex=0",generateRealFull(_,_) )

  val nOrderBF16ex1 = 0
  val adrWBF16ex1 = 7
  val extraBitsBF16ex1 = 1

  runtest(RealSpec.BFloat16Spec, nOrderBF16ex1, adrWBF16ex1, extraBitsBF16ex1, MathFuncPipelineConfig.none,
    n, r, "Test Within (-128,128) BF16, ex=1",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16ex1, adrWBF16ex1, extraBitsBF16ex1, MathFuncPipelineConfig.none,
    n, r, "Test All range BF16, ex=1",generateRealFull(_,_) )

  runtest(RealSpec.BFloat16Spec, nOrderBF16ex1, adrWBF16ex1, extraBitsBF16ex1, simplePipeline,
    n, r, "Test Within (-128,128) BF16, ex=1",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16ex1, adrWBF16ex1, extraBitsBF16ex1, simplePipeline,
    n, r, "Test All range BF16, ex=1",generateRealFull(_,_) )

  runtest(RealSpec.BFloat16Spec, nOrderBF16ex1, adrWBF16ex1, extraBitsBF16ex1, complexPipeline,
    n, r, "Test Within (-128,128) BF16, ex=1",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16ex1, adrWBF16ex1, extraBitsBF16ex1, complexPipeline,
    n, r, "Test All range BF16, ex=1",generateRealFull(_,_) )

}
