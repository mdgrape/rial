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
// Testing ATan2Phase1 using ChiselTest
//

class ATan2Phase1Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test atan2Phase1"

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
      generatorX      : ( (RealSpec, Random) => RealGeneric),
      generatorYoverX : ( (RealSpec, Random) => RealGeneric),
      fncfg: MathFuncConfig = MathFuncConfig.all
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"atan2Phase1(x, y) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(fncfg, spec, nOrder, adrW, extraBits, stage, None, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          import FuncKind._

          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage
          val recTable   = ReciprocalSim.reciprocalTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit) )
          val reference  = ATan2Phase1Sim.atan2Phase1SimGeneric(recTable, _, _, false)

          val q  = new Queue[(BigInt,BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
//             println("-----------------------------")
            val xi = generatorX(spec,r)
            val yi = generatorYoverX(spec,r)
            val z0r= reference(xi, yi)._1
            q += ((xi.value.toBigInt, yi.value.toBigInt, z0r.value.toBigInt))
            c.io.sel.poke(fncfg.signal(ATan2Phase1))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            c.io.y.get.poke(yi.value.toBigInt.U(spec.W.W))
            val zi = c.io.z.peek().litValue.toBigInt
            if (i > nstage) {
              val (xid, yid, z0d) = q.dequeue()

              val xidsgn = bit(spec.W-1, xid).toInt
              val xidexp = slice(spec.manW, spec.exW, xid)
              val xidman = xid & maskSL(spec.manW)

              val yidsgn = bit(spec.W-1, yid).toInt
              val yidexp = slice(spec.manW, spec.exW, yid)
              val yidman = yid & maskSL(spec.manW)

              val zisgn = bit(spec.W-1, zi).toInt
              val ziexp = slice(spec.manW, spec.exW, zi)
              val ziman = zi & maskSL(spec.manW)

              val z0dsgn = bit(spec.W-1, z0d).toInt
              val z0dexp = slice(spec.manW, spec.exW, z0d)
              val z0dman = z0d & maskSL(spec.manW)

              val x = new RealGeneric(spec, xidsgn, xidexp.toInt, xidman)
              val y = new RealGeneric(spec, yidsgn, yidexp.toInt, yidman)
              val z = new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman)

//               println(f"x                 = ${xid.toLong.toBinaryString}(${x.toDouble})")
//               println(f"y                 = ${yid.toLong.toBinaryString}(${y.toDouble})")
//               println(f"min(x,y)/max(x,y) = ${if(x.toDouble < y.toDouble) {x.toDouble / y.toDouble} else {y.toDouble / x.toDouble}}")
//               println(f"z                 = ${z0d.toLong.toBinaryString}(${z.toDouble})")

              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
                                f"y = (${yidsgn}|${yidexp}(${yidexp - spec.exBias})|${yidman.toLong.toBinaryString}), " +
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
      true, true, true)

  val complexPipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(3),
      PipelineStageConfig.atOut(2),
      true, true, true)

  val nOrderFP32 = 2
  val adrWFP32   = 8
  val extraBitsFP32 = 3

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^24  < y/x <  inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^24  > y/x > -inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^24  < y/x <  inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^24  > y/x > -inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^12  < y/x <  2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^12  > y/x > -2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^12  < y/x <  2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^12  > y/x > -2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  1     < y/x <  2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -1     > y/x > -2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  1     < y/x <  2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -1     > y/x > -2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^-12 < y/x <  1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^-12 > y/x > -1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^-12 < y/x <  1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^-12 > y/x > -1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within 0     < y/x <  2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within 0     > y/x > -2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))



  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^24  < y/x <  inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^24  > y/x > -inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^24  < y/x <  inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^24  > y/x > -inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^12  < y/x <  2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^12  > y/x > -2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^12  < y/x <  2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^12  > y/x > -2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  1     < y/x <  2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -1     > y/x > -2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  1     < y/x <  2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -1     > y/x > -2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^-12 < y/x <  1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^-12 > y/x > -1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^-12 < y/x <  1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^-12 > y/x > -1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within 0     < y/x <  2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within 0     > y/x > -2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))



  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^24  < y/x <  inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^24  > y/x > -inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^24  < y/x <  inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^24  > y/x > -inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^12  < y/x <  2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^12  > y/x > -2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^12  < y/x <  2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^12  > y/x > -2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  1     < y/x <  2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -1     > y/x > -2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  1     < y/x <  2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -1     > y/x > -2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^-12 < y/x <  1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^-12 > y/x > -1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^-12 < y/x <  1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^-12 > y/x > -1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within 0     < y/x <  2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within 0     > y/x > -2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))

  val nOrderBF16 = 0
  val adrWBF16   = 7
  val extraBitsBF16 = 1


  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^24  < y/x <  inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^24  > y/x > -inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^24  < y/x <  inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^24  > y/x > -inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^12  < y/x <  2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^12  > y/x > -2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^12  < y/x <  2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^12  > y/x > -2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  1     < y/x <  2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -1     > y/x > -2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  1     < y/x <  2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -1     > y/x > -2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^-12 < y/x <  1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^-12 > y/x > -1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^-12 < y/x <  1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^-12 > y/x > -1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within 0     < y/x <  2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within 0     > y/x > -2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))



  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^24  < y/x <  inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^24  > y/x > -inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^24  < y/x <  inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^24  > y/x > -inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^12  < y/x <  2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^12  > y/x > -2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^12  < y/x <  2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^12  > y/x > -2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  1     < y/x <  2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -1     > y/x > -2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  1     < y/x <  2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -1     > y/x > -2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^-12 < y/x <  1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^-12 > y/x > -1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^-12 < y/x <  1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^-12 > y/x > -1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within 0     < y/x <  2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within 0     > y/x > -2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))



  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^24  < y/x <  inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^24  > y/x > -inf", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^24  < y/x <  inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^24  > y/x > -inf with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^12  < y/x <  2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^12  > y/x > -2^24",generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^12  < y/x <  2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^12  > y/x > -2^24 with large x",generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12), _,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  1     < y/x <  2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -1     > y/x > -2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  1     < y/x <  2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -1     > y/x > -2^12 with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0), _,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^-12 < y/x <  1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^-12 > y/x > -1",     generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^-12 < y/x <  1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^-12 > y/x > -1 with large x",     generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12), _,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within 0     < y/x <  2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within 0     > y/x > -2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(         0.0 , pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),          0.0,_,_))


}

