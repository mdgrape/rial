import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._



import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math.ReciprocalSim
import rial.math.ATan2Sim
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
// Testing ATan2Stage2 using ChiselTest
//

class MathFuncATan2Stage2Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test atan2Stage2"

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
      generatorYoverX : ( (RealSpec, Random) => RealGeneric)
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"atan2Stage2(x, y) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(spec, nOrder, adrW, extraBits, stage, false, false)) { c =>
        {
          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage

          val recTable  = ReciprocalSim.reciprocalTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))
          val atanTable = ATan2Sim.atanTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))

          val refstage1  = ATan2Stage1Sim.atan2Stage1SimGeneric(recTable, _, _)
          val refstage2  = ATan2Stage2Sim.atan2Stage2SimGeneric(atanTable, _, _, _, _)

          val generatorX = generateRealWithin(-1.0, 1.0, _, _)

          val q  = new Queue[(BigInt,BigInt,BigInt)]
          for(i <- 1 to n+2*nstage) {
//             println("test: -----------------------------------------------------")
            val xi = generatorX(spec,r)
            val yi = generatorYoverX(spec,r)
            val z1r= refstage1(yi, xi) // XXX order is inverted: atan2(y, x)
            val z2r= refstage2(z1r._1, z1r._2, z1r._3, z1r._4)

//             println(f"test: z1r = (${z1r._1.toDouble}, ${z1r._2}, ${z1r._3}, ${z1r._4})")
//             println(f"test: xi = ${xi.toDouble}, yi = ${yi.toDouble}, atan2(yi, xi) = ${atan2(yi.toDouble, xi.toDouble)}, z2r = ${z2r.toDouble}")

            q += ((xi.value.toBigInt, yi.value.toBigInt, z2r.value.toBigInt))

            c.io.sel.poke(SelectFunc.ATan2Stage1)
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            c.io.y.poke(yi.value.toBigInt.U(spec.W.W))

            for(j <- 0 until max(1, stage.total)) {
              c.clock.step(1)
            }
            val z1i = c.io.z.peek().litValue.toBigInt

            assert(z1r._1.value.toBigInt == z1i, f"z1r(${z1r._1.value.toLong.toBinaryString}) != z1i(${z1i.toLong.toBinaryString})");
//             println(f"test:z1i = ${z1i.toLong.toBinaryString}")

            c.io.sel.poke(SelectFunc.ATan2Stage2)
            c.io.x.poke(z1i.U(spec.W.W))
            c.io.y.poke(0.U(spec.W.W))

            for(j <- 0 until max(1, stage.total)) {
              c.clock.step(1)
            }
            val z2i = c.io.z.peek().litValue.toBigInt

            val (xid, yid, z0d) = q.dequeue()

            val xidsgn = bit(spec.W-1, xid).toInt
            val xidexp = slice(spec.manW, spec.exW, xid)
            val xidman = xid & maskSL(spec.manW)

            val yidsgn = bit(spec.W-1, yid).toInt
            val yidexp = slice(spec.manW, spec.exW, yid)
            val yidman = yid & maskSL(spec.manW)

            val zisgn = bit(spec.W-1, z2i).toInt
            val ziexp = slice(spec.manW, spec.exW, z2i)
            val ziman = z2i & maskSL(spec.manW)

            val z0dsgn = bit(spec.W-1, z0d).toInt
            val z0dexp = slice(spec.manW, spec.exW, z0d)
            val z0dman = z0d & maskSL(spec.manW)

            val x = new RealGeneric(spec, xidsgn, xidexp.toInt, xidman)
            val y = new RealGeneric(spec, yidsgn, yidexp.toInt, yidman)
            val z = new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman)

            if(z2i != z0d) {
              println(f"test:x          = ${xid.toLong.toBinaryString}(${x.toDouble})")
              println(f"test:y          = ${yid.toLong.toBinaryString}(${y.toDouble})")
              println(f"test:atan2(x,y) = ${atan2(y.toDouble, x.toDouble)}")
              println(f"test:zref       = ${z0d.toLong.toBinaryString}(${z.toDouble})")
              println(f"test:zsim       = ${z2i.toLong.toBinaryString}(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble})")
            }

            assert(z2i == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
                               f"y = (${yidsgn}|${yidexp}(${yidexp - spec.exBias})|${yidman.toLong.toBinaryString}), " +
                               f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != " +
                               f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
          }
        }
      }
    }
  }

  val nOrder = 2
  val adrW   = 8
  val extraBits = 3

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none, n, r, "Test Within 2^24  < y/x < inf",   generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none, n, r, "Test Within 2^12  < y/x < 2^24",  generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none, n, r, "Test Within 1     < y/x < 2^12",  generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none, n, r, "Test Within 2^-12 < y/x < 1",     generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none, n, r, "Test Within 0     < y/x < 2^-12", generateRealWithin(         0.0 , pow(2.0, -12),_,_))

  val simplePipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      true, true)

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r, "Test Within 2^24  < y/x < inf",   generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r, "Test Within 2^12  < y/x < 2^24",  generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r, "Test Within 1     < y/x < 2^12",  generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r, "Test Within 2^-12 < y/x < 1",     generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r, "Test Within 0     < y/x < 2^-12", generateRealWithin(         0.0 , pow(2.0, -12),_,_))

  val complexPipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(3),
      PipelineStageConfig.atOut(2),
      true, true)

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r, "Test Within 2^24  < y/x < inf",   generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r, "Test Within 2^12  < y/x < 2^24",  generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r, "Test Within 1     < y/x < 2^12",  generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r, "Test Within 2^-12 < y/x < 1",     generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r, "Test Within 0     < y/x < 2^-12", generateRealWithin(         0.0 , pow(2.0, -12),_,_))
}

