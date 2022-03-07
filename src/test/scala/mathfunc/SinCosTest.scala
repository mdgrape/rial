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

import rial.mathfunc._
import rial.arith._
import rial.table._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

class MathFuncSinTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test sin"

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

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  private def runtest (spec : RealSpec,
      nOrder : Int, adrW : Int, extraBits : Int, stage: MathFuncPipelineConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (RealSpec, Random) => RealGeneric)
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"sin(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(spec, nOrder, adrW, extraBits, stage, false, false)) { c =>
        {
          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage

          val sinF32TableI = MathFuncSinSim.sinTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits,
            Some(maxCalcW), Some(maxCbit) )
          val reference  = MathFuncSinSim.sinSimGeneric(sinF32TableI, _ )

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
//             println("------------------")
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(SelectFunc.Sin)
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

              if (zi != z0d) {
                c.io.x.poke(xid.U(64.W))
                for(i <- 1 to nstage) c.clock.step(1)
                c.clock.step(1)
              }

              assert(zi == z0d, f"x =  ${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}(${new RealGeneric(spec, xidsgn, xidexp.toInt, xidman).toDouble}), " +
                                f"test ${zisgn }|${ziexp }(${ziexp  - spec.exBias})|${ziman.toLong.toBinaryString }(${new RealGeneric(spec, zisgn , ziexp .toInt, ziman ).toDouble}) != " +
                                f"ref  ${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString}(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
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

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (-2pi, -pi)",     generateRealWithin(-2.0*Pi,-1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (-pi, 0)",     generateRealWithin(-1.0*Pi,0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (     0, 2^-6)", generateRealWithin(0.0,pow(2.0, -6),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within ( 2^-6, pi/2)",   generateRealWithin(pow(2.0, -6),0.5*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (pi/2, pi)",      generateRealWithin(0.5*Pi,1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (pi, 2pi)",      generateRealWithin(1.0*Pi,2.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (2pi, 10pi)", generateRealWithin(2.0*Pi, 10.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (-10pi, -pi)", generateRealWithin(-10.0*Pi, -1.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (16pi, 32pi)", generateRealWithin(16.0*Pi, 32.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (32pi, 64pi)", generateRealWithin(32.0*Pi, 64.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (-32pi, -16pi)", generateRealWithin(-32.0*Pi, -16.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (-64pi, -32pi)", generateRealWithin(-64.0*Pi, -32.0*Pi,_,_))

  val simplePipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      true, true)

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (-2pi, -pi)",     generateRealWithin(-2.0*Pi,-1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (-pi, 0)",     generateRealWithin(-1.0*Pi,0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (     0, 2^-6)", generateRealWithin(0.0,pow(2.0, -6),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within ( 2^-6, pi/2)",   generateRealWithin(pow(2.0, -6),0.5*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (pi/2, pi)",      generateRealWithin(0.5*Pi,1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (pi, 2pi)",      generateRealWithin(1.0*Pi,2.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (2pi, 10pi)", generateRealWithin(2.0*Pi, 10.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (-10pi, -pi)", generateRealWithin(-10.0*Pi, -1.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (16pi, 32pi)", generateRealWithin(16.0*Pi, 32.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (32pi, 64pi)", generateRealWithin(32.0*Pi, 64.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (-32pi, -16pi)", generateRealWithin(-32.0*Pi, -16.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (-64pi, -32pi)", generateRealWithin(-64.0*Pi, -32.0*Pi,_,_))

  val complexPipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(3),
      PipelineStageConfig.atOut(2),
      true, true)

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (-2pi, -pi)",     generateRealWithin(-2.0*Pi,-1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (-pi, 0)",     generateRealWithin(-1.0*Pi,0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (     0, 2^-6)", generateRealWithin(0.0,pow(2.0, -6),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within ( 2^-6, pi/2)",   generateRealWithin(pow(2.0, -6),0.5*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (pi/2, pi)",      generateRealWithin(0.5*Pi,1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (pi, 2pi)",      generateRealWithin(1.0*Pi,2.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (2pi, 10pi)", generateRealWithin(2.0*Pi, 10.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (-10pi, -pi)", generateRealWithin(-10.0*Pi, -1.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (16pi, 32pi)", generateRealWithin(16.0*Pi, 32.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (32pi, 64pi)", generateRealWithin(32.0*Pi, 64.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (-32pi, -16pi)", generateRealWithin(-32.0*Pi, -16.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (-64pi, -32pi)", generateRealWithin(-64.0*Pi, -32.0*Pi,_,_))
}


class MathFuncCosTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test cos"

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

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  private def runtest (spec : RealSpec,
      nOrder : Int, adrW : Int, extraBits : Int, stage: MathFuncPipelineConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (RealSpec, Random) => RealGeneric)
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"cos(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(spec, nOrder, adrW, extraBits, stage, false, false)) { c =>
        {
          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage

          val sinF32TableI = MathFuncSinSim.sinTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits,
            Some(maxCalcW), Some(maxCbit) )
          val reference  = MathFuncCosSim.cosSimGeneric(sinF32TableI, _ )

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(SelectFunc.Cos)
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

              assert(zi == z0d, f"x =  ${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}(${new RealGeneric(spec, xidsgn, xidexp.toInt, xidman).toDouble}), " +
                                f"test ${zisgn }|${ziexp }(${ziexp  - spec.exBias})|${ziman.toLong.toBinaryString }(${new RealGeneric(spec, zisgn , ziexp .toInt, ziman ).toDouble}) != " +
                                f"ref  ${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString}(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
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

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (-2pi, -pi)",     generateRealWithin(-2.0*Pi,-1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (-pi, 0)",     generateRealWithin(-1.0*Pi,0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (     0, 2^-6)", generateRealWithin(0.0,pow(2.0, -6),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within ( 2^-6, pi/2)",   generateRealWithin(pow(2.0, -6),0.5*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (pi/2, pi)",      generateRealWithin(0.5*Pi,1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Within (pi, 2pi)",      generateRealWithin(1.0*Pi,2.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (2pi, 10pi)", generateRealWithin(2.0*Pi, 10.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (-10pi, -pi)", generateRealWithin(-10.0*Pi, -1.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (16pi, 32pi)", generateRealWithin(16.0*Pi, 32.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (32pi, 64pi)", generateRealWithin(32.0*Pi, 64.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (-32pi, -16pi)", generateRealWithin(-32.0*Pi, -16.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Out of range (-64pi, -32pi)", generateRealWithin(-64.0*Pi, -32.0*Pi,_,_))

  val simplePipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      true, true)

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (-2pi, -pi)",     generateRealWithin(-2.0*Pi,-1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (-pi, 0)",     generateRealWithin(-1.0*Pi,0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (     0, 2^-6)", generateRealWithin(0.0,pow(2.0, -6),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within ( 2^-6, pi/2)",   generateRealWithin(pow(2.0, -6),0.5*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (pi/2, pi)",      generateRealWithin(0.5*Pi,1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Within (pi, 2pi)",      generateRealWithin(1.0*Pi,2.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (2pi, 10pi)", generateRealWithin(2.0*Pi, 10.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (-10pi, -pi)", generateRealWithin(-10.0*Pi, -1.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (16pi, 32pi)", generateRealWithin(16.0*Pi, 32.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (32pi, 64pi)", generateRealWithin(32.0*Pi, 64.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (-32pi, -16pi)", generateRealWithin(-32.0*Pi, -16.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, simplePipeline, n, r,
    "Test Out of range (-64pi, -32pi)", generateRealWithin(-64.0*Pi, -32.0*Pi,_,_))

  val complexPipeline = new MathFuncPipelineConfig(
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(3),
      PipelineStageConfig.atOut(2),
      true, true)

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (-2pi, -pi)",     generateRealWithin(-2.0*Pi,-1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (-pi, 0)",     generateRealWithin(-1.0*Pi,0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (     0, 2^-6)", generateRealWithin(0.0,pow(2.0, -6),_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within ( 2^-6, pi/2)",   generateRealWithin(pow(2.0, -6),0.5*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (pi/2, pi)",      generateRealWithin(0.5*Pi,1.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Within (pi, 2pi)",      generateRealWithin(1.0*Pi,2.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (2pi, 10pi)", generateRealWithin(2.0*Pi, 10.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (-10pi, -pi)", generateRealWithin(-10.0*Pi, -1.0*Pi,_,_))

  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (16pi, 32pi)", generateRealWithin(16.0*Pi, 32.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (32pi, 64pi)", generateRealWithin(32.0*Pi, 64.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (-32pi, -16pi)", generateRealWithin(-32.0*Pi, -16.0*Pi,_,_))
  runtest(RealSpec.Float32Spec, nOrder, adrW, extraBits, complexPipeline, n, r,
    "Test Out of range (-64pi, -32pi)", generateRealWithin(-64.0*Pi, -32.0*Pi,_,_))
}
