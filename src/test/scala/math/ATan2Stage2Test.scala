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
// Testing ATan2Stage2 using ChiselTest
//

class ATan2Stage2Test extends AnyFlatSpec
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

  var counter = 0
  val specialValues = Seq(
      ( 1.0,  1.0),
      ( 1.0, -1.0),
      (-1.0,  1.0),
      (-1.0, -1.0),

      ( 100.0,  100.0),
      ( 100.0, -100.0),
      (-100.0,  100.0),
      (-100.0, -100.0),

      (Double.PositiveInfinity,  1.0),
      (Double.PositiveInfinity, -1.0),
      (Double.NegativeInfinity,  1.0),
      (Double.NegativeInfinity, -1.0),

      ( 1.0, Double.PositiveInfinity),
      (-1.0, Double.PositiveInfinity),
      ( 1.0, Double.NegativeInfinity),
      (-1.0, Double.NegativeInfinity),

      (Double.PositiveInfinity, Double.PositiveInfinity),
      (Double.PositiveInfinity, Double.NegativeInfinity),
      (Double.NegativeInfinity, Double.PositiveInfinity),
      (Double.NegativeInfinity, Double.NegativeInfinity),


    )
  def generateSpecialValuesX( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx)._1)
  }
  def generateSpecialValuesY( spec: RealSpec, r: Random ) = {
    val idx = counter
    counter += 1
    if(counter >= specialValues.length) {
      counter = 0
    }
    new RealGeneric(spec, specialValues(idx)._2)
  }


  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  private def runtest ( spec : RealSpec,
      nOrder : Int, adrW : Int, extraBits : Int, stage: MathFuncPipelineConfig,
      n : Int, r : Random, generatorStr : String,
      generatorX : ( (RealSpec, Random) => RealGeneric),
      generatorY : ( (RealSpec, Random) => RealGeneric),
      fncfg: MathFuncConfig = MathFuncConfig.all
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"atan2Stage2(x, y) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(fncfg, spec, nOrder, adrW, extraBits, stage, None, false, false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          import FuncKind._

          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage

          val recTable  = ReciprocalSim.reciprocalTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))
          val atanTable = ATan2Stage2Sim.atanTableGeneration(
            nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit))

          val refstage1  = ATan2Stage1Sim.atan2Stage1SimGeneric(recTable, _, _, false)
          val refstage2  = ATan2Stage2Sim.atan2Stage2SimGeneric(atanTable, _, _, _, _)

          val q  = new Queue[(BigInt,BigInt,BigInt)]
          for(i <- 1 to n+2*nstage) {
//             println("test: stage1 -----------------------------------------------------")
            val xi = generatorX(spec,r)
            val yi = generatorY(spec,r)
            val z1r= refstage1(yi, xi) // XXX order is inverted: atan2(y, x)
            val z2r= refstage2(z1r._1, z1r._2, z1r._3, z1r._4)

//             println(f"test: z1r = (${z1r._1.toDouble}, ${z1r._2}, ${z1r._3}, ${z1r._4})")
//             println(f"test: xi = ${xi.toDouble}, yi = ${yi.toDouble}, atan2(yi, xi) = ${atan2(yi.toDouble, xi.toDouble)}, z2r = ${z2r.toDouble}")

            q += ((xi.value.toBigInt, yi.value.toBigInt, z2r.value.toBigInt))

            c.io.sel.poke(fncfg.signal(ATan2Phase1))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            c.io.y.poke(yi.value.toBigInt.U(spec.W.W))

            if(stage.total > 0) {
              c.clock.step(1)
              c.io.sel.poke(fncfg.signalNone())
              c.io.x.poke(0.U(spec.W.W))
              c.io.y.poke(0.U(spec.W.W))
              for(j <- 1 until max(1, stage.total)) {
//                 println(f"out = ${c.io.z.peek().litValue.toBigInt.toLong.toBinaryString}")
                c.clock.step(1)
              }
            }
            val z1i = c.io.z.peek().litValue.toBigInt
//             println("test: stage2 -----------------------------------------------------")


            if(stage.total == 0) {
              c.clock.step(1)
            }

            if(z1i == 0) {
              if(min(xi.toDouble.abs, yi.toDouble.abs) / max(xi.toDouble.abs, yi.toDouble.abs) > pow(2.0, -126)) {
                println(f"z1 is zero. |x|/|y| should be smaller than ${pow(2.0, -126)}. sim out = ${z1r._1.value.toBigInt.toLong.toBinaryString}")
                println(f"x = ${xi.sgn}|${xi.ex}|${xi.man.toLong.toBinaryString}(${xi.toDouble})")
                println(f"y = ${yi.sgn}|${yi.ex}|${yi.man.toLong.toBinaryString}(${yi.toDouble})")
                println(f"min(x,y)/max(x,y) = ${min(xi.toDouble.abs, yi.toDouble.abs) / max(xi.toDouble.abs, yi.toDouble.abs)}")
                ATan2Stage1Sim.atan2Stage1SimGeneric( recTable, yi, xi, true)
              }
              assert(z1r._3 != 0, "if z1 == 0, special value flag should be nonzero")
              if (!xi.isInfinite && !xi.isNaN && !yi.isInfinite && !yi.isNaN) {
                assert(min(xi.toDouble.abs, yi.toDouble.abs) / max(xi.toDouble.abs, yi.toDouble.abs) <= pow(2.0, -126))
              }
            }

            // XXX
            if (!xi.isInfinite && !xi.isNaN && !yi.isInfinite && !yi.isNaN) {
              if(z1r._1.value.toBigInt != z1i) {
                println(f"x/y = ${min(xi.toDouble, yi.toDouble) / max(xi.toDouble, yi.toDouble)}")
                println(f"out = ${new RealGeneric(spec, z1i.toSafeLong).toDouble}")
                println(f"z1r = ${z1r._1.value.toLong.toBinaryString}")
                println(f"z1i = ${z1i.toLong.toBinaryString}")
              }
              assert(z1r._1.value.toBigInt == z1i,
                f"x = ${xi.toDouble}, yi = ${yi.toDouble}, " +
                f"z1r(${z1r._1.value.toLong.toBinaryString}) != z1i(${z1i.toLong.toBinaryString})");
            }
//             println(f"test:z1i = ${z1i.toLong.toBinaryString}")
//             println(f"test:z1.special = ${z1r._3}")

            c.io.sel.poke(fncfg.signal(ATan2Phase2))
            c.io.x.poke(z1i.U(spec.W.W))
            c.io.y.poke(0.U(spec.W.W))

            val z1real = new RealGeneric(spec, z1i.toSafeLong)
//             println(f"z1i = ${z1real.sgn}|${z1real.ex}|${z1real.man.toLong.toBinaryString}(${z1real.toDouble})")

            if(stage.total > 0) {
              c.clock.step(1)
              c.io.sel.poke(fncfg.signalNone())
              c.io.x.poke(0.U(spec.W.W))
              c.io.y.poke(0.U(spec.W.W))
              for(j <- 1 until max(1, stage.total)) {
//                 println(f"out = ${c.io.z.peek().litValue.toBigInt.toLong.toBinaryString}")
                c.clock.step(1)
              }
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
              println(f"test:zcir       = ${z2i.toLong.toBinaryString}(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble})")
            }

            val ref = new RealGeneric(spec, atan2(y.toDouble, x.toDouble))
            assert((ref.value.toLong - z0d.toLong).abs < 4)

            assert(z2i == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
                               f"y = (${yidsgn}|${yidexp}(${yidexp - spec.exBias})|${yidman.toLong.toBinaryString}), " +
                               f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != " +
                               f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
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

//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^24  < y/x <  inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^24  > y/x > -inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^24  < y/x <  inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^24  > y/x > -inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
// 
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^12  < y/x <  2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^12  > y/x > -2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^12  < y/x <  2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^12  > y/x > -2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))
// 
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  1     < y/x <  2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -1     > y/x > -2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  1     < y/x <  2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -1     > y/x > -2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))
// 
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^-12 < y/x <  1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^-12 > y/x > -1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within  2^-12 < y/x <  1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within -2^-12 > y/x > -1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))
// 
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within 0     < y/x <  2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within 0     > y/x > -2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
//   runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, MathFuncPipelineConfig.none, n, r, "Test Special Values", generateSpecialValuesX(_,_), generateSpecialValuesY(_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^24  < y/x <  inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^24  > y/x > -inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^24  < y/x <  inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^24  > y/x > -inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^12  < y/x <  2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^12  > y/x > -2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^12  < y/x <  2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^12  > y/x > -2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  1     < y/x <  2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -1     > y/x > -2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  1     < y/x <  2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -1     > y/x > -2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^-12 < y/x <  1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^-12 > y/x > -1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within  2^-12 < y/x <  1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within -2^-12 > y/x > -1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within 0     < y/x <  2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within 0     > y/x > -2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, simplePipeline, n, r, "Test Special Values", generateSpecialValuesX(_,_), generateSpecialValuesY(_,_))


  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^24  < y/x <  inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^24  > y/x > -inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^24  < y/x <  inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^24  > y/x > -inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^12  < y/x <  2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^12  > y/x > -2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^12  < y/x <  2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^12  > y/x > -2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  1     < y/x <  2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -1     > y/x > -2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  1     < y/x <  2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -1     > y/x > -2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^-12 < y/x <  1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^-12 > y/x > -1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within  2^-12 < y/x <  1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within -2^-12 > y/x > -1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within 0     < y/x <  2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within 0     > y/x > -2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))

  runtest(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, complexPipeline, n, r, "Test Special Values", generateSpecialValuesX(_,_), generateSpecialValuesY(_,_))

  val nOrderBF16 = 0
  val adrWBF16   = 7
  val extraBitsBF16 = 1

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^24  < y/x <  inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^24  > y/x > -inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^24  < y/x <  inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^24  > y/x > -inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^12  < y/x <  2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^12  > y/x > -2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^12  < y/x <  2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^12  > y/x > -2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  1     < y/x <  2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -1     > y/x > -2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  1     < y/x <  2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -1     > y/x > -2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^-12 < y/x <  1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^-12 > y/x > -1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within  2^-12 < y/x <  1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within -2^-12 > y/x > -1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within 0     < y/x <  2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within 0     > y/x > -2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, MathFuncPipelineConfig.none, n, r, "Test Special Values", generateSpecialValuesX(_,_), generateSpecialValuesY(_,_))


  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^24  < y/x <  inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^24  > y/x > -inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^24  < y/x <  inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^24  > y/x > -inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^12  < y/x <  2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^12  > y/x > -2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^12  < y/x <  2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^12  > y/x > -2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  1     < y/x <  2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -1     > y/x > -2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  1     < y/x <  2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -1     > y/x > -2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^-12 < y/x <  1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^-12 > y/x > -1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within  2^-12 < y/x <  1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within -2^-12 > y/x > -1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within 0     < y/x <  2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within 0     > y/x > -2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, simplePipeline, n, r, "Test Special Values", generateSpecialValuesX(_,_), generateSpecialValuesY(_,_))


  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^24  < y/x <  inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^24  > y/x > -inf",               generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^24  < y/x <  inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^24  > y/x > -inf with large x",  generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, 128),-pow(2.0,  24),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^12  < y/x <  2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^12  > y/x > -2^24",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^12  < y/x <  2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^12  > y/x > -2^24 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  24),-pow(2.0,  12),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  1     < y/x <  2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -1     > y/x > -2^12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  1     < y/x <  2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -1     > y/x > -2^12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,  12),-pow(2.0,   0),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^-12 < y/x <  1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^-12 > y/x > -1",                 generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within  2^-12 < y/x <  1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin( pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within -2^-12 > y/x > -1 with large x",    generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0,   0),-pow(2.0, -12),_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within 0     < y/x <  2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within 0     > y/x > -2^-12",              generateRealWithin(-1.0,           1.0,          _,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within 0     < y/x <  2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(           0.0, pow(2.0, -12),_,_))
  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Within 0     > y/x > -2^-12 with large x", generateRealWithin(-pow(2.0, 100), pow(2.0, 100),_,_), generateRealWithin(-pow(2.0, -12),           0.0,_,_))

  runtest(RealSpec.BFloat16Spec, nOrderBF16, adrWBF16, extraBitsBF16, complexPipeline, n, r, "Test Special Values", generateSpecialValuesX(_,_), generateSpecialValuesY(_,_))
}

