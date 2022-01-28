import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._
//import org.scalatest.flatspec.AnyFlatSpec
//import org.scalatest.matchers.should.Matchers
//import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math.ReciprocalSim
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
// Testing ATan2Stage1 using ChiselTest
//

class MathFuncATan2Stage1Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test atan2Stage1"

  var n = 10000

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
      nOrder : Int, adrW : Int, extraBits : Int, stage : PipelineStageConfig,
      n : Int, r : Random, generatorStr : String,
      generatorYoverX : ( (RealSpec, Random) => RealGeneric)
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"atan2Stage1(x, y) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(spec, nOrder, adrW, extraBits, stage, false, false)) { c =>
        {
          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage
          val reference  = ATan2Stage1Sim.atan2Stage1SimGeneric(
            ReciprocalSim.reciprocalTableGeneration( nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit) ),
            _, _)

          val generatorX = generateRealWithin(-1.0, 1.0, _, _)

          val q  = new Queue[(BigInt,BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
//             println("-----------------------------")
            val xi = generatorX(spec,r)
            val yi = generatorYoverX(spec,r)
            val z0r= reference(xi, yi)._1
            q += ((xi.value.toBigInt, yi.value.toBigInt, z0r.value.toBigInt))
            c.io.sel.poke(SelectFunc.ATan2Stage1)
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            c.io.y.poke(yi.value.toBigInt.U(spec.W.W))
            val zi = c.io.z.peek.litValue.toBigInt
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

              if (zi != z0d) {
                c.io.x.poke(xid.U(64.W))
                for(i <- 1 to nstage) c.clock.step(1)
                c.clock.step(1)
              }

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

//   runtest(RealSpec.BFloat16Spec, 0, 7, 0, PipelineStageConfig.none(),
//     n, r, "Test Within (-128,128)",generateRealWithin(128.0,_,_))
//   runtest(RealSpec.BFloat16Spec, 0, 7, 0, PipelineStageConfig.none(),
//     n, r, "Test All range",generateRealFull(_,_) )

  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r, "Test Within 2^24  < y/x < inf",   generateRealWithin(pow(2.0,  24), pow(2.0, 128),_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r, "Test Within 2^12  < y/x < 2^24",  generateRealWithin(pow(2.0,  12), pow(2.0,  24),_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r, "Test Within 1     < y/x < 2^12",  generateRealWithin(pow(2.0,   0), pow(2.0,  12),_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r, "Test Within 2^-12 < y/x < 1",     generateRealWithin(pow(2.0, -12), pow(2.0,   0),_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r, "Test Within 0     < y/x < 2^-12", generateRealWithin(         0.0 , pow(2.0, -12),_,_))
}

