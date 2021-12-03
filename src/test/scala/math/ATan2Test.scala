import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._
import org.scalatest.FlatSpec
import org.scalatest.Matchers
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
// Testing ATan2 using ChiselTest
//
class ATan2Test extends FlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test atan2"

  var n = 2000

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

  private def runtest ( spec : RealSpec, stage : PipelineStageConfig,
    n : Int, r : Random, tableRec : FuncTableInt, tableATan : Seq[FuncTableInt],
    generatorStr : String,
    generatorX      : ( (RealSpec, Random) => RealGeneric),
    generatorYoverX : ( (RealSpec, Random) => RealGeneric)
    ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    val reference  = ATan2Sim.atan2SimGeneric( tableRec, tableATan, _, _ )
    it should f"atan2(y, x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new ATan2Generic(spec, tableRec.nOrder,     tableRec.adrW,     tableRec.bp-spec.manW,
                                   tableATan(0).nOrder, tableATan(0).adrW, tableATan(0).bp-spec.manW,
                             stage, false, false)) { c =>
        {
          val nstage = c.getStage
          val q  = new Queue[(BigInt,BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generatorX(spec,r)
            val yi = generatorYoverX(spec,r)
            val z0r= reference(yi, xi)
            q += ((xi.value.toBigInt, yi.value.toBigInt, z0r.value.toBigInt))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            c.io.y.poke(yi.value.toBigInt.U(spec.W.W))
            val zi = c.io.z.peek.litValue.toBigInt
            if (i > nstage) {
              val (xid,yid,z0d) = q.dequeue

              val zisgn = bit(spec.W-1, zi).toInt
              val ziexp = slice(spec.manW, spec.exW, zi)
              val ziman = zi & maskSL(spec.manW)

              val z0dsgn = bit(spec.W-1, z0d).toInt
              val z0dexp = slice(spec.manW, spec.exW, z0d)
              val z0dman = z0d & maskSL(spec.manW)

              if (zi != z0d) {
                c.io.x.poke(xid.U(64.W))
                c.io.y.poke(yid.U(64.W))
                for(i <- 1 to nstage) c.clock.step(1)
                c.clock.step(1)
              }
              assert(zi == z0d, f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  val atan2BF16ReciprocalTableI = ATan2Sim.atan2BF16ReciprocalTableI
  val atan2BF16ATanTableI       = ATan2Sim.atan2BF16ATanTableI

  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(), n, r, atan2BF16ReciprocalTableI, atan2BF16ATanTableI,
    "Test Within y/x > 2^24", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 127),_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(), n, r, atan2BF16ReciprocalTableI, atan2BF16ATanTableI,
    "Test Within y/x > 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(), n, r, atan2BF16ReciprocalTableI, atan2BF16ATanTableI,
    "Test Within 1 < y/x < 2^12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(), n, r, atan2BF16ReciprocalTableI, atan2BF16ATanTableI,
    "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(), n, r, atan2BF16ReciprocalTableI, atan2BF16ATanTableI,
    "Test Within y/x < 2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_))

  val atan2F32ReciprocalTableI = ATan2Sim.atan2F32ReciprocalTableI
  val atan2F32ATanTableI       = ATan2Sim.atan2F32ATanTableI

  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(), n, r, atan2F32ReciprocalTableI, atan2F32ATanTableI,
    "Test Within y/x > 2^24", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 24), pow(2.0, 127),_,_))
  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(), n, r, atan2F32ReciprocalTableI, atan2F32ATanTableI,
    "Test Within y/x > 2^12",  generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(pow(2.0, 12), pow(2.0, 24),_,_))
  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(), n, r, atan2F32ReciprocalTableI, atan2F32ATanTableI,
    "Test Within 1 < y/x < 2^12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_))
  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(), n, r, atan2F32ReciprocalTableI, atan2F32ATanTableI,
    "Test Within 2^-12 < y/x < 1", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(1.0, pow(2.0, 8),_,_))
  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(), n, r, atan2F32ReciprocalTableI, atan2F32ATanTableI,
    "Test Within y/x < 2^-12", generateRealWithin(-1.0, 1.0,_,_), generateRealWithin(0.0, pow(2.0, -12),_,_))

}
