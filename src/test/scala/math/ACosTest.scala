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

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  private def runtest ( spec : RealSpec, stage : PipelineStageConfig,
    n : Int, r : Random, table : Seq[FuncTableInt],
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    val reference  = ACosSim.acosSimGeneric( table, _ )
    it should f"acos(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new ACosGeneric(spec, table(0).nOrder, table(0).adrW, table(0).bp-spec.manW, stage, false, false)) { c =>
        {
          val nstage = c.getStage
          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
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

              if (zi != z0d) {
                c.io.x.poke(xid.U(64.W))
                for(i <- 1 to nstage) c.clock.step(1)
                c.clock.step(1)
              }

              val ref = new RealGeneric(spec, acos(xi.toDouble))
              println(f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) ?= ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman}), test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  val acosBF16TableI = ACosSim.acosTableGeneration( 0, 7, 7, 7+2 )

  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(),
    n, r, acosBF16TableI, "Test Within (     0, 2^-13)", generateRealWithin(0.0,pow(2.0, -13),_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(),
    n, r, acosBF16TableI, "Test Within ( 2^-13, 0.5)",     generateRealWithin(pow(2.0, -13),0.5,_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(),
    n, r, acosBF16TableI, "Test Within ( 0.5, 1.0)",     generateRealWithin(0.5,1.0,_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(),
    n, r, acosBF16TableI, "Test Within (-1.0, 0.0)",     generateRealWithin(-1.0,0.0,_,_))

//   val acosF32TableI  = ACosSim.acosTableGeneration( 2, 8, 23, 23+2 )
//
//   runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
//     n, r, acosF32TableI, "Test Within (-0.95,  -2^-8)", generateRealWithin(-0.95, -pow(2.0, -8),_,_))
//   runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
//     n, r, acosF32TableI, "Test Within (-2^-8,  -2^-23)", generateRealWithin(-pow(2.0, -8),-pow(2.0, -23),_,_))
//   runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
//     n, r, acosF32TableI, "Test Within (-2^-23, 0)", generateRealWithin(-pow(2.0, -23),0.0,_,_))
//   runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
//     n, r, acosF32TableI, "Test Within (     0, 2^-23)", generateRealWithin(0.0,pow(2.0, -23),_,_))
//   runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
//     n, r, acosF32TableI, "Test Within ( 2^-23, 2^-8)", generateRealWithin(pow(2.0, -23),pow(2.0, -8),_,_))
//   runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
//     n, r, acosF32TableI, "Test Within ( 2^-8,  0.95)",     generateRealWithin(pow(2.0, -8), 0.95,_,_))
}
