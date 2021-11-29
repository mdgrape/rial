import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._
//import org.scalatest.flatspec.AnyFlatSpec
//import org.scalatest.matchers.should.Matchers
//import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}
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
// Testing InvSqrt using ChiselTest
//

class InvSqrtTest extends FlatSpec
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
  def generateReal1to4( spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()*3.0)+1.0
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

  private def runtest ( spec : RealSpec, stage : PipelineStageConfig,
    n : Int, r : Random, table : FuncTableInt,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    val reference  = InvSqrtSim.invsqrtSimGeneric( table, _ )
    it should f"1/sqrt(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new InvSqrtGeneric(spec, table.nOrder, table.adrW-1, table.bp-spec.manW, stage, false, false)) { c =>
        {
          val nstage = c.getStage
          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            val zi = c.io.z.peek.litValue.toBigInt
            if (i > nstage) {
              val (xid,z0d) = q.dequeue

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

//               println(f"x    = ${new RealGeneric(spec, xid).toDouble},  1/sqrt(x) = ${1.0 / sqrt(new RealGeneric(spec, xid).toDouble)}")
//               println(f"test = ${new RealGeneric(spec, zi ).toDouble} = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}(${ziman}))")
//               println(f"ref  = ${new RealGeneric(spec, z0d).toDouble} = (${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString}(${z0dman}))")

              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman}), test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman}) != ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman})")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  val invsqrtBF16TableI = InvSqrtSim.invsqrtTableGeneration( 0, 7, 7, 7 )
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(), n, r, invsqrtBF16TableI,
    "Test Within (1,4)",generateReal1to4(_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(), n, r, invsqrtBF16TableI,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(), n, r, invsqrtBF16TableI,
    "Test All range",generateRealFull(_,_) )

  val invsqrtF32TableI = InvSqrtSim.invsqrtTableGeneration( 2, 8, 23, 23+2 )
  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(), n, r, invsqrtF32TableI,
    "Test Within (1,4)",generateReal1to4(_,_))
  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(), n, r, invsqrtF32TableI,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(), n, r, invsqrtF32TableI,
    "Test All range",generateRealFull(_,_) )
}
