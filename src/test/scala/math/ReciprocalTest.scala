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
// Testing Reciprocal using ChiselTest
//

class ReciprocalTest extends FlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test reciprocal"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

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

  private def runtest ( spec : RealSpec, stage : PipelineStageConfig,
    n : Int, r : Random,
    table : FuncTableInt,
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    val reference  = ReciprocalSim.reciprocalSimGeneric( table, _ )
    it should f"1/x pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new ReciprocalGeneric(spec, table.nOrder, table.adrW, table.bp-spec.manW, stage, false, false)) { c =>
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
              assert(zi == z0d, f"x=$xid%x $zi%x!=$z0d%x")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  val reciprocalBF16TableI = ReciprocalSim.reciprocalTableGeneration( 0, 7, 7, 7 )
  val reciprocalF32TableI = ReciprocalSim.reciprocalTableGeneration( 2, 8, 23, 23+2 )

  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(),
    n, r, reciprocalBF16TableI,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(),
    n, r, reciprocalBF16TableI,
    "Test All range",generateRealFull(_,_) )

  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
    n, r, reciprocalF32TableI,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
    n, r, reciprocalF32TableI,
    "Test All range",generateRealFull(_,_) )

}

