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

import rial.math.SinPiSim
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
// Testing SinPi using ChiselTest
//

class MathFuncSinPiTest extends FlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test sinPi"

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

  private def runtest ( spec : RealSpec,
      nOrder : Int, adrW : Int, extraBits : Int, stage : PipelineStageConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (RealSpec, Random) => RealGeneric)
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"sinPi(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(spec, nOrder, adrW, extraBits, stage, false, false)) { c =>
        {
          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage
          val reference  = SinPiSim.sinPiSimGeneric( SinPiSim.sinPiTableGeneration( nOrder, adrW, spec.manW, spec.manW+extraBits, Some(maxCalcW), Some(maxCbit) ), _, false )

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(SelectFunc.selectSinPi)
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            c.io.y.poke(0.U(spec.W.W))
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

              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r,
    "Test Within (    -1, 0)",     generateRealWithin(-1.0,0.0,_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r,
    "Test Within (     0, 2^-13)", generateRealWithin(0.0,pow(2.0, -13) - pow(2.0, -36),_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r,
    "Test Within ( 2^-13, 0.5)",   generateRealWithin(pow(2.0, -13),0.5,_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r,
    "Test Within (0.5, 1.0)",      generateRealWithin(0.5,1.0,_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r,
    "Test Within (1.0, 2.0)",      generateRealWithin(1.0,2.0,_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r,
    "Test Out of range (2.0, inf)", generateRealWithin(2.0, 1e34,_,_))
  runtest(RealSpec.Float32Spec, 2, 8, 2, PipelineStageConfig.none(), n, r,
    "Test Out of range (-inf, -1)", generateRealWithin(-1e34, -1.0,_,_))
}

