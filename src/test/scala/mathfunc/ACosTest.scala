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

//
// Testing ACos using ChiselTest
//

class MathFuncACosTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test acos"

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
      ts : Seq[FuncTableInt], tSqrt : FuncTableInt,
      nOrder : Int, adrW : Int, extraBits : Int, stage : PipelineStageConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (RealSpec, Random) => RealGeneric)
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"acos(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(spec, nOrder, adrW, extraBits, stage, false, false)) { c =>
        {
          val maxCbit    = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW
          val nstage     = c.getStage
          val reference  = MathFuncACosSim.acosSimGeneric( ts, tSqrt, _ )
          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(SelectFunc.ACos)
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

              assert(zi == z0d, f"x = (${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  val acosF32TableI = MathFuncACosSim.acosTableGeneration(
    2, 8, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+2,
    Some(Seq(27, 24, 24)), Some(Seq(27, 24, 24)))

  val sqrtF32TableI = MathFuncACosSim.sqrtTableGeneration(
    2, 8, RealSpec.Float32Spec.manW, RealSpec.Float32Spec.manW+2,
    Some(Seq(27, 24, 24)), Some(Seq(27, 24, 24)))

  runtest(RealSpec.Float32Spec, acosF32TableI, sqrtF32TableI, 2, 8, 2, PipelineStageConfig.none(),
    n, r, "Test Within (-1, -1+2^-4)",     generateRealWithin(-1.0, -1.0+pow(2.0, -4)-pow(2.0, -23),_,_))
  runtest(RealSpec.Float32Spec, acosF32TableI, sqrtF32TableI, 2, 8, 2, PipelineStageConfig.none(),
    n, r, "Test Within (-1+2^-4,  -2^-4)",    generateRealWithin(-1.0+pow(2.0, -4), -pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, acosF32TableI, sqrtF32TableI, 2, 8, 2, PipelineStageConfig.none(),
    n, r, "Test Within (-2^-4, 0)",   generateRealWithin(-pow(2.0, -4), 0,_,_))
  runtest(RealSpec.Float32Spec, acosF32TableI, sqrtF32TableI, 2, 8, 2, PipelineStageConfig.none(),
    n, r, "Test Within ( 0, 2^-4)",     generateRealWithin(0,pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, acosF32TableI, sqrtF32TableI, 2, 8, 2, PipelineStageConfig.none(),
    n, r, "Test Within ( 2^-4,  1-2^-4)",     generateRealWithin(pow(2.0, -4), 1.0 - pow(2.0, -4),_,_))
  runtest(RealSpec.Float32Spec, acosF32TableI, sqrtF32TableI, 2, 8, 2, PipelineStageConfig.none(),
    n, r, "Test Within (1-2^-4, 1)",     generateRealWithin(1.0-pow(2.0, -4)+pow(2.0, -23), 1.0,_,_))
}
