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

import rial.math.Log2Sim
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
// Testing Log2 using ChiselTest
//

class MathFuncLog2Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test log2"

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

  private def runtest (spec : RealSpec,
      nOrder : Int, adrW : Int, extraBits : Int, stage: MathFuncPipelineConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (RealSpec, Random) => RealGeneric),
      disableTimeout : Boolean = false
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"log2(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new MathFunctions(spec, nOrder, adrW, extraBits, stage, false, false)) { c =>
        {
          // since table result depends on these values, it is unavoidable to
          // construct tables multiple times.
          val maxCbits   = c.getMaxCbit
          val maxCalcW   = c.getMaxCalcW

          val log2F32TableI = MathFuncLog2Sim.log2NormalTableGeneration(
            RealSpec.Float32Spec, nOrder, adrW, extraBits, Some(maxCalcW), Some(maxCbits))
          val log2F32SmallPositiveTableI = MathFuncLog2Sim.log2SmallPositiveTableGeneration(
            RealSpec.Float32Spec, nOrder, adrW, extraBits, Some(maxCalcW), Some(maxCbits))
          val log2F32SmallNegativeTableI = MathFuncLog2Sim.log2SmallNegativeTableGeneration(
            RealSpec.Float32Spec, nOrder, adrW, extraBits, Some(maxCalcW), Some(maxCbits))

          val reference  = MathFuncLog2Sim.log2SimGeneric(
            log2F32TableI, log2F32SmallPositiveTableI, log2F32SmallNegativeTableI, _ )

          // To avoid timeoutException while testing z == neg.
          // Detailed explanation follows.
          if(disableTimeout) { c.clock.setTimeout(0) }
          // Detailed explanation:
          //     Here, we disable timeout checking in chiseltest. By default,
          // chiseltest fails after 1000 idle cycles. This is a reasonable
          // behavior in most cases. However, unfortunately, the definition of
          // "an idle cycle" is like: a cycle where input or output does not
          // change. Since we are testing log2 function, if we take a negative
          // value, the result becomes nan. So, the test case [-inf, 0] fails
          // after 1000 cycles. The behavior, returning nan, is the correct,
          // expected behavior. But it is indistinguishable from doing nothing
          // without any calculation. To avoid this false positive, we manually
          // disable timeout setting if we test negative values that always
          // return nan. Note that, assertion and `io.expect` still works so
          // there is no problem in the testing functionality.

          val nstage     = c.getStage

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
//             println("---------------------------------------------------------")
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.sel.poke(SelectFunc.Log2)
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            c.io.y.poke(0.U(spec.W.W))
            val zi = c.io.z.peek().litValue.toBigInt
            c.clock.step(1)
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

              assert(zi == z0d, f"x = ${new RealGeneric(spec, xidsgn, xidexp.toInt, xidman).toDouble}(${xidsgn}|${xidexp}(${xidexp - spec.exBias})|${xidman.toLong.toBinaryString}), " +
                                f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) != " +
                                f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})")
              c.io.z.expect(z0d.U(spec.W.W))
            }
          }
        }
      }
    }
  }

  val extraBits = 3

  runtest(RealSpec.Float32Spec, 2, 8, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Large More Than 1 [2, inf]", generateRealWithin(2.0, pow(2.0, 128.0),_,_))
  runtest(RealSpec.Float32Spec, 2, 8, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Small More Than 1 [1+2^-8, 2]",   generateRealWithin(1.0 + pow(2.0, -8), 2.0,_,_))
  runtest(RealSpec.Float32Spec, 2, 8, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Taylor More Than 1 [1, 1+2^-8]",   generateRealWithin(1.0, 1.0 + pow(2.0, -8),_,_))

  runtest(RealSpec.Float32Spec, 2, 8, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Taylor Less Than 1 [1-2^-8, 1]", generateRealWithin(1.0-pow(2.0, -8.0),1.0,_,_))
  runtest(RealSpec.Float32Spec, 2, 8, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Small Less Than 1 [0.5, 1-2^-8]", generateRealWithin(0.5,1.0-pow(2.0, -8.0),_,_))
  runtest(RealSpec.Float32Spec, 2, 8, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Large Less Than 1 [0, 0.5]", generateRealWithin(0.0,0.5-pow(2.0, -24),_,_))

  runtest(RealSpec.Float32Spec, 2, 8, extraBits, MathFuncPipelineConfig.none(), n, r,
    "Test Any Negative [-inf, 0]", generateRealWithin(-pow(2.0, 128), 0.0,_,_),
    /*disableTimeout = */ true)
}

