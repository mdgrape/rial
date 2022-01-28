import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}
import rial.arith._
import rial.fixed._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

import spire.math.SafeLong
import spire.implicits._

//
// Testing Float -> Fixed using ChiselTest
//

class FloatToFixedTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test float -> fixed"

  var n = 10000

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

  def convTest(xSpec : RealSpec, zSpec : FixedSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig, range : Double ) = {
      test(new FloatToFixedGeneric( xSpec, zSpec, roundSpec, stage )) { c =>
        {
          var q  = new Queue[(BigInt,BigInt)]
          val delta  = pow(2.0, -zSpec.fracW)
          val nstage = c.getStage
          val gen = generateRealWithin(range,_,_)
          println(f"Test Within (-${range}, ${range})")
          for(i <- 1 to n+nstage) {
            val x0 = gen(xSpec, r).toDouble
            val z  = (x0 / delta).toLong
            val x  = new RealGeneric(xSpec, z * delta)

            val xi = x.value.toBigInt
            val ziref = if(z >= 0) { z.toBigInt } else { (z & maskL(zSpec.W)).toBigInt }

            q += ((xi, ziref))
            c.io.x.poke(xi.U(xSpec.W.W))
            val zi = c.io.z.peek().litValue.toBigInt
            c.clock.step(1)
            if (i > nstage) {
              val (xid,zid) = q.dequeue()
              assert(zi == zid, f"x=${xid.toLong.toBinaryString}(${new RealGeneric(xSpec, xid).toDouble}) actual(${zi.toLong.toBinaryString}) != expect(${zid.toLong.toBinaryString})")
            }
          }
          q.clear()
        }
      }
  }

  it should f"Converter signed Q24 to Float with pipereg 0" in {
    convTest( RealSpec.Float32Spec, new FixedSpec(32, 24, true, true),
      RoundSpec.roundToEven, n, PipelineStageConfig.none(), pow(2.0, 32-24-1))
  }
  it should f"Converter signed Q30 to Float with pipereg 0" in {
    convTest( RealSpec.Float32Spec, new FixedSpec(32, 30, true, true),
      RoundSpec.roundToEven, n, PipelineStageConfig.none(), pow(2.0, 32-30-1))
  }
  it should f"Converter signed Q16 to Float with pipereg 0" in {
    convTest( RealSpec.Float32Spec, new FixedSpec(32, 16, true, true),
      RoundSpec.roundToEven, n, PipelineStageConfig.none(), pow(2.0, 32-16-1))
  }
}

