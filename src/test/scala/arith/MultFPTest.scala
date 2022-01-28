
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
import rial.math._
import rial.arith._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

import spire.math.SafeLong
import spire.implicits._

//
// Testing Pow2F32 using ChiselTest
//

class MultFPTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test x*y"

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

  def multTest(xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig ) = {
      test(new MultFPGeneric( xSpec, ySpec, zSpec, roundSpec, stage )) { c =>
        {
          var q  = new Queue[(BigInt,BigInt,BigInt)]
          val nstage = c.getStage
          for (gen <- List( ("Test Within (-128,128)",generateRealWithin(128.0,_,_)),
                            ("Test All range",generateRealFull(_,_)) ) ) {
            println(gen._1)
            for(i <- 1 to n+nstage) {
              val xr = gen._2(xSpec, r)
              val yr = gen._2(ySpec, r)
              val xi = xr.value.toBigInt
              val yi = yr.value.toBigInt
              val zr = xr.multiply(zSpec, roundSpec, yr)
              val z0i= zr.value.toBigInt
              q += ((xi,yi,z0i))
              c.io.x.poke(xi.U(64.W))
              c.io.y.poke(yi.U(64.W))
              val zi = c.io.z.peek.litValue.toBigInt
              c.clock.step(1)
              if (i > nstage) {
                val (xid,yid,z0d) = q.dequeue()
                assert(zi == z0d, f"x=$xid%16x y=$yid%16x $zi%16x!=$z0d%16x")
              }
            }
            q.clear
          }
        }
      }
  }

  it should f"Multiplier Double with pipereg 0" in {
    multTest( RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }
  it should f"Multiplier Float with pipereg 0" in {
    multTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }
  it should f"Multiplier Double*Float->Double with pipereg 0" in {
    multTest( RealSpec.Float64Spec, RealSpec.Float32Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }
  //runtest(n, PipelineStageConfig.default(2))
}

