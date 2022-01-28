
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
import rial.simd._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

import spire.math.SafeLong
import spire.implicits._

//
// Testing CrossMultPackedFP using ChiselTest
//

class CrossMultPackedFPTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test x (cross*) y"

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
    stage : PipelineStageConfig ) = {
      test(new CrossMultPackedFPGeneric( xSpec, ySpec, zSpec, roundSpec, stage )) { c =>
        {
          var q  = new Queue[((BigInt,BigInt,BigInt,BigInt),
                              (BigInt,BigInt,BigInt,BigInt),
                              (BigInt,BigInt,BigInt,BigInt))]

          val nstage = c.getStage
          for (gen <- List( ("Test Within (-128,128)",generateRealWithin(128.0,_,_)),
                            ("Test All range",generateRealFull(_,_)) ) ) {
            println(gen._1)
            for(i <- 1 to n+nstage) {
              val x0r = gen._2(xSpec, r)
              val x1r = gen._2(xSpec, r)
              val x2r = gen._2(xSpec, r)
              val x3r = gen._2(xSpec, r)

              val y0r = gen._2(ySpec, r)
              val y1r = gen._2(ySpec, r)
              val y2r = gen._2(ySpec, r)
              val y3r = gen._2(ySpec, r)

              val x0i = x0r.value.toBigInt
              val x1i = x1r.value.toBigInt
              val x2i = x2r.value.toBigInt
              val x3i = x3r.value.toBigInt

              val y0i = y0r.value.toBigInt
              val y1i = y1r.value.toBigInt
              val y2i = y2r.value.toBigInt
              val y3i = y3r.value.toBigInt

              // cross product:
              // z0 = x1 * y2
              // z1 = x2 * y0
              // z2 = x0 * y1
              // z3 = x3 * y3

              val z0r = x1r.multiply(zSpec, roundSpec, y2r)
              val z1r = x2r.multiply(zSpec, roundSpec, y0r)
              val z2r = x0r.multiply(zSpec, roundSpec, y1r)
              val z3r = x3r.multiply(zSpec, roundSpec, y3r)

              val zref0i = z0r.value.toBigInt
              val zref1i = z1r.value.toBigInt
              val zref2i = z2r.value.toBigInt
              val zref3i = z3r.value.toBigInt

              q += (((x0i, x1i, x2i, x3i), (y0i, y1i, y2i, y3i), (zref0i, zref1i, zref2i, zref3i)))

              c.io.x(0).poke(x0i.U(64.W))
              c.io.x(1).poke(x1i.U(64.W))
              c.io.x(2).poke(x2i.U(64.W))
              c.io.x(3).poke(x3i.U(64.W))

              c.io.y(0).poke(y0i.U(64.W))
              c.io.y(1).poke(y1i.U(64.W))
              c.io.y(2).poke(y2i.U(64.W))
              c.io.y(3).poke(y3i.U(64.W))

              val z0i = c.io.z(0).peek.litValue.toBigInt
              val z1i = c.io.z(1).peek.litValue.toBigInt
              val z2i = c.io.z(2).peek.litValue.toBigInt
              val z3i = c.io.z(3).peek.litValue.toBigInt

              c.clock.step(1)
              if (i > nstage) {
                val ((x0id, x1id, x2id, x3id),(y0id,y1id,y2id,y3id), (z0d, z1d, z2d, z3d)) = q.dequeue()
                assert(z0i == z0d, f"x=$x0id%16x y=$y0id%16x $z0i%16x!=$z0d%16x")
                assert(z1i == z1d, f"x=$x1id%16x y=$y1id%16x $z1i%16x!=$z1d%16x")
                assert(z2i == z2d, f"x=$x0id%16x y=$y0id%16x $z0i%16x!=$z0d%16x")
                assert(z3i == z3d, f"x=$x1id%16x y=$y1id%16x $z1i%16x!=$z1d%16x")
              }
            }
            q.clear
          }
        }
      }
  }

  it should f"Multiplier Double with pipereg 0" in {
    multTest( RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, PipelineStageConfig.none())
  }
  it should f"Multiplier Float with pipereg 0" in {
    multTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
      RoundSpec.roundToEven, PipelineStageConfig.none())
  }
  it should f"Multiplier Double*Float->Double with pipereg 0" in {
    multTest( RealSpec.Float64Spec, RealSpec.Float32Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, PipelineStageConfig.none())
  }
//   runtest(n, PipelineStageConfig.default(2))
}

