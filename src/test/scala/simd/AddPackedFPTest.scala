
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
import rial.math._
import rial.arith._
import rial.simd._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig
import rial.util.DebugControlMaster

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

import spire.math.SafeLong
import spire.implicits._

//
// Testing FP Adder
//

class AddPackedFPTest extends FlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test xs+ys"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19930514)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def generateRealWithinPair( p : Double, xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    val x = generateRealWithin( p, xSpec, r)
    val y = generateRealWithin( p, ySpec, r)
    (x,y)
  }

  def generateRealFullPair ( xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    val x = generateRealFull(xSpec, r)
    val y = generateRealFull(ySpec, r)
    (x,y)
  }

  // p : x range
  // q : result scale min
  def generateRealSmallDifferencePair ( p : Double, q: Double, xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    val x = generateRealWithin( p, xSpec, r )
    val y = generateRealWithin( p, ySpec, r )
    val scale = r.nextDouble * (q+1.0)
    if (r.nextInt(2)==0) {
      val z =
        if (scale>=q) {
          new RealGeneric(xSpec, 0.0)
        } else {
          x.scalbn( -scale.toInt )
        }
      (y.negate.add(ySpec, RoundSpec.roundToEven,z), y)
    } else {
      val z =
        if (scale>=q) {
          new RealGeneric(ySpec, 0.0)
        } else {
          y.scalbn( -scale.toInt )
        }
      (x,x.negate.add(ySpec, RoundSpec.roundToEven,z))
    }
  }

  def addTest(xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig ) = {
    test(new AddPackedFPGeneric(2, xSpec, ySpec, zSpec, roundSpec, stage, true ) with DebugControlMaster) { c =>
      {
        c.debugEnableIO.poke(false.B)

        var q  = new Queue[((BigInt, BigInt),(BigInt, BigInt),(BigInt, BigInt))] // x, y, z
        val nstage = c.getStage
        for (gen <- List(
          ("Test Within (-128,128)",generateRealWithinPair(128.0,_,_,_)),
          ("Test All range",generateRealFullPair(_,_,_)),
          ("Test small result",generateRealSmallDifferencePair(128.0,60.0,_,_,_))
        ) ) {
          println(gen._1)
          for(i <- 1 to n+nstage) {

            val (xr0,yr0) = gen._2(xSpec, ySpec, r)
            val (xr1,yr1) = gen._2(xSpec, ySpec, r)
            val zr0 = xr0.add(zSpec, roundSpec, yr0)
            val zr1 = xr1.add(zSpec, roundSpec, yr1)

            val xr0i = xr0.value.toBigInt
            val yr0i = yr0.value.toBigInt
            val zr0i = zr0.value.toBigInt

            val xr1i = xr1.value.toBigInt
            val yr1i = yr1.value.toBigInt
            val zr1i = zr1.value.toBigInt

            q += (((xr0i,xr1i), (yr0i, yr1i), (zr0i, zr1i)))

            c.io.x(0).poke(xr0i.U(64.W))
            c.io.x(1).poke(xr1i.U(64.W))
            c.io.y(0).poke(yr0i.U(64.W))
            c.io.y(1).poke(yr1i.U(64.W))
            val z0i = c.io.z(0).peek.litValue.toBigInt
            val z1i = c.io.z(1).peek.litValue.toBigInt

            //if (zi != z0d) c.debugControlIO.poke(true.B)
            c.clock.step(1)
            if (i > nstage) {
              val ((x0id, x1id),(y0id, y1id),(z0d, z1d)) = q.dequeue
              if (z0i != z0d || z1i != z1d) {
                c.io.x(0).poke(x0id.U(64.W))
                c.io.x(1).poke(x1id.U(64.W))
                c.io.y(0).poke(y0id.U(64.W))
                c.io.y(1).poke(y1id.U(64.W))
                for(i <- 1 to nstage) c.clock.step(1) // step `nstage` clocks
                c.debugEnableIO.poke(true.B)
                c.clock.step(1)
              }
              assert(z0i == z0d, f"x=$x0id%16x y=$y0id%16x $z0i%16x!=$z0d%16x")
              assert(z1i == z1d, f"x=$x1id%16x y=$y1id%16x $z1i%16x!=$z1d%16x")
            }
          }
          q.clear
        }
      }
    }
  }

  it should f"Add Double with pipereg 0" in {
    addTest( RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }

}

