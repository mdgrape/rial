
import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._



import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
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

class Sum4PackedFPTest extends AnyFlatSpec
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
    val scale = r.nextDouble() * (q+1.0)
    if (r.nextInt(2)==0) {
      val z =
        if (scale>=q) {
          new RealGeneric(xSpec, 0.0)
        } else {
          x.scalbn( -scale.toInt )
        }
      (y.negate().add(ySpec, RoundSpec.roundToEven,z), y)
    } else {
      val z =
        if (scale>=q) {
          new RealGeneric(ySpec, 0.0)
        } else {
          y.scalbn( -scale.toInt )
        }
      (x,x.negate().add(ySpec, RoundSpec.roundToEven,z))
    }
  }

  def addTest(xSpec : RealSpec, mask : Seq[Boolean], zSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig ) = {
    test(new Sum4PackedFPGeneric(xSpec, zSpec, roundSpec, stage, true ) with DebugControlMaster) { c =>
      {
        c.debugEnableIO.poke(false.B)

        var q  = new Queue[(Seq[BigInt], Seq[Boolean], BigInt)] // x, mask, z
        val nstage = c.getStage
        for (gen <- List(
          ("Test Within (-128,128)",generateRealWithinPair(128.0,_,_,_)),
          ("Test All range",generateRealFullPair(_,_,_)),
          ("Test small result",generateRealSmallDifferencePair(128.0,60.0,_,_,_))
        ) ) {
          println(gen._1)
          for(i <- 1 to n+nstage) {

            val xrs = Seq.tabulate(4)(i => {gen._2(xSpec, xSpec, r)._1})
            val xis = Seq.tabulate(4)(i => {xrs(i).value.toBigInt});
            val xmasked = Seq.tabulate(4)(i => {
              if (mask(i)) {
                new RealGeneric(xSpec, 0.0)
              } else {
                xrs(i)
              }
            })
            val tmp1 = xmasked(0).add(xSpec, roundSpec, xmasked(1))
            val tmp2 = xmasked(2).add(xSpec, roundSpec, xmasked(3))
            val zr   = tmp1.add(zSpec, roundSpec, tmp2)
            val zi   = zr.value.toBigInt

            q += ((xis, mask, zi))

            for(i <- 0 until 4) {
              c.io.x(i).poke(xis(i).U(64.W))
              c.io.mask(i).poke(mask(i).B)
            }
            val z0i = c.io.z.peek().litValue.toBigInt

            //if (zi != z0d) c.debugControlIO.poke(true.B)
            c.clock.step(1)
            if (i > nstage) {
              val (xids, masks, zd) = q.dequeue()
              val x0 = xids(0)
              val x1 = xids(1)
              val x2 = xids(2)
              val x3 = xids(3)
              assert(z0i == zd, f"xs=[$x0%16x,$x1%16x,$x2%16x,$x3%16x], $z0i%16x!=$zd%16x")
            }
          }
          q.clear()
        }
      }
    }
  }

  it should f"Add 4xDouble with pipereg 0" in {
    addTest( RealSpec.Float64Spec, Seq(false, false, false, false), RealSpec.Float64Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }

  it should f"Add masked 4xDouble with pipereg 0" in {
    addTest( RealSpec.Float64Spec, Seq(true, false, true, false), RealSpec.Float64Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }

}

