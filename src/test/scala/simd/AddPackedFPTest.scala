
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

class AddPackedFPTest extends AnyFlatSpec
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

  def addTest(len : Int, xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig ) = {
    test(new AddPackedFPGeneric(len, xSpec, ySpec, zSpec, roundSpec, stage, true ) with DebugControlMaster) { c =>
      {
        c.debugEnableIO.poke(false.B)

        var q  = new Queue[(Seq[BigInt], Seq[BigInt], Seq[BigInt])] // x, y, z
        val nstage = c.getStage
        for (gen <- List(
          ("Test Within (-128,128)",generateRealWithinPair(128.0,_,_,_)),
          ("Test All range",generateRealFullPair(_,_,_)),
          ("Test small result",generateRealSmallDifferencePair(128.0,60.0,_,_,_))
        ) ) {
          println(gen._1)
          for(i <- 1 to n+nstage) {

            var xyrs = Seq.tabulate(len)(i => {gen._2(xSpec, ySpec, r)})
            val zrs  = Seq.tabulate(len)(i => {xyrs(i)._1.add(zSpec, roundSpec, xyrs(i)._2)})
            val xis  = Seq.tabulate(len)(i => {xyrs(i)._1.value.toBigInt})
            val yis  = Seq.tabulate(len)(i => {xyrs(i)._2.value.toBigInt})
            val z0is = Seq.tabulate(len)(i => {zrs(i).value.toBigInt})

            q += ((xis, yis, z0is))

            for (i <- 0 until len) {
              c.io.x(i).poke(xis(i).U(64.W))
              c.io.y(i).poke(yis(i).U(64.W))
            }
            val zis = Seq.tabulate(len)(i => {c.io.z(i).peek().litValue.toBigInt})

            //if (zi != z0d) c.debugControlIO.poke(true.B)
            c.clock.step(1)
            if (i > nstage) {
              val (xids,yids,z0ds) = q.dequeue()
              val z_eq  = Seq
                .tabulate(len)(i => {(zis(i), z0ds(i))})
                .foldLeft(true)((b : Boolean, zz : (BigInt, BigInt)) => {b && (zz._1 == zz._2)})

              if ( ! z_eq) {
                for (i <- 0 until len) {
                  c.io.x(i).poke(xids(i).U(64.W))
                  c.io.y(i).poke(yids(i).U(64.W))
                }
                for(i <- 1 to nstage) c.clock.step(1) // step `nstage` clocks
                c.debugEnableIO.poke(true.B)
                c.clock.step(1)
              }
              for (i <- 0 until len) {
                val (xid, yid, zi, z0d) = (xids(i), yids(i), zis(i), z0ds(i))
                assert(zis == z0ds, f"x=$xid%16x y=$yid%16x $zi%16x!=$z0d%16x")
              }
            }
          }
          q.clear()
        }
      }
    }
  }

  it should f"Add 4xDouble with pipereg 0" in {
    addTest(4, RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }
  it should f"Add 3xDouble with pipereg 0" in {
    addTest(3, RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }

}

