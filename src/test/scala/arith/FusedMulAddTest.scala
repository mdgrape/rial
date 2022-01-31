
import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._



import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}
import rial.math._
import rial.arith._
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
// Testing FP FMA
//

class FMAFPTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test w=x*y+z"

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

  def generateRealExWithin( from : Double, to : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()*(to - from)) + from
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def generateRealWithinPair( p : Double, xSpec: RealSpec, ySpec: RealSpec, zSpec: RealSpec, r : Random ) = {
    val x = generateRealWithin( p, xSpec, r)
    val y = generateRealWithin( p, ySpec, r)
    val z = generateRealWithin( p, zSpec, r)
    (x,y,z)
  }

  def generateRealFullPair ( xSpec: RealSpec, ySpec: RealSpec, zSpec: RealSpec, r : Random ) = {
    val x = generateRealFull(xSpec, r)
    val y = generateRealFull(ySpec, r)
    val z = generateRealFull(zSpec, r)
    (x,y,z)
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

  def generateClosePathTriple ( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = generateRealExWithin(  64.0, 128.0, xSpec, r )
    val y = generateRealExWithin(   0.5,   2.0, ySpec, r )
    val z = generateRealExWithin(-128.0, -64.0, zSpec, r )
    (x, y, z)
  }

  def generateCloseAddPathTriple ( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = generateRealExWithin(64.0, 128.0, xSpec, r )
    val y = generateRealExWithin( 0.5,   2.0, ySpec, r )
    val z = generateRealExWithin(64.0, 128.0, zSpec, r )
    (x, y, z)
  }

  def generateFarProdPathTriple ( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = generateRealExWithin( 64.0, 128.0, xSpec, r )
    val y = generateRealExWithin( 64.0, 128.0, ySpec, r )
    val z = generateRealExWithin(  1.0,   2.0, zSpec, r )
    (x, y, z)
  }

  def generateFarAddendPathTriple ( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = generateRealExWithin(  1.0,   2.0, xSpec, r )
    val y = generateRealExWithin(  1.0,   2.0, ySpec, r )
    val z = generateRealExWithin( 64.0, 128.0, zSpec, r )
    (x, y, z)
  }

  // TODO: make it efficient
  def generateEdgeCaseFloatTriple1( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = new RealGeneric(xSpec, 0, xSpec.exBias,    0)
    val y = new RealGeneric(ySpec, 0, ySpec.exBias,    maskL(ySpec.manW))
    val z = new RealGeneric(zSpec, 0, zSpec.exBias-1-ySpec.manW, maskL(zSpec.manW))
    (x, y, z)
  }
  def generateEdgeCaseFloatTriple2( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = new RealGeneric(xSpec, 0, xSpec.exBias,    0)
    val y = new RealGeneric(ySpec, 0, ySpec.exBias-1-zSpec.manW, maskL(ySpec.manW))
    val z = new RealGeneric(zSpec, 0, zSpec.exBias, maskL(zSpec.manW))
    (x, y, z)
  }
  def generateEdgeCaseFloatTriple3( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = new RealGeneric(xSpec, 0, xSpec.exBias-1, 0)
    val y = new RealGeneric(ySpec, 0, ySpec.exBias-1, Integer.parseInt(  "01010101010101010101011", 2))
    val z = new RealGeneric(zSpec, 0, zSpec.exBias, Integer.parseInt("11101010101010101010101", 2))
    (x, y, z)
  }

  def fmaTest(xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, wSpec : RealSpec,
    roundSpec : RoundSpec, n : Int, stage : PipelineStageConfig ) = {
    test(new FusedMulAddFPGeneric( xSpec, ySpec, zSpec, wSpec, roundSpec, stage, true ) with DebugControlMaster) { c =>
      {
        c.debugEnableIO.poke(false.B)

        var q  = new Queue[(BigInt,BigInt,BigInt,BigInt)]
        val nstage = c.getStage
        for (gen <- List(
          ("Test close path", generateClosePathTriple(_,_,_,_)),
          ("Test close add path", generateCloseAddPathTriple(_,_,_,_)),
          ("Test far prod path", generateFarProdPathTriple(_,_,_,_)),
          ("Test far addend path", generateFarAddendPathTriple(_,_,_,_)),
          ("Test Within (-128,128)",generateRealWithinPair(128.0,_,_,_,_)),
          ("Test All range",generateRealFullPair(_,_,_,_)),
          ("Test edge case",generateEdgeCaseFloatTriple1(_,_,_,_)),
          ("Test edge case",generateEdgeCaseFloatTriple2(_,_,_,_)),
          ("Test edge case",generateEdgeCaseFloatTriple3(_,_,_,_))
        ) ) {
          println(gen._1)
          for(i <- 1 to n+nstage) {
            val (xr, yr, zr) = gen._2(xSpec, ySpec, zSpec, r)

            val xi = xr.value.toBigInt
            val yi = yr.value.toBigInt
            val zi = zr.value.toBigInt

            val wr = xr.fmadd(wSpec, roundSpec, yr, zr)
            val wi = wr.value.toBigInt

            q += ((xi,yi,zi,wi))

            c.io.x.poke(xi.U(xSpec.W.W))
            c.io.y.poke(yi.U(ySpec.W.W))
            c.io.z.poke(zi.U(zSpec.W.W))

            val wi0 = c.io.w.peek().litValue.toBigInt

            c.clock.step(1)
            if (i > nstage) {
              val (xid,yid,zid,wid) = q.dequeue()
              if (wi0 != wid) {
                c.io.x.poke(xid.U(xSpec.W.W))
                c.io.y.poke(yid.U(ySpec.W.W))
                c.io.z.poke(zid.U(zSpec.W.W))
                for(i <- 1 to nstage) c.clock.step(1)
                c.debugEnableIO.poke(true.B)
                c.clock.step(1)
              }

              val wi0sgn = bit(wSpec.W-1, wi0).toInt
              val wi0exp = slice(wSpec.manW, wSpec.exW, wi0)
              val wi0man = wi0 & maskSL(wSpec.manW)

              val widsgn = bit(wSpec.W-1, wid).toInt
              val widexp = slice(wSpec.manW, wSpec.exW, wid)
              val widman = wid & maskSL(wSpec.manW)
              assert(wi0 == wid, f"test(${wi0sgn}|${wi0exp}|${wi0man}) != ref(${widsgn}|${widexp}|${widman})")
            }
          }
          q.clear()
        }
      }
    }
  }

  it should f"FMA float32 with pipereg 0" in {
    fmaTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }
  it should f"FMA float64 with pipereg 0" in {
    fmaTest( RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }
}

