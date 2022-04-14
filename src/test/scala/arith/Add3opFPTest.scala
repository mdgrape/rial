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

class Add3opFPTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test x+y+z"

  var n = 1000

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

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def generateRealWithinTriple( p : Double, xSpec: RealSpec, ySpec: RealSpec, zSpec: RealSpec, r : Random ) = {
    val x = generateRealWithin(p, xSpec, r)
    val y = generateRealWithin(p, ySpec, r)
    val z = generateRealWithin(p, zSpec, r)
    (x,y,z)
  }

  def generateRealFullTriple ( xSpec: RealSpec, ySpec: RealSpec, zSpec: RealSpec, r : Random ) = {
    val x = generateRealFull(xSpec, r)
    val y = generateRealFull(ySpec, r)
    val z = generateRealFull(zSpec, r)
    (x,y,z)
  }

  def generateRealZeroNonZeroTriple ( xSpec: RealSpec, ySpec: RealSpec, zSpec: RealSpec, r : Random ) = {
    val x = if(r.nextBoolean) {RealGeneric.zero(xSpec)} else {generateRealFull(xSpec, r)}
    val y = if(r.nextBoolean) {RealGeneric.zero(ySpec)} else {generateRealFull(ySpec, r)}
    val z = if(r.nextBoolean) {RealGeneric.zero(zSpec)} else {generateRealFull(zSpec, r)}
    (x,y,z)
  }

  def add3opTest(xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, wSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig, tolerance: Int ) = {
    test(new ThreeOpAddFPGeneric( xSpec, ySpec, zSpec, wSpec, roundSpec, stage, true ) with DebugControlMaster) { c =>
      {
        c.debugEnableIO.poke(false.B)

        var q  = new Queue[(BigInt,BigInt,BigInt,BigInt)]
        val nstage = c.getStage
        for (gen <- List(
          ("Test Within (-128,128)",generateRealWithinTriple(128.0,_,_,_,_)),
          ("Test All range",generateRealFullTriple(_,_,_,_)),
          ("Test triple includes zero",generateRealZeroNonZeroTriple(_,_,_,_)),
        ) ) {
          println(gen._1)
          for(i <- 1 to n+nstage) {
            val (xr,yr,zr) = gen._2(xSpec, ySpec, zSpec, r)
            val xi = xr.value.toBigInt
            val yi = yr.value.toBigInt
            val zi = zr.value.toBigInt
            val wr = xr.add3op(zSpec, roundSpec, yr, zr)
            val w0i= wr.value.toBigInt
            q += ((xi,yi,zi,w0i))
            c.debugEnableIO.poke(false.B)
            c.io.x.poke(xi.U(xSpec.W.W))
            c.io.y.poke(yi.U(ySpec.W.W))
            c.io.z.poke(zi.U(zSpec.W.W))
            val wi = c.io.w.peek().litValue.toBigInt
            //if (zi != z0d) c.debugControlIO.poke(true.B)
            c.clock.step(1)
            if (i > nstage) {
              val (xid,yid,zid,w0d) = q.dequeue()
              if (wi != w0d) {
                c.io.x.poke(xid.U(xSpec.W.W))
                c.io.y.poke(yid.U(ySpec.W.W))
                c.io.z.poke(zid.U(zSpec.W.W))
                for(i <- 1 to nstage) c.clock.step(1)
                c.debugEnableIO.poke(true.B)
                c.clock.step(1)
              }

              val xsgn = bit(xSpec.W-1, xid)
              val xex  = slice(xSpec.manW, xSpec.exW, xid).toInt
              val xman = slice(0, xSpec.manW, xid)

              val ysgn = bit(ySpec.W-1, yid)
              val yex  = slice(ySpec.manW, ySpec.exW, yid).toInt
              val yman = slice(0, ySpec.manW, yid)

              val zsgn = bit(zSpec.W-1, zi)
              val zex  = slice(zSpec.manW, zSpec.exW, zi).toInt
              val zman = slice(0, zSpec.manW, zi)

              val wsgn = bit(wSpec.W-1, wi)
              val wex  = slice(wSpec.manW, wSpec.exW, wi).toInt
              val wman = slice(0, wSpec.manW, wi)

              val rsgn = bit(wSpec.W-1, w0d)
              val rex  = slice(wSpec.manW, wSpec.exW, w0d).toInt
              val rman = slice(0, wSpec.manW, w0d)

              val diff = (wi - w0d).abs
              if(wi != w0d) {
                println(f"diff = ${wi - w0d}.abs = ${diff}, " +
                        f"x=${xsgn}|${xex-xSpec.exBias}|${xman.toLong.toBinaryString}(${new RealGeneric(xSpec, xsgn, xex, xman).toDouble}), " +
                        f"y=${ysgn}|${yex-ySpec.exBias}|${yman.toLong.toBinaryString}(${new RealGeneric(ySpec, ysgn, yex, yman).toDouble}), " +
                        f"z=${zsgn}|${zex-zSpec.exBias}|${zman.toLong.toBinaryString}(${new RealGeneric(zSpec, zsgn, zex, zman).toDouble}), " +
                        f"test=${wsgn}|${wex-wSpec.exBias}|${wman.toLong.toBinaryString}(${new RealGeneric(wSpec, wsgn, wex, wman).toDouble}) != " +
                        f"wref=${rsgn}|${rex-wSpec.exBias}|${rman.toLong.toBinaryString}(${new RealGeneric(wSpec, rsgn, rex, rman).toDouble})")
              }

              // 3op adder has better accuracy than 2x 2op adder because it uses wider intermediate value.
              // so the result might differ from 2x 2op adder result at the LSBs.
              // to ignore those small differences, we set tolerance.
              assert(diff <= tolerance,
                     f"x=${xsgn}|${xex-xSpec.exBias}|${xman.toLong.toBinaryString}(${new RealGeneric(xSpec, xsgn, xex, xman).toDouble}), " +
                     f"y=${ysgn}|${yex-ySpec.exBias}|${yman.toLong.toBinaryString}(${new RealGeneric(ySpec, ysgn, yex, yman).toDouble}), " +
                     f"z=${zsgn}|${zex-zSpec.exBias}|${zman.toLong.toBinaryString}(${new RealGeneric(zSpec, zsgn, zex, zman).toDouble}), " +
                     f"test=${wsgn}|${wex-wSpec.exBias}|${wman.toLong.toBinaryString}(${new RealGeneric(wSpec, wsgn, wex, wman).toDouble}) != " +
                     f"wref=${rsgn}|${rex-wSpec.exBias}|${rman.toLong.toBinaryString}(${new RealGeneric(wSpec, rsgn, rex, rman).toDouble})")
            }
          }
          q.clear()
        }
      }
    }
  }

  it should f"Add Double with pipereg 0" in {
    add3opTest( RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none, 3)
  }
  it should f"Add Float with pipereg 0" in {
    add3opTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none, 3)
  }
}

