
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

class MultFPIntTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test x:Float * y:Int"

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

  def generateLongWithin( width : Int, signed : Boolean, r : Random ) = {
    if (signed) {
      val sgn = if ( r.nextBoolean() ) { -1 } else { 1 }
      val rI : Long = r.nextLong() & maskL(width-1)
      rI * sgn
    } else {
      val rI : Long = r.nextLong() & maskL(width)
      rI
    }
  }

  def multTest(xSpec : RealSpec, yWidth : Int, ySigned : Boolean, zSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig ) = {
      test(new MultFPIntGeneric( xSpec, yWidth, ySigned, zSpec, roundSpec, stage )) { c =>
        {
          var q  = new Queue[(BigInt,BigInt,BigInt)]
          val nstage = c.getStage
          for (gen <- List( ("Test Within (-128,128)",generateRealWithin(128.0,_,_)),
                            ("Test All range",        generateRealFull(_,_)) ) ) {
            println(gen._1)
            for(i <- 1 to n+nstage) {
              val xr = gen._2(xSpec, r)
              val yr = generateLongWithin(yWidth, ySigned, r)
              val xi = xr.value.toBigInt
              val yi = yr.toBigInt
              val zr = new RealGeneric(zSpec, xr.toDouble * yr) // TODO write simulator
              val z0i= zr.value.toBigInt
              println(f"x = ${xi}(${xr.toDouble}), y = ${yr}, z = ${zr.value}")
              q += ((xi,yi,z0i))
              c.io.x.poke(xi.U(xSpec.W.W))
              if(ySigned) {
                c.io.y.poke((yi & maskL(yWidth)).U(yWidth.W))
              } else {
                c.io.y.poke(yi.U(yWidth.W))
              }
              val zi = c.io.z.peek.litValue.toBigInt
              c.clock.step(1)
              if (i > nstage) {
                val (xid,yid,z0d) = q.dequeue()

                val zisgn  = bit(zSpec.W-1, zi).toInt
                val ziexp  = slice(zSpec.manW, zSpec.exW, zi)
                val ziman  = zi & maskSL(zSpec.manW)
                val z0dsgn = bit(zSpec.W-1, z0d).toInt
                val z0dexp = slice(zSpec.manW, zSpec.exW, z0d)
                val z0dman = z0d & maskSL(zSpec.manW)
                assert(zisgn == z0dsgn && ziexp == z0dexp && ziman == z0dman,
                  f"test(${zisgn}|${ziexp}|${ziman}) != ref(${z0dsgn}|${z0dexp}|${z0dman})")
              }
            }
            q.clear()
          }
        }
      }
  }

  it should f"Multiplier Float * UInt8 with pipereg 0" in {
    multTest( RealSpec.Float32Spec, 8, false, RealSpec.Float32Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }
  it should f"Multiplier Float * SInt8 with pipereg 0" in {
    multTest( RealSpec.Float32Spec, 8, true, RealSpec.Float32Spec,
      RoundSpec.roundToEven, n, PipelineStageConfig.none())
  }
}

