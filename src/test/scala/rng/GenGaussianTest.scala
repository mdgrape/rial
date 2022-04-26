import org.scalatest._

import chisel3._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.arith._
import rial.math._
import rial.rng._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

//
// Test GenRandomFloat12 using chi-squared test.
//

class SinCos2PiTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test uniform real distribution [1, 2)"

  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRandomUInt(width: Int, r: Random) = {
    val rnd = SafeLong(r.nextLong)
    rnd & maskSL(width)
  }

  private def runtest ( isSin: Boolean, rndW: Int, spec : RealSpec, polySpec: PolynomialSpec,
      stage: PipelineStageConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (Int, Random) => SafeLong),
      tolerance : Int
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"sin/cos(2pi*x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new SinCos2Pi(rndW, spec, polySpec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val nstage = stage.total

          val reference = (x: SafeLong) => {
            val xr = x.toDouble / (SafeLong(1) << rndW).toDouble
            val z  = if(isSin) {
              if(x == (SafeLong(0) << (rndW-2)) || x == (SafeLong(2) << (rndW-2))) {0.0} else {sin(2 * Pi * xr)}
            } else {
              if(x == (SafeLong(1) << (rndW-2)) || x == (SafeLong(3) << (rndW-2))) {0.0} else {cos(2 * Pi * xr)}
            }
            new RealGeneric(spec, z)
          }

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(rndW, r)
            val z0r= reference(xi)
            q += ((xi.toBigInt, z0r.value.toBigInt))

            c.io.isSin.poke(isSin.B)
            c.io.x.poke(xi.toBigInt.U(spec.W.W))
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

              val diff = (zi - z0d).abs

              val xr = xid.toDouble / (SafeLong(1) << rndW).toDouble
              val xc = if(xr < 0.25) {xr} else if(xr < 0.5) {0.5 - xr} else if(xr < 0.75) {xr - 0.5} else {1.0 - xr}

              assert(diff <= tolerance,
                     f"x = ${xid.toLong.toBinaryString}(${xr} -> ${xc}), z = ${sin(2 * Pi * xr)}, "+
                     f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble}) != " +
                     f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
            }
          }

          // special values: 0, 0.25, 0.5, 0.75, 1-eps
          for(xid <- Seq(SafeLong(0), SafeLong(1)<<(rndW-2), SafeLong(2)<<(rndW-2), SafeLong(3)<<(rndW-2), maskSL(rndW) )) {
            val z0r = reference(xid)
            val z0d = z0r.value.toBigInt

            c.io.isSin.poke(isSin.B)
            c.io.x.poke(xid.toBigInt.U(spec.W.W))
            for(j <- 0 until nstage) {
              c.clock.step(1)
            }
            val zi = c.io.z.peek().litValue.toBigInt

            val diff = (zi - z0d).abs

            val xr = xid.toDouble / (SafeLong(1) << rndW).toDouble
            val xc = if(xr < 0.25) {xr} else if(xr < 0.5) {0.5 - xr} else if(xr < 0.75) {xr - 0.5} else {1.0 - xr}
            val zr = if(isSin) {sin(2 * Pi * xr)} else {cos(2 * Pi * xr)}

            val xidsgn = bit(spec.W-1, xid).toInt
            val xidexp = slice(spec.manW, spec.exW, xid)
            val xidman = xid & maskSL(spec.manW)

            val zisgn = bit(spec.W-1, zi).toInt
            val ziexp = slice(spec.manW, spec.exW, zi)
            val ziman = zi & maskSL(spec.manW)

            val z0dsgn = bit(spec.W-1, z0d).toInt
            val z0dexp = slice(spec.manW, spec.exW, z0d)
            val z0dman = z0d & maskSL(spec.manW)

            if(diff > tolerance) {
              c.clock.step(1) // run printf
            }

            assert(diff <= tolerance,
                   f"x = ${xid.toLong.toBinaryString}(${xr} -> ${xc}), z(scala.math.FP64) = ${zr}, "+
                   f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble}) != " +
                   f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
          }
        }
      }
    }
  }
  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3
  val polySpecFP32  = new PolynomialSpec(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)

  runtest(true,  32, RealSpec.Float32Spec, polySpecFP32, PipelineStageConfig.none,
    n, r, "Test FP32 sin",generateRandomUInt, 3)
  runtest(false, 32, RealSpec.Float32Spec, polySpecFP32, PipelineStageConfig.none,
    n, r, "Test FP32 cos",generateRandomUInt, 3)
}
