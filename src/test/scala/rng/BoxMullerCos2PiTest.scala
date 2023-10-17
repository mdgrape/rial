import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import rial.arith._
import rial.math._
import rial.rng._

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

class HTBoxMullerCosTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Box-Muller cos(2pix)"

  val r = new Random(123456789)

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  private def runtest (
    cfg:         BoxMullerConfig,
    n:           Int,
    description: String
  ) = {
    it should f"cos(2pix) ${description}" in {
      test( new HTBoxMullerSinCos2Pi(cfg, isSin=false)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val spec = cfg.realSpec
          val nstage = c.nStage

          val q  = new Queue[(Double, Long)]

          assert(cfg.rndW < 64)

//           val specialValues = Seq(1.0, 1.0-pow(2.0, -spec.manW))
//           var specialIdx = 0

          for(i <- 1 to n+nstage) {
            val xi0 = r.nextDouble
            val xi  = if(xi0 == 0.0) {1.0} else {xi0}
            val x   = new RealGeneric(spec, xi)
            val xr  = x.toFloat

            val pif = (new RealGeneric(spec, Pi)).toFloat
            val z0r = cos(2 * pif * xr)
            val z0  = new RealGeneric(spec, z0r)

            val y0  = new RealGeneric(spec, z0r / xr)
            println("============================================================")
            println(f"x = ${xi}, cos(2pix) = ${z0r}")

            q += ((xi, z0.value.toLong))

            c.io.en    .poke(true.B)
            c.io.x.sgn .poke(0.U)
            c.io.x.ex  .poke(x.ex.toBigInt.U)
            c.io.x.man .poke(x.man.toBigInt.U)
            c.io.x.zero.poke(x.isZero.B)
            c.io.x.inf .poke(false.B)
            c.io.x.nan .poke(false.B)

            val zi = c.io.z.peek().litValue.toLong

            val zisgn = bit(spec.W-1, zi).toInt
            val ziexp = slice(spec.manW, spec.exW, zi)
            val ziman = zi & maskSL(spec.manW)

            c.clock.step(1)

//             println(f"current zi = (${zisgn }|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn.toInt, ziexp.toInt, ziman).toDouble})")

            if (i > nstage) {

              val (xref, zref) = q.dequeue()

              val rel = Seq(xref, abs(xref-0.25), abs(xref-0.5), abs(xref-0.75)).min

              val zrefsgn = bit(spec.W-1, zref).toInt
              val zrefexp = slice(spec.manW, spec.exW, zref)
              val zrefman = zref & maskSL(spec.manW)

              val log_threshold = floor(-log(rel) / log(2.0)) + 1
              val threshold = pow(2.0, log_threshold)
              println(f"x = ${xref}, diff = ${rel}, log(threshold) = ${log_threshold}, threshold = ${threshold}, diff = ${abs(zi-zref)}")

              assert(abs(zi - zref) <= threshold, f"x = ${xref}, "+
                  f"test(${zisgn }|${ziexp  }(${ziexp   - spec.exBias})|${ziman  .toLong.toBinaryString})(${new RealGeneric(spec, zisgn.toInt,   ziexp.toInt,   ziman).toDouble}) != " +
                  f"ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString})(${new RealGeneric(spec, zrefsgn.toInt, zrefexp.toInt, zrefman).toFloat})")
            }
          }
        }
      }
    }
  }

  val rndW          = 32
  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3

  val cfgNone = BoxMullerConfig(
      rndW,
      realSpec = RealSpec.Float32Spec,
      polySpec = new PolynomialSpec(RealSpec.Float32Spec.manW, nOrderFP32, adrWFP32, extraBitsFP32),
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      false, false, false,
      PipelineStageConfig.none
    )

  val cfgSimple = BoxMullerConfig(
      rndW,
      realSpec = RealSpec.Float32Spec,
      polySpec = new PolynomialSpec(RealSpec.Float32Spec.manW, nOrderFP32, adrWFP32, extraBitsFP32),
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      PipelineStageConfig.none,
      true,
      true,
      true,
      PipelineStageConfig.none
    )

  val cfgFull = BoxMullerConfig(
      rndW,
      realSpec = RealSpec.Float32Spec,
      polySpec = new PolynomialSpec(RealSpec.Float32Spec.manW, nOrderFP32, adrWFP32, extraBitsFP32),
      PipelineStageConfig.atOut(1),
      PipelineStageConfig.atOut(3),
      PipelineStageConfig.atOut(1),
      true, true, true,
      PipelineStageConfig.atOut(3)
    )

  runtest(cfgNone,   1000, "Test w/ no registers")
  runtest(cfgSimple, 1000, "Test w/ simple registers")
  runtest(cfgFull,   1000, "Test w/ full registers")
}
