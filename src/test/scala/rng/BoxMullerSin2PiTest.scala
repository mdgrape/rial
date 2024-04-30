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

class HTBoxMullerSinTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Box-Muller sin(x)"

  val r = new Random(123456789)

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  private def runtest (
    cfg:         HTBoxMullerConfig,
    n:           Int,
    description: String
  ) = {
    it should f"sin(2pix) ${description}" in {
      test( new HTBoxMullerSinCos2Pi(cfg, isSin=true)).
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

            val xr = if(x.toDouble < 0.25) {
              if(x.toDouble == 0.0) {
                0.0
              } else {
                x.toDouble
              }
            } else if (x.toDouble < 0.5) {
              if(x.toDouble == 0.25) {
                0.25
              } else {
                new RealGeneric(spec, 0.5 - x.toDouble).toDouble
              }
            } else if (x.toDouble < 0.75) {
              new RealGeneric(spec, x.toDouble - 0.5).toDouble * (-1.0)
            } else {
              new RealGeneric(spec, 1.0 - x.toDouble).toDouble * (-1.0)
            }

            val z0r = sin(2 * Pi * xr)
            val z0  = new RealGeneric(spec, z0r)

            println("============================================================")
            println(f"x = ${xi}(-> ${x.toDouble}), sin(2pix) = ${z0r}")

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

              val zrefsgn = bit(spec.W-1, zref).toInt
              val zrefexp = slice(spec.manW, spec.exW, zref)
              val zrefman = zref & maskSL(spec.manW)

              assert(abs(zi - zref) <= 7, f"x = ${new RealGeneric(spec, xref).toDouble}, "+
                  f"test(${zisgn }|${ziexp  }(${ziexp   - spec.exBias})|${ziman  .toLong.toBinaryString})(${new RealGeneric(spec, zisgn.toInt,   ziexp.toInt,   ziman  ).toDouble}) != " +
                  f"ref(${zrefsgn}|${zrefexp}(${zrefexp - spec.exBias})|${zrefman.toLong.toBinaryString})(${new RealGeneric(spec, zrefsgn.toInt, zrefexp.toInt, zrefman).toFloat })")
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

  val cfgNone = HTBoxMullerConfig(
      rndW,
      realSpec = RealSpec.Float32Spec,
      polySpec = new PolynomialSpec(RealSpec.Float32Spec.manW, nOrderFP32, adrWFP32, extraBitsFP32),
      int2floatStge = PipelineStageConfig.none,
      polyPreStage  = PipelineStageConfig.none,
      polyCalcStage = PipelineStageConfig.none,
      polyPostStage = PipelineStageConfig.none,
      i2fPolyGap    = false,
      preCalcGap    = false,
      tableCalcGap  = false,
      calcPostGap   = false,
      mulStage = PipelineStageConfig.none
    )

  val cfgSimple = HTBoxMullerConfig(
      rndW,
      realSpec = RealSpec.Float32Spec,
      polySpec = new PolynomialSpec(RealSpec.Float32Spec.manW, nOrderFP32, adrWFP32, extraBitsFP32),
      int2floatStge = PipelineStageConfig.none,
      polyPreStage  = PipelineStageConfig.none,
      polyCalcStage = PipelineStageConfig.none,
      polyPostStage = PipelineStageConfig.none,
      i2fPolyGap    = true,
      preCalcGap    = true,
      tableCalcGap  = true,
      calcPostGap   = true,
      mulStage = PipelineStageConfig.none
    )

  val cfgFull = HTBoxMullerConfig(
      rndW,
      realSpec = RealSpec.Float32Spec,
      polySpec = new PolynomialSpec(RealSpec.Float32Spec.manW, nOrderFP32, adrWFP32, extraBitsFP32),
      int2floatStge = PipelineStageConfig.atOut(1),
      polyPreStage  = PipelineStageConfig.atOut(1),
      polyCalcStage = PipelineStageConfig.atOut(3),
      polyPostStage = PipelineStageConfig.atOut(1),
      i2fPolyGap    = true,
      preCalcGap    = true,
      tableCalcGap  = true,
      calcPostGap   = true,
      mulStage = PipelineStageConfig.atOut(3)
    )

  runtest(cfgNone,   1000, "Test w/ no registers")
  runtest(cfgSimple, 1000, "Test w/ simple registers")
  runtest(cfgFull,   1000, "Test w/ full registers")

  val nOrderBF16    = 0
  val adrWBF16      = 7
  val extraBitsBF16 = 1

  val cfgBF16 = HTBoxMullerConfig(
      rndW,
      realSpec = RealSpec.BFloat16Spec,
      polySpec = new PolynomialSpec(RealSpec.BFloat16Spec.manW, nOrderBF16, adrWBF16, extraBitsBF16),
      int2floatStge = PipelineStageConfig.atOut(0),
      polyPreStage  = PipelineStageConfig.atOut(0),
      polyCalcStage = PipelineStageConfig.atOut(0),
      polyPostStage = PipelineStageConfig.atOut(0),
      i2fPolyGap    = true,
      preCalcGap    = true,
      tableCalcGap  = true,
      calcPostGap   = true,
      mulStage = PipelineStageConfig.atOut(1)
    )

  runtest(cfgBF16,   1000, "Test BF16")
}
