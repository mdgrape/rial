package rial.tests

import org.scalatest._
import org.scalatest.funspec.AnyFunSpec

import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import rial.arith._
import rial.fixed._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue

import spire.math.SafeLong
import spire.implicits._

//
// Testing Fixed -> Float using ChiselTest
//

class ScaledFixedToFloatTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 10000)"))
  def n: Int = getOption[Int]("n").getOrElse(10000)
  val r = new Random(123456789)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def convTest(xSpec : FixedSpec, ySpec : RealSpec, zSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig, xrange : Double, yrange : Double ) = {
      simulate(new ScaledFixedToFloatGeneric( xSpec, ySpec, zSpec, roundSpec, stage )) { c =>
        {
          var q  = new Queue[(BigInt, BigInt, BigInt)]
          val delta  = pow(2.0, -xSpec.fracW)
          val nstage = c.getStage
          val gen_x = generateRealWithin(xrange,_,_)
          val gen_y = generateRealWithin(yrange,_,_)
          // println(f"Test x within (-${xrange}, ${xrange})")
          for(i <- 1 to n+nstage) {
            val x0 = gen_x(zSpec, r).toDouble / delta
            val xi = if(x0 >= 0) {
              x0.toBigInt
            } else {
              (x0.toLong & maskL(xSpec.W)).toBigInt
            }
            val y  = gen_y(ySpec, r)
            val yi = y.value.toBigInt
            val z  = new RealGeneric(zSpec, ((x0.toLong * delta) * y.toDouble))
            val zi0 = (z.value.toLong & maskL(zSpec.W)).toBigInt

            // println(f"x = ${xi}(${x0 * delta}), y = ${y.toDouble}, z = ${z.toDouble}")

            q += ((xi,yi,zi0))
            c.io.x.poke(xi.U(xSpec.W.W))
            c.io.y.poke(yi.U(ySpec.W.W))
            val zi = c.io.z.peek().litValue.toBigInt
            c.clock.step(1)
            if (i > nstage) {
              val (xid,yid,zid) = q.dequeue()

              val zisgn = bit(zSpec.W-1,  zi)
              val ziex  = slice(zSpec.manW, zSpec.exW,  zi)
              val ziman = slice(0,          zSpec.manW, zi)
              val zdsgn = bit(zSpec.W-1,  zid)
              val zdex  = slice(zSpec.manW, zSpec.exW,  zid)
              val zdman = slice(0,          zSpec.manW, zid)

              assert(zi == zid, f"x=${xid.toLong.toBinaryString} ${zisgn}|${ziex}|${ziman.toLong.toBinaryString} != ${zdsgn}|${zdex}|${zdman.toLong.toBinaryString}")
            }
          }
          q.clear()
        }
      }
  }

  describe("Test fixed -> float") {
    it(f"Scale signed Q24 by large Float with pipereg 0") {
      convTest( new FixedSpec(32, 24, true, true), RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-24-1), pow(2.0, 128.0))
    }
    it(f"Scale signed Q24 by comparable Float with pipereg 0") {
      convTest( new FixedSpec(32, 24, true, true), RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-24-1), pow(2.0, 7.0))
    }

    it(f"Scale signed Q30 by large Float with pipereg 0") {
      convTest( new FixedSpec(32, 30, true, true), RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-30-1), pow(2.0, 128.0))
    }
    it(f"Scale signed Q30 by comparable Float with pipereg 0") {
      convTest( new FixedSpec(32, 30, true, true), RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-30-1), pow(2.0, 7.0))
    }

    it(f"Scale signed Q16 by large Float with pipereg 0") {
      convTest( new FixedSpec(32, 16, true, true), RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-16-1), pow(2.0, 128.0))
    }
    it(f"Scale signed Q16 to comparable Float with pipereg 0") {
      convTest( new FixedSpec(32, 16, true, true), RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-16-1), pow(2.0, 7.0))
    }
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by ScaledFixedToFloatGeneric itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("Test fixed -> float with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.default(3),
      PipelineStageConfig.specify("norm" -> 1, "mul" -> 1, "round" -> 1))) {
      it(f"Scale signed Q24 by comparable Float with pipereg [${stage.getString}]") {
        convTest( new FixedSpec(32, 24, true, true), RealSpec.Float32Spec, RealSpec.Float32Spec,
          RoundSpec.roundToEven, n, stage, pow(2.0, 32-24-1), pow(2.0, 7.0))
      }
    }
    it(f"Scale signed Q24 by comparable BFloat16 with pipereg [${PipelineStageConfig.default(3).getString}]") {
      convTest( new FixedSpec(32, 24, true, true), RealSpec.BFloat16Spec, RealSpec.BFloat16Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.default(3), pow(2.0, 32-24-1), pow(2.0, 7.0))
    }
  }
}

