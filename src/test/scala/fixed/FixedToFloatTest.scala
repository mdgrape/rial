package rial.tests

import org.scalatest._
import org.scalatest.funspec.AnyFunSpec

import chisel3._
import chisel3.util._
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

class FixedToFloatTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 10000)"))
  def n: Int = getOption[Int]("n").getOrElse(10000)
  val r = new Random(123456789)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def convTest(xSpec : FixedSpec, bp: Option[(Int, Int)], zSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig, range : Double ) = {
      simulate(new FixedToFloatGeneric( xSpec, bp.map(x => x._1), zSpec, roundSpec, stage )) { c =>
        {
          var q  = new Queue[(BigInt,BigInt)]
          val fracW  = if(bp.isEmpty) {xSpec.fracW} else {bp.get._2}
          val delta  = pow(2.0, -fracW)
          val nstage = c.getStage
          val gen = generateRealWithin(range,_,_)
          // println(f"Test Within (-${range}, ${range})")
          for(i <- 1 to n+nstage) {
            val x0 = gen(zSpec, r).toDouble / delta
            val z  = new RealGeneric(zSpec, (x0.toLong * delta))
            val xi = if(x0 >= 0) {
              x0.toBigInt
            } else {
              (x0.toLong & maskL(xSpec.W)).toBigInt
            }
            q += ((xi,z.value.toBigInt))
            c.io.x.poke(xi.U(xSpec.W.W))
            if (bp.isDefined) {
              val yW   = bp.get._1
              val yval = bp.get._2
              c.io.y.get.poke(yval.U(yW.W))
            }
            val zi = c.io.z.peek().litValue.toBigInt
            c.clock.step(1)
            if (i > nstage) {
              val (xid,zid) = q.dequeue()
              assert(zi == zid, f"x=${xid.toLong.toBinaryString} ${zi.toLong.toBinaryString} != ${zid.toLong.toBinaryString}")
            }
          }
          q.clear()
        }
      }
  }

  describe("Test fixed -> float") {
    it(f"Converter signed Q24 to Float with pipereg 0") {
      convTest( new FixedSpec(32, 24, true, true), None, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-24-1))
    }
    it(f"Converter signed Qd24 to Float with pipereg 0") {
      convTest( new FixedSpec(32, 24, true, true), Some(5, 24), RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-24-1))
    }
    it(f"Converter signed Qd20 to Float with pipereg 0") {
      convTest( new FixedSpec(32, 24, true, true), Some(5, 20), RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-20-1))
    }
    it(f"Converter signed Q20 to Float with pipereg 0") {
      convTest( new FixedSpec(32, 20, true, true), None, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-20-1))
    }
    it(f"Converter int to Float with pipereg 0") {
      convTest( new FixedSpec(32, 24, true, true), Some(5, 0), RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-1))
    }

    it(f"Converter signed Q30 to Float with pipereg 0") {
      convTest( new FixedSpec(32, 30, true, true), None, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-30-1))
    }
    it(f"Converter signed Q16 to Float with pipereg 0") {
      convTest( new FixedSpec(32, 16, true, true), None, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-16-1))
    }
    it(f"Converter signed int32 to Float with pipereg 0") {
      convTest( new FixedSpec(32, 0, true, true), None, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-1))

    }
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by FixedToFloatGeneric itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("Test fixed -> float with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.default(3),
      PipelineStageConfig.specify("abs" -> 1, "norm" -> 1, "round" -> 1))) {
      it(f"Converter signed Q24 to Float with pipereg [${stage.getString}]") {
        convTest( new FixedSpec(32, 24, true, true), None, RealSpec.Float32Spec,
          RoundSpec.roundToEven, n, stage, pow(2.0, 32-24-1))
      }
    }
    it(f"Converter signed Q24 to BFloat16 with pipereg [${PipelineStageConfig.default(3).getString}]") {
      convTest( new FixedSpec(32, 24, true, true), None, RealSpec.BFloat16Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.default(3), pow(2.0, 32-24-1))
    }
  }
}
