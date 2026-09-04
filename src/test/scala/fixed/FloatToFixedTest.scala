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
// Testing Float -> Fixed using ChiselTest
//

class FloatToFixedTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 10000)"))
  def n: Int = getOption[Int]("n").getOrElse(10000)
  val r = new Random(123456789)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def convTest(xSpec : RealSpec, bp : Option[(Int, Int)], zSpec : FixedSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig, range : Double, checkHold : Boolean = false ) = {
      simulate(new FloatToFixedGeneric( xSpec, bp.map(x => x._1), zSpec, roundSpec, stage )) { c =>
        {
          // the pipeline registers advance while in.valid is high; it stays
          // high for the whole run but for the hold check at the end
          c.io.in.valid.poke(true.B)
          var q  = new Queue[(BigInt,BigInt)]
          val fracW  = if(bp.isEmpty) {zSpec.fracW} else {bp.get._2}
          val delta  = pow(2.0, -fracW)
          val nstage = c.getStage
          val gen = generateRealWithin(range,_,_)
          // clock edges since the start of the run: out.valid is in.valid
          // delayed by nstage of them
          var cycle = 0
          // println(f"Test Within (-${range}, ${range})")
          for(i <- 1 to n+nstage) {
            val x0 = gen(xSpec, r).toDouble
            val z  = (x0 / delta).toLong
            val x  = new RealGeneric(xSpec, z * delta)

            val xi = x.value.toBigInt
            val ziref = if(z >= 0) { z.toBigInt } else { (z & maskL(zSpec.W)).toBigInt }

            q += ((xi, ziref))
            c.io.in.bits.x.poke(xi.U(xSpec.W.W))
            if (bp.isDefined) {
              val yW   = bp.get._1
              val yval = bp.get._2
              c.io.in.bits.y.get.poke(yval.U(yW.W))
            }
            val zi = c.io.out.bits.peek().litValue.toBigInt
            assert(c.io.out.valid.peek().litToBoolean == (cycle >= nstage),
              f"out.valid is ${c.io.out.valid.peek().litToBoolean}%b after $cycle%d " +
              f"cycle(s) of a $nstage%d stage pipeline")
            c.clock.step(1)
            cycle += 1
            if (i > nstage) {
              val (xid,zid) = q.dequeue()
              assert(zi == zid, f"x=${xid.toLong.toBinaryString}(${new RealGeneric(xSpec, xid).toDouble}) actual(${zi.toLong.toBinaryString}) != expect(${zid.toLong.toBinaryString})")
            }
          }
          q.clear()

          // With valid low the pipeline registers do not advance, so the
          // output holds whatever the datapath carried, whatever the input does.
          if (checkHold) {
            assert(nstage > 0, "the hold check needs at least one pipeline register")
            def pokeX() = {
              val xi = gen(xSpec, r).value.toBigInt
              c.io.in.bits.x.poke(xi.U(xSpec.W.W))
            }
            // fill every stage with the same input, so draining the pipeline
            // does not move the output by itself
            pokeX()
            c.clock.step(nstage)
            val zHold = c.io.out.bits.peek().litValue

            c.io.in.valid.poke(false.B)
            pokeX()
            for(i <- 1 to nstage+2) {
              c.clock.step(1)
              assert(c.io.out.bits.peek().litValue == zHold,
                f"out.bits moved while in.valid was low after $i%d cycle(s)")
            }
            // the low valid has travelled through the whole chain by now
            assert(!c.io.out.valid.peek().litToBoolean,
              "out.valid is still high after in.valid was low for the whole latency")
            c.io.in.valid.poke(true.B)
          }
        }
      }
  }

  describe( "Test float -> fixed") {
    it(f"Converter signed Q24 to Float with pipereg 0") {
      convTest( RealSpec.Float32Spec, None, new FixedSpec(32, 24, true, true),
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-24-1))
    }
    it(f"Converter signed Qd24 to Float with pipereg 0") {
      convTest( RealSpec.Float32Spec, Some((5, 24)), new FixedSpec(32, 24, true, true),
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-24-1))
    }
    it(f"Converter signed Qd20 to Float with pipereg 0") {
      convTest( RealSpec.Float32Spec, Some((5, 20)), new FixedSpec(32, 24, true, true),
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-20-1))
    }
    it(f"Converter signed Q20 to Float with pipereg 0") {
      convTest( RealSpec.Float32Spec, None, new FixedSpec(32, 20, true, true),
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-20-1))
    }

    it(f"Converter signed Q30 to Float with pipereg 0") {
      convTest( RealSpec.Float32Spec, None, new FixedSpec(32, 30, true, true),
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-30-1))
    }
    it(f"Converter signed Q16 to Float with pipereg 0") {
      convTest( RealSpec.Float32Spec, None, new FixedSpec(32, 16, true, true),
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-16-1))
    }
    it(f"Converter signed int32 to Float with pipereg 0") {
      convTest( RealSpec.Float32Spec, None, new FixedSpec(32, 0, true, true),
        RoundSpec.roundToEven, n, PipelineStageConfig.none, pow(2.0, 32-1))
    }
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by FloatToFixedGeneric itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("Test float -> fixed with pipeline registers") {
    // the configuration that also carries the hold check of `valid`
    val holdStage = PipelineStageConfig.specify("unpack" -> 1, "align" -> 1, "round" -> 1)
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.default(3),
      holdStage)) {
      it(f"Converter signed Q24 to Float with pipereg [${stage.getString}]") {
        convTest( RealSpec.Float32Spec, None, new FixedSpec(32, 24, true, true),
          RoundSpec.roundToEven, n, stage, pow(2.0, 32-24-1),
          checkHold = (stage == holdStage))
      }
    }
    it(f"Converter BFloat16 to signed Q24 with pipereg [${PipelineStageConfig.default(3).getString}]") {
      convTest( RealSpec.BFloat16Spec, None, new FixedSpec(32, 24, true, true),
        RoundSpec.roundToEven, n, PipelineStageConfig.default(3), pow(2.0, 32-24-1))
    }
  }
}

