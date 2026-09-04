package rial.tests

import org.scalatest._
import org.scalatest.funspec.AnyFunSpec

import chisel3._
import chisel3.util._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import rial.arith._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.math._
import scala.util.Random
import scala.collection.mutable.Queue

import spire.math.SafeLong
import spire.implicits._

//
// Testing mul using ChiselTest
//
class MultFPTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000)
  val r = new Random(19660809)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def multTest(xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig ) = {
      simulate(new MultFPGeneric( xSpec, ySpec, zSpec, roundSpec, stage )) { c =>
        {
          // the pipeline registers advance while in.valid is high; it stays high
          c.io.in.valid.poke(true.B)
          var q  = new Queue[(BigInt,BigInt,BigInt)]
          val nstage = c.getStage
          // clock edges since the start of the run: out.valid is in.valid
          // delayed by nstage of them
          var cycle = 0
          for (gen <- List( ("Test Within (-128,128)",generateRealWithin(128.0,_,_)),
                            ("Test All range",generateRealFull(_,_)) ) ) {
            println(gen._1)
            for(i <- 1 to n+nstage) {
              val xr = gen._2(xSpec, r)
              val yr = gen._2(ySpec, r)
              val xi = xr.value.toBigInt
              val yi = yr.value.toBigInt
              val zr = xr.multiply(zSpec, roundSpec, yr)
              val z0i= zr.value.toBigInt
              q += ((xi,yi,z0i))
              c.io.in.bits.x.poke(xi.U(xSpec.W.W))
              c.io.in.bits.y.poke(yi.U(ySpec.W.W))
              val zi = c.io.out.bits.peek().litValue.toBigInt
              assert(c.io.out.valid.peek().litToBoolean == (cycle >= nstage),
                f"out.valid is ${c.io.out.valid.peek().litToBoolean}%b after $cycle%d " +
                f"cycle(s) of a $nstage%d stage pipeline")
              c.clock.step(1)
              cycle += 1
              if (i > nstage) {
                val (xid,yid,z0d) = q.dequeue()

                val xsgn = bit(xSpec.W-1, xid)
                val xex  = slice(xSpec.manW, xSpec.exW, xid).toInt
                val xman = slice(0, xSpec.manW, xid)

                val ysgn = bit(ySpec.W-1, yid)
                val yex  = slice(ySpec.manW, ySpec.exW, yid).toInt
                val yman = slice(0, ySpec.manW, yid)

                val zsgn = bit(zSpec.W-1, zi)
                val zex  = slice(zSpec.manW, zSpec.exW, zi).toInt
                val zman = slice(0, zSpec.manW, zi)

                val rsgn = bit(zSpec.W-1, z0d)
                val rex  = slice(zSpec.manW, zSpec.exW, z0d).toInt
                val rman = slice(0, zSpec.manW, z0d)

                assert(zi == z0d, f"x=${xsgn}|${xex-xSpec.exBias}|${xman.toLong.toBinaryString}(${new RealGeneric(xSpec, xsgn, xex, xman).toDouble}), " +
                                  f"y=${ysgn}|${yex-ySpec.exBias}|${yman.toLong.toBinaryString}(${new RealGeneric(ySpec, ysgn, yex, yman).toDouble}), " +
                                  f"test=${zsgn}|${zex-zSpec.exBias}(${zex})|${zman.toLong.toBinaryString}(${new RealGeneric(zSpec, zsgn, zex, zman).toDouble}) != " +
                                  f"zref=${rsgn}|${rex-zSpec.exBias}(${rex})|${rman.toLong.toBinaryString}(${new RealGeneric(zSpec, rsgn, rex, rman).toDouble})")
              }
            }
            q.clear()
          }
        }
      }
  }

  describe("Test x*y") {
    it(f"Multiplier Double with pipereg 0") {
      multTest( RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
    it(f"Multiplier Float with pipereg 0") {
      multTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
    it(f"Multiplier Double*Float->Double with pipereg 0") {
      multTest( RealSpec.Float64Spec, RealSpec.Float32Spec, RealSpec.Float64Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
    it(f"Multiplier Float*Double->Double with pipereg 0") {
      multTest( RealSpec.Float32Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
    it(f"Multiplier Float*Float->Double with pipereg 0") {
      multTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float64Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by MultFPGeneric itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("Test x*y with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.default(3),
      PipelineStageConfig.specify("unpack" -> 1, "mul" -> 1, "norm" -> 1,
                                  "round" -> 1))) {
      it(f"Multiplier Float with pipereg [${stage.getString}]") {
        multTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
          RoundSpec.roundToEven, n, stage)
      }
    }
    it(f"Multiplier BFloat16 with pipereg [${PipelineStageConfig.default(3).getString}]") {
      multTest( RealSpec.BFloat16Spec, RealSpec.BFloat16Spec, RealSpec.BFloat16Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.default(3))
    }
  }
}

