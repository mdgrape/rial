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
// Testing Pow2F32 using ChiselTest
//

class MultFPIntTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000)
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
      simulate(new MultFPIntGeneric( xSpec, yWidth, ySigned, zSpec, roundSpec, stage )) { c =>
        {
          // the pipeline registers advance while in.valid is high; it stays high
          c.io.in.valid.poke(true.B)
          var q  = new Queue[(BigInt,BigInt,BigInt)]
          val nstage = c.getStage
          // clock edges since the start of the run: out.valid is in.valid
          // delayed by nstage of them
          var cycle = 0
          for (gen <- List( ("Test Within (-128,128)",generateRealWithin(128.0,_,_)),
                            ("Test All range",        generateRealFull(_,_)) ) ) {
            println(gen._1)
            for(i <- 1 to n+nstage) {
              val xr = gen._2(xSpec, r)
              val yr = generateLongWithin(yWidth, ySigned, r)
              val xi = xr.value.toBigInt
              val yi = yr.toBigInt
              val zr = new RealGeneric(zSpec, xr.toDouble * yr.toDouble) // TODO write simulator
              val z0i= zr.value.toBigInt
              // println(f"x = ${xi}(${xr.toDouble}), y = ${yr}, z = ${zr.value}")
              q += ((xi,yi,z0i))
              c.io.in.bits.x.poke(xi.U(xSpec.W.W))
              if(ySigned) {
                c.io.in.bits.y.poke((yi & maskL(yWidth)).U(yWidth.W))
              } else {
                c.io.in.bits.y.poke(yi.U(yWidth.W))
              }
              val zi = c.io.out.bits.peek().litValue.toBigInt
              assert(c.io.out.valid.peek().litToBoolean == (cycle >= nstage),
                f"out.valid is ${c.io.out.valid.peek().litToBoolean}%b after $cycle%d " +
                f"cycle(s) of a $nstage%d stage pipeline")
              c.clock.step(1)
              cycle += 1
              if (i > nstage) {
                val (xid,yid,z0d) = q.dequeue()

                val zisgn  = bit(zSpec.W-1, zi).toInt
                val ziexp  = slice(zSpec.manW, zSpec.exW, zi)
                val ziman  = zi & maskSL(zSpec.manW)
                val z0dsgn = bit(zSpec.W-1, z0d).toInt
                val z0dexp = slice(zSpec.manW, zSpec.exW, z0d)
                val z0dman = z0d & maskSL(zSpec.manW)
                if(ziexp != 0) {
                  assert(zisgn == z0dsgn && ziexp == z0dexp && ziman == z0dman,
                    f"test(${zisgn}|${ziexp}|${ziman}) != ref(${z0dsgn}|${z0dexp}|${z0dman})")
                } else {
                  assert(ziexp == z0dexp && ziman == z0dman,
                    f"test(${zisgn}|${ziexp}|${ziman}) != ref(${z0dsgn}|${z0dexp}|${z0dman})")
                }
              }
            }
            q.clear()
          }
        }
      }
  }

  describe("Test x:Float * y:Int") {
    it(f"Multiplier Float * UInt8 with pipereg 0") {
      multTest( RealSpec.Float32Spec, 8, false, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
    it(f"Multiplier Float * SInt8 with pipereg 0") {
      multTest( RealSpec.Float32Spec, 8, true, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by MultFPIntGeneric itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("Test x:Float * y:Int with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.default(3),
      PipelineStageConfig.specify("unpack" -> 1, "mul" -> 1, "norm" -> 1,
                                  "round" -> 1))) {
      it(f"Multiplier Float * SInt8 with pipereg [${stage.getString}]") {
        multTest( RealSpec.Float32Spec, 8, true, RealSpec.Float32Spec,
          RoundSpec.roundToEven, n, stage)
      }
    }
    it(f"Multiplier BFloat16 * SInt8 with pipereg [${PipelineStageConfig.default(3).getString}]") {
      multTest( RealSpec.BFloat16Spec, 8, true, RealSpec.BFloat16Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.default(3))
    }
  }
}

