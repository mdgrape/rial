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
// Testing FP Adder
//
class AddFPTest extends AnyFunSpec with ChiselSim with EmitVcd {

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

  def generateRealWithinPair( p : Double, xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    val x = generateRealWithin( p, xSpec, r)
    val y = generateRealWithin( p, ySpec, r)
    (x,y)
  }

  def generateRealFullPair ( xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    val x = generateRealFull(xSpec, r)
    val y = generateRealFull(ySpec, r)
    (x,y)
  }

  def generateRealZeroNonZeroPair ( xIsZero: Boolean, xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    if(xIsZero) {
      (RealGeneric.zero(xSpec), generateRealFull(ySpec, r))
    } else {
      (generateRealFull(xSpec, r), RealGeneric.zero(ySpec))
    }
  }

  def generateRealSmallDifferencePair ( p : Double, xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    val x = generateRealWithin( p, xSpec, r )
    val coef = ((r.nextDouble()) * 3.0 + 1.0) / 2.0 // [0,1) -> [0,3) -> [1,4) -> [1/2, 2)
    val yd = -1 * x.toDouble * coef

    (x, new RealGeneric(ySpec, yd))
  }

  def generateRealFiniteInfPair ( positiveInf: Boolean, p : Double, xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    val x = generateRealWithin( p, xSpec, r )
    val sgn = if(positiveInf) {0} else {1}
    val y = RealGeneric.inf(ySpec, sgn)

    (x, y)
  }
  def generateRealInfPair ( xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    val x = RealGeneric.inf(xSpec, (if(r.nextBoolean()) {0} else {1}))
    val y = RealGeneric.inf(ySpec, (if(r.nextBoolean()) {0} else {1}))
    (x, y)
  }

  def addTest(xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig, useParallelLZA : Boolean = false,
    enableDebug : Boolean = true ) = {
    simulate(new AddFPGeneric( xSpec, ySpec, zSpec, roundSpec, stage, enableDebug,
      useParallelLZA ) ) { c =>
      {
        c.debugEnableIO.foreach(_.poke(false.B))
        var q  = new Queue[(BigInt,BigInt,BigInt)]
        val nstage = c.getStage
        for (gen <- List(
          ("Test Within (-128,128)",generateRealWithinPair(128.0,_,_,_)),
          ("Test All range",generateRealFullPair(_,_,_)),
          ("Test near path result",generateRealSmallDifferencePair(128.0,_,_,_)),
          ("Test zero+y result",generateRealZeroNonZeroPair(true,_,_,_)),
          ("Test x+zero result",generateRealZeroNonZeroPair(false,_,_,_)),
          ("Test x+inf result",generateRealFiniteInfPair(true,128.0,_,_,_)),
          ("Test x-inf result",generateRealFiniteInfPair(false,128.0,_,_,_)),
          ("Test inf+/-inf result",generateRealInfPair(_,_,_))
        ) ) {
          println(gen._1)
          for(i <- 1 to n+nstage) {
            val (xr,yr) = gen._2(xSpec, ySpec, r)
            val xi = xr.value.toBigInt
            val yi = yr.value.toBigInt
            val zr = xr.add(zSpec, roundSpec, yr)
            val z0i= zr.value.toBigInt
            q += ((xi,yi,z0i))
            assert(xi < (BigInt(1) << xSpec.W))
            assert(yi < (BigInt(1) << ySpec.W))
            c.io.x.poke(xi.U(xSpec.W.W))
            c.io.y.poke(yi.U(ySpec.W.W))
            val zi = c.io.z.peek().litValue.toBigInt
            //if (zi != z0d) c.debugControlIO.poke(true.B)
            c.clock.step(1)
            if (i > nstage) {
              val (xid,yid,z0d) = q.dequeue()
              if (zi != z0d) {
                c.io.x.poke(xid.U(xSpec.W.W))
                c.io.y.poke(yid.U(ySpec.W.W))
                for(i <- 1 to nstage) c.clock.step(1)
                c.debugEnableIO.foreach(_.poke(false.B))
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

              val rsgn = bit(zSpec.W-1, z0d)
              val rex  = slice(zSpec.manW, zSpec.exW, z0d).toInt
              val rman = slice(0, zSpec.manW, z0d)

              assert(zi == z0d, f"x=${xsgn}|${xex-xSpec.exBias}|${xman.toLong.toBinaryString}(${new RealGeneric(xSpec, xsgn, xex, xman).toDouble}), " +
                                f"y=${ysgn}|${yex-ySpec.exBias}|${yman.toLong.toBinaryString}(${new RealGeneric(ySpec, ysgn, yex, yman).toDouble}), " +
                                f"test=${zsgn}|${zex-zSpec.exBias}|${zman.toLong.toBinaryString}(${new RealGeneric(zSpec, zsgn, zex, zman).toDouble}) != " +
                                f"zref=${rsgn}|${rex-zSpec.exBias}|${rman.toLong.toBinaryString}(${new RealGeneric(zSpec, rsgn, rex, rman).toDouble})")

            }
          }
          q.clear()
        }
      }
    }
  }

  describe("Test x+y") {
    it(f"Add Double with pipereg 0") {
      addTest( RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
    it(f"Add Float with pipereg 0") {
      addTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
    it(f"Add Float+Double=Double with pipereg 0") {
      addTest( RealSpec.Float32Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
    it(f"Add Double+Float=Double with pipereg 0") {
      addTest( RealSpec.Float64Spec, RealSpec.Float32Spec, RealSpec.Float64Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
  }

  // Output precision different from both inputs: wider (no rounding in the
  // near path) and narrower (rounding of the near-path sum).
  describe("Test x+y with a wider or narrower output") {
    for ((xSpec, ySpec, zSpec, name) <- Seq(
      (RealSpec.BFloat16Spec, RealSpec.BFloat16Spec, RealSpec.Float32Spec, "BFloat16+BFloat16=Float"),
      (RealSpec.Float16Spec,  RealSpec.Float16Spec,  RealSpec.Float32Spec, "Half+Half=Float"),
      (RealSpec.BFloat16Spec, RealSpec.Float16Spec,  RealSpec.Float32Spec, "BFloat16+Half=Float"),
      (RealSpec.Float32Spec,  RealSpec.Float32Spec,  RealSpec.Float64Spec, "Float+Float=Double"),
      (RealSpec.Float32Spec,  RealSpec.Float32Spec,  RealSpec.BFloat16Spec, "Float+Float=BFloat16"),
      (RealSpec.Float64Spec,  RealSpec.Float64Spec,  RealSpec.Float32Spec, "Double+Double=Float"))) {
      it(f"Add $name with pipereg 0") {
        addTest( xSpec, ySpec, zSpec, RoundSpec.roundToEven, n, PipelineStageConfig.none)
      }
    }
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by AddFPGeneric itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the datapath
  describe("Test x+y with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.default(3),
      PipelineStageConfig.specify("unpack" -> 1, "swap" -> 1, "align" -> 1,
                                  "add" -> 1, "norm" -> 1, "round" -> 1),
      PipelineStageConfig.specify("in" -> 1, "swap" -> 1, "add" -> 1,
                                  "norm" -> 1, "out" -> 1))) {
      it(f"Add Float with pipereg [${stage.getString}]") {
        addTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
          RoundSpec.roundToEven, n, stage)
      }
    }
    it(f"Add BFloat16 with pipereg [${PipelineStageConfig.default(3).getString}]") {
      addTest( RealSpec.BFloat16Spec, RealSpec.BFloat16Spec, RealSpec.BFloat16Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.default(3))
    }
    // The leading zero anticipator runs in parallel with the near path
    // subtraction, so it is cut at the same points.
    // enableDebug is on: the LZA assertion that compares the anticipated shift
    // with the reference leading-zero count now excludes the exact
    // cancellation case, where no leading one exists and the shift is unused.
    it(f"Add Float with parallel LZA and pipereg [${PipelineStageConfig.default(3).getString}]") {
      addTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.default(3), true, true)
    }
    it(f"Add BFloat16+BFloat16=Float with parallel LZA and pipereg [${PipelineStageConfig.default(3).getString}]") {
      addTest( RealSpec.BFloat16Spec, RealSpec.BFloat16Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.default(3), true, true)
    }
  }
}
