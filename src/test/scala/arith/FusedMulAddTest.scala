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
// Testing FP FMA
//
class FMAFPTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n", "number of test cases (default 1000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000)
  val r = new Random(19660809)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealExWithin( from : Double, to : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()*(to - from)) + from
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def generateRealWithinPair( p : Double, xSpec: RealSpec, ySpec: RealSpec, zSpec: RealSpec, r : Random ) = {
    val x = generateRealWithin( p, xSpec, r)
    val y = generateRealWithin( p, ySpec, r)
    val z = generateRealWithin( p, zSpec, r)
    (x,y,z)
  }

  def generateRealFullPair ( xSpec: RealSpec, ySpec: RealSpec, zSpec: RealSpec, r : Random ) = {
    val x = generateRealFull(xSpec, r)
    val y = generateRealFull(ySpec, r)
    val z = generateRealFull(zSpec, r)
    (x,y,z)
  }

  // p : x range
  // q : result scale min
  def generateRealSmallDifferencePair ( p : Double, q: Double, xSpec: RealSpec, ySpec: RealSpec, r : Random ) = {
    val x = generateRealWithin( p, xSpec, r )
    val y = generateRealWithin( p, ySpec, r )
    val scale = r.nextDouble() * (q+1.0)
    if (r.nextInt(2)==0) {
      val z =
        if (scale>=q) {
          new RealGeneric(xSpec, 0.0)
        } else {
          x.scalbn( -scale.toInt )
        }
      (y.negate().add(ySpec, RoundSpec.roundToEven,z), y)
    } else {
      val z =
        if (scale>=q) {
          new RealGeneric(ySpec, 0.0)
        } else {
          y.scalbn( -scale.toInt )
        }
      (x,x.negate().add(ySpec, RoundSpec.roundToEven,z))
    }
  }

  def generateClosePathTriple ( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = generateRealExWithin(  64.0, 128.0, xSpec, r )
    val y = generateRealExWithin(   0.5,   2.0, ySpec, r )
    val z = generateRealExWithin(-128.0, -64.0, zSpec, r )
    (x, y, z)
  }

  def generateCloseAddPathTriple ( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = generateRealExWithin(64.0, 128.0, xSpec, r )
    val y = generateRealExWithin( 0.5,   2.0, ySpec, r )
    val z = generateRealExWithin(64.0, 128.0, zSpec, r )
    (x, y, z)
  }

  def generateFarProdPathTriple ( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = generateRealExWithin( 64.0, 128.0, xSpec, r )
    val y = generateRealExWithin( 64.0, 128.0, ySpec, r )
    val z = generateRealExWithin(  1.0,   2.0, zSpec, r )
    (x, y, z)
  }

  def generateFarAddendPathTriple ( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = generateRealExWithin(  1.0,   2.0, xSpec, r )
    val y = generateRealExWithin(  1.0,   2.0, ySpec, r )
    val z = generateRealExWithin( 64.0, 128.0, zSpec, r )
    (x, y, z)
  }

  // TODO: make it efficient
  def generateEdgeCaseFloatTriple1( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = new RealGeneric(xSpec, 0, xSpec.exBias,    0)
    val y = new RealGeneric(ySpec, 0, ySpec.exBias,    maskL(ySpec.manW))
    val z = new RealGeneric(zSpec, 0, zSpec.exBias-1-ySpec.manW, maskL(zSpec.manW))
    (x, y, z)
  }
  def generateEdgeCaseFloatTriple2( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = new RealGeneric(xSpec, 0, xSpec.exBias,    0)
    val y = new RealGeneric(ySpec, 0, ySpec.exBias-1-zSpec.manW, maskL(ySpec.manW))
    val z = new RealGeneric(zSpec, 0, zSpec.exBias, maskL(zSpec.manW))
    (x, y, z)
  }
  def generateEdgeCaseFloatTriple3( xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, r : Random ) = {
    val x = new RealGeneric(xSpec, 0, xSpec.exBias-1, 0)
    val y = new RealGeneric(ySpec, 0, ySpec.exBias-1, Integer.parseInt(  "01010101010101010101011", 2))
    val z = new RealGeneric(zSpec, 0, zSpec.exBias, Integer.parseInt("11101010101010101010101", 2))
    (x, y, z)
  }

  def fmaTest(xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, wSpec : RealSpec,
    roundSpec : RoundSpec, n : Int, stage : PipelineStageConfig ) = {
    simulate(new FusedMulAddFPGeneric( xSpec, ySpec, zSpec, wSpec, roundSpec, stage, true )) { c =>
      {
        c.debugEnableIO.get.poke(false.B)
        // the pipeline registers advance while in.valid is high; it stays high
        c.io.in.valid.poke(true.B)

        var q  = new Queue[(BigInt,BigInt,BigInt,BigInt)]
        val nstage = c.getStage
        // clock edges since the start of the run: out.valid is in.valid
        // delayed by nstage of them
        var cycle = 0
        for (gen <- List(
          ("Test close path", generateClosePathTriple(_,_,_,_)),
          ("Test close add path", generateCloseAddPathTriple(_,_,_,_)),
          ("Test far prod path", generateFarProdPathTriple(_,_,_,_)),
          ("Test far addend path", generateFarAddendPathTriple(_,_,_,_)),
          ("Test Within (-128,128)",generateRealWithinPair(128.0,_,_,_,_)),
          ("Test All range",generateRealFullPair(_,_,_,_)),
          ("Test edge case",generateEdgeCaseFloatTriple1(_,_,_,_)),
          ("Test edge case",generateEdgeCaseFloatTriple2(_,_,_,_)),
          ("Test edge case",generateEdgeCaseFloatTriple3(_,_,_,_))
        ) ) {
          println(gen._1)
          for(i <- 1 to n+nstage) {
            val (xr, yr, zr) = gen._2(xSpec, ySpec, zSpec, r)

            val xi = xr.value.toBigInt
            val yi = yr.value.toBigInt
            val zi = zr.value.toBigInt

            val wr = xr.fmadd(wSpec, roundSpec, yr, zr)
            val wi = wr.value.toBigInt

            q += ((xi,yi,zi,wi))

            c.io.in.bits.x.poke(xi.U(xSpec.W.W))
            c.io.in.bits.y.poke(yi.U(ySpec.W.W))
            c.io.in.bits.z.poke(zi.U(zSpec.W.W))

            val wi0 = c.io.out.bits.peek().litValue.toBigInt
            assert(c.io.out.valid.peek().litToBoolean == (cycle >= nstage),
              f"out.valid is ${c.io.out.valid.peek().litToBoolean}%b after $cycle%d " +
              f"cycle(s) of a $nstage%d stage pipeline")

            c.clock.step(1)
            cycle += 1
            if (i > nstage) {
              val (xid,yid,zid,wid) = q.dequeue()
              if (wi0 != wid) {
                c.io.in.bits.x.poke(xid.U(xSpec.W.W))
                c.io.in.bits.y.poke(yid.U(ySpec.W.W))
                c.io.in.bits.z.poke(zid.U(zSpec.W.W))
                for(i <- 1 to nstage) c.clock.step(1)
                c.debugEnableIO.get.poke(true.B)
                c.clock.step(1)
              }

              val wi0sgn = bit(wSpec.W-1, wi0).toInt
              val wi0exp = slice(wSpec.manW, wSpec.exW, wi0)
              val wi0man = wi0 & maskSL(wSpec.manW)

              val widsgn = bit(wSpec.W-1, wid).toInt
              val widexp = slice(wSpec.manW, wSpec.exW, wid)
              val widman = wid & maskSL(wSpec.manW)
              assert(wi0 == wid, f"test(${wi0sgn}|${wi0exp}|${wi0man}) != ref(${widsgn}|${widexp}|${widman})")
            }
          }
          q.clear()
        }
      }
    }
  }

  describe("Test x*y+z") {
    it(f"FMA float32 with pipereg 0") {
      fmaTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
    it(f"FMA float64 with pipereg 0") {
      fmaTest( RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec, RealSpec.Float64Spec,
        RoundSpec.roundToEven, n, PipelineStageConfig.none)
    }
  }

  // Pipelined configurations. A representative subset only: one verilator
  // build per configuration, so this is deliberately not a sweep.
  //  - default(n) : the registers are placed by FusedMulAddFPGeneric itself
  //  - specify    : one register at every declared point, which exercises
  //                 every single pipeline register of the three paths
  describe("Test x*y+z with pipeline registers") {
    for (stage <- Seq(
      PipelineStageConfig.default(1),
      PipelineStageConfig.default(2),
      PipelineStageConfig.default(3),
      PipelineStageConfig.specify("unpack" -> 1, "mul" -> 1, "align" -> 1,
                                  "add" -> 1, "norm" -> 1, "round" -> 1))) {
      it(f"FMA float32 with pipereg [${stage.getString}]") {
        fmaTest( RealSpec.Float32Spec, RealSpec.Float32Spec, RealSpec.Float32Spec,
          RealSpec.Float32Spec, RoundSpec.roundToEven, n, stage)
      }
    }
    it(f"FMA BFloat16 with pipereg [${PipelineStageConfig.default(3).getString}]") {
      fmaTest( RealSpec.BFloat16Spec, RealSpec.BFloat16Spec, RealSpec.BFloat16Spec,
        RealSpec.BFloat16Spec, RoundSpec.roundToEven, n, PipelineStageConfig.default(3))
    }
  }

  // Every sign combination of the cases where the result is a zero whose sign
  // is decided by IEEE 754 (2019) 6.3 rather than by the datapath:
  //  - a zero product plus a zero addend, of every sign combination
  //  - a product that underflows to zero plus a zero addend
  //  - a zero product plus a normal addend
  def zeroSignTriples( spec : RealSpec ) = {
    val sgns = Seq(0, 1)
    def zero( s : Int ) = {
      new RealGeneric(spec, s, 0, 0)
    }
    def a( s : Int ) = {
      new RealGeneric(spec, s, spec.exBias+2, maskSL(spec.manW-1))
    }
    def c( s : Int ) = {
      new RealGeneric(spec, s, spec.exBias-5, maskSL(spec.manW)-1)
    }
    // the smallest normal : its square underflows to zero
    def tiny( s : Int ) = {
      new RealGeneric(spec, s, 1, 0)
    }
    val triples = for(xs <- sgns; ys <- sgns; zs <- sgns) yield {
      Seq((zero(xs), zero(ys), zero(zs)), // (+-0) * (+-0) + (+-0)
          (zero(xs), a   (ys), zero(zs)), // (+-0) * (+-a) + (+-0)
          (a   (xs), zero(ys), zero(zs)), // (+-a) * (+-0) + (+-0)
          (tiny(xs), tiny(ys), zero(zs)), // underflowing product  + (+-0)
          (zero(xs), zero(ys), c   (zs))) // (+-0) * (+-0) + (+-c)
    }
    triples.flatten
  }

  def fmaZeroSignTest( spec : RealSpec, roundSpec : RoundSpec ) = {
    simulate(new FusedMulAddFPGeneric( spec, spec, spec, spec, roundSpec,
      PipelineStageConfig.none, true )) { c =>
      {
        c.debugEnableIO.get.poke(false.B)
        // the pipeline registers advance while in.valid is high; it stays high
        c.io.in.valid.poke(true.B)

        for((xr, yr, zr) <- zeroSignTriples(spec)) {
          val xi = xr.value.toBigInt
          val yi = yr.value.toBigInt
          val zi = zr.value.toBigInt

          val wid = xr.fmadd(spec, roundSpec, yr, zr).value.toBigInt

          c.io.in.bits.x.poke(xi.U(spec.W.W))
          c.io.in.bits.y.poke(yi.U(spec.W.W))
          c.io.in.bits.z.poke(zi.U(spec.W.W))

          val wi0 = c.io.out.bits.peek().litValue.toBigInt
          // combinational (PipelineStageConfig.none): out.valid is in.valid
          assert(c.io.out.valid.peek().litToBoolean, "out.valid is low")

          assert(wi0 == wid, s"x=${xi.toString(16)} y=${yi.toString(16)} " +
            s"z=${zi.toString(16)} : test(${wi0.toString(16)}) != ref(${wid.toString(16)})")

          c.clock.step(1)
        }
      }
    }
  }

  // One verilator build per rounding mode. `truncate` is this library's
  // stand-in for roundTowardNegative, so it flips the sign of the zero.
  describe("Test the sign of a zero x*y+z") {
    for(roundSpec <- Seq(RoundSpec.roundToEven, RoundSpec.truncate)) {
      it(f"FMA float32 zero signs with ${roundSpec}") {
        fmaZeroSignTest( RealSpec.Float32Spec, roundSpec )
      }
    }
  }
}

