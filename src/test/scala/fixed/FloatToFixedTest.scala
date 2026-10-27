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

  // a random float in (-p, p): rD fixes the exponent, the mantissa is random
  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) : RealGeneric = {
    val rD : Double = (r.nextDouble()*2.0-1.0)*p
    val x = new RealGeneric(spec, rD)
    val y = new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
    // rD can round up to p in spec
    if (abs(y.toDouble) < p) {
      y
    } else {
      generateRealWithin(p, spec, r)
    }
  }

  // x * 2^fracW rounded to the nearest integer (ties to even), saturated to
  // the range of zSpec and masked to zSpec.W bits
  def convRef( x : RealGeneric, fracW : Int, zSpec : FixedSpec ) : BigInt = {
    val zMax = if (zSpec.signed) {
      (BigInt(1) << (zSpec.W-1)) - 1
    } else {
      (BigInt(1) << zSpec.W) - 1
    }
    val zMin = if (zSpec.signed) {
      -(BigInt(1) << (zSpec.W-1))
    } else {
      BigInt(0)
    }
    val z = if (x.isNaN) {
      zMax
    } else if (x.isInfinite) {
      if (x.sgn == 0) {
        zMax
      } else {
        zMin
      }
    } else {
      val zr = BigDecimal(Math.rint(x.toDouble * pow(2.0, fracW))).toBigInt
      zr.max(zMin).min(zMax)
    }
    z & ((BigInt(1) << zSpec.W) - 1)
  }

  // feeds the inputs one per cycle and checks each output nstage cycles later
  // against its reference
  def runConv( c : FloatToFixedGeneric, xSpec : RealSpec, bp : Option[(Int, Int)],
    xs : Seq[(BigInt, BigInt)] ) = {
    val nstage = c.getStage
    val q      = new Queue[(BigInt, BigInt)]
    val errors = new scala.collection.mutable.ArrayBuffer[String]
    if (bp.isDefined) {
      c.io.in.bits.y.get.poke(bp.get._2.U(bp.get._1.W))
    }
    // clock edges since the start of the run: out.valid is in.valid delayed
    // by nstage of them
    var cycle = 0
    for(i <- 0 until xs.length + nstage) {
      val (xi, ziref) = xs(min(i, xs.length - 1))
      if (i < xs.length) {
        q += ((xi, ziref))
      }
      c.io.in.bits.x.poke(xi.U(xSpec.W.W))
      val zi = c.io.out.bits.peek().litValue
      assert(c.io.out.valid.peek().litToBoolean == (cycle >= nstage),
        f"out.valid is ${c.io.out.valid.peek().litToBoolean}%b after $cycle%d " +
        f"cycle(s) of a $nstage%d stage pipeline")
      c.clock.step(1)
      cycle += 1
      if (i >= nstage) {
        val (xid, zid) = q.dequeue()
        if (zi != zid) {
          errors += f"x=${xid.toString(16)}(${new RealGeneric(xSpec, SafeLong(xid)).toDouble}) " +
            f"actual(${zi.toString(16)}) != expect(${zid.toString(16)})"
        }
      }
    }
    if (errors.nonEmpty) {
      fail(f"${errors.length}%d mismatch(es):\n" + errors.take(20).mkString("\n"))
    }
  }

  def convTest(xSpec : RealSpec, bp : Option[(Int, Int)], zSpec : FixedSpec, roundSpec : RoundSpec,
    n : Int, stage : PipelineStageConfig, range : Double, checkHold : Boolean = false ) = {
      simulate(new FloatToFixedGeneric( xSpec, bp.map(x => x._1), zSpec, roundSpec, stage )) { c =>
        {
          // the pipeline registers advance while in.valid is high; it stays
          // high for the whole run but for the hold check at the end
          c.io.in.valid.poke(true.B)
          val fracW  = if(bp.isEmpty) {zSpec.fracW} else {bp.get._2}
          val delta  = pow(2.0, -fracW)
          val nstage = c.getStage
          val gen = generateRealWithin(range,_,_)
          // inputs on the grid of the fixed point number
          val xs = Seq.fill(n) {
            val x0 = gen(xSpec, r).toDouble
            val z  = (x0 / delta).toLong
            val x  = new RealGeneric(xSpec, z * delta)
            (x.value.toBigInt, (z & maskL(zSpec.W)).toBigInt)
          }
          runConv(c, xSpec, bp, xs)

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

  // checks the given inputs against convRef, with rounding to even
  def rangeTest(xSpec : RealSpec, bp : Option[(Int, Int)], zSpec : FixedSpec,
    stage : PipelineStageConfig, xs : Seq[RealGeneric] ) = {
    simulate(new FloatToFixedGeneric( xSpec, bp.map(x => x._1), zSpec, RoundSpec.roundToEven, stage )) { c =>
      {
        c.io.in.valid.poke(true.B)
        val fracW = bp.map(_._2).getOrElse(zSpec.fracW)
        runConv(c, xSpec, bp, xs.map(x => (x.value.toBigInt, convRef(x, fracW, zSpec))))
      }
    }
  }

  // random floats out to 2^(intW+1), off the grid; non-negative only for an
  // unsigned zSpec
  def beyondRange(xSpec : RealSpec, bp : Option[(Int, Int)], zSpec : FixedSpec, n : Int) = {
    val intW = zSpec.W - bp.map(_._2).getOrElse(zSpec.fracW)
    Seq.fill(n) {
      val x = generateRealWithin(pow(2.0, intW+1), xSpec, r)
      if (zSpec.signed) {
        x
      } else {
        new RealGeneric(xSpec, x.value & maskSL(xSpec.exW + xSpec.manW))
      }
    }
  }

  // the edges of the range of zSpec, and the special values
  def boundaries(xSpec : RealSpec, bp : Option[(Int, Int)], zSpec : FixedSpec) = {
    val fracW = bp.map(_._2).getOrElse(zSpec.fracW)
    val intW  = zSpec.W - fracW
    val lsb   = pow(2.0, -fracW)
    def edges(p : Double) = {
      // p - lsb is the nearest float to it if it has more bits than xSpec
      Seq(p, p - lsb, Math.nextDown(p.toFloat).toDouble)
    }
    val positive = edges(pow(2.0, intW-1)) ++ Seq(pow(2.0, intW), pow(2.0, intW+3)) ++
      (if (zSpec.signed) {
        Seq()
      } else {
        edges(pow(2.0, intW))
      })
    val values = if (zSpec.signed) {
      positive ++ positive.map(-_)
    } else {
      positive
    }
    values.map(d => new RealGeneric(xSpec, d)) ++ Seq(
      RealGeneric.zero(xSpec, 0), RealGeneric.zero(xSpec, 1), RealGeneric.nan(xSpec),
      RealGeneric.inf(xSpec, 0), RealGeneric.inf(xSpec, 1))
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

  describe("Test float -> fixed at and beyond the range") {
    val specs = Seq(
      ("signed Q24 in 32",  None,          new FixedSpec(32, 24, true,  true)),
      ("signed Q14 in 18",  None,          new FixedSpec(18, 14, true,  true)),
      ("signed int32",      None,          new FixedSpec(32,  0, true,  true)),
      ("signed Qd24 in 32", Some((5, 24)), new FixedSpec(32, 24, true,  true)),
      ("unsigned Q14 in 18", None,         new FixedSpec(18, 14, false, true)))
    for ((name, bp, zSpec) <- specs) {
      it(f"Converter Float to $name beyond the range with pipereg 0") {
        rangeTest(RealSpec.Float32Spec, bp, zSpec, PipelineStageConfig.none,
          beyondRange(RealSpec.Float32Spec, bp, zSpec, n))
      }
      it(f"Converter Float to $name at the boundaries with pipereg 0") {
        rangeTest(RealSpec.Float32Spec, bp, zSpec, PipelineStageConfig.none,
          boundaries(RealSpec.Float32Spec, bp, zSpec))
      }
    }
    for ((name, bp, zSpec) <- specs.filter(_._3.W == 18)) {
      it(f"Converter Float to $name at the boundaries with pipereg [${PipelineStageConfig.default(3).getString}]") {
        rangeTest(RealSpec.Float32Spec, bp, zSpec, PipelineStageConfig.default(3),
          boundaries(RealSpec.Float32Spec, bp, zSpec))
      }
    }
  }
}
