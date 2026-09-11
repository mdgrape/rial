package rial.tests

import chisel3._
import chisel3.simulator.scalatest.ChiselSim

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._
import rial.arith._
import rial.util.{PipelineStageConfig, PipelineStages}

import scala.collection.mutable.ArrayBuffer
import scala.math._
import scala.util.Random

/** [[rial.math.MathFunctions]] with one latency per function: the latencies it
 *  reports, the cycles in which it lets a function be issued, and the results
 *  and tags it hands back for a stream of interleaved operations.
 */
class PerFunctionLatencyTest extends AnyFunSpec with ChiselSim {

  val r = new Random(20260910)

  /** width of the tag every configuration under test rides along with. */
  val tagW = 8

  /** how many cycles each configuration is issued into. */
  val cycles = 3000

  /** one operation: what to poke and what to expect back. */
  private case class Op(x: BigInt, y: BigInt, z: BigInt)

  private def within(spec: RealSpec, from: Double, to: Double): RealGeneric = {
    new RealGeneric(spec, from + r.nextDouble() * (to - from))
  }

  /** a random argument for the second phase of atan2, i.e. the result the
   *  first phase gives a random (y, x).
   */
  private def atan2Phase1(cfg: MathFuncConfig, spec: RealSpec):
      (RealGeneric, RealGeneric, RealGeneric) = {
    val x = within(spec, -1.0, 1.0)
    val y = new RealGeneric(spec,
      x.toDouble * within(spec, pow(2.0, -12), pow(2.0, 12)).toDouble)
    (x, y, ATan2Phase1SimGeneric(cfg.coeffs(Reciprocal), y, x))
  }

  /** a random operation of each configured function, with the result its
   *  software model gives it -- the same generators and models the tests of
   *  the single functions use.
   */
  private def operations(cfg: MathFuncConfig, spec: RealSpec):
      Map[FuncKind, () => Op] = {
    val logCoeffs = if (cfg.has(Log)) {
      val t = cfg.tables(Log).asInstanceOf[LogTable]
      Some((t.genNormalCoeffs(cfg.bitW), t.genSmallPosCoeffs(cfg.bitW),
            t.genSmallNegCoeffs(cfg.bitW)))
    } else {
      None
    }

    def op(x: RealGeneric, z: RealGeneric): Op = {
      Op(x.value.toBigInt, BigInt(0), z.value.toBigInt)
    }

    cfg.funcs.map(fn => fn -> (() => fn match {
      case Sqrt => {
        val x = within(spec, -128.0, 128.0)
        op(x, SqrtSimGeneric(cfg.coeffs(Sqrt), x))
      }
      case InvSqrt => {
        val x = within(spec, -128.0, 128.0)
        op(x, InvSqrtSimGeneric(cfg.coeffs(InvSqrt), x))
      }
      case Reciprocal => {
        val x = within(spec, -128.0, 128.0)
        op(x, ReciprocalSimGeneric(cfg.coeffs(Reciprocal), x))
      }
      case Sin => {
        val x = within(spec, -2.0 * Pi, 2.0 * Pi)
        op(x, SinCosSimGeneric(true, cfg.coeffs(Sin), x))
      }
      case Cos => {
        val x = within(spec, -2.0 * Pi, 2.0 * Pi)
        op(x, SinCosSimGeneric(false, cfg.coeffs(Cos), x))
      }
      case ACosPhase1 => {
        val x = within(spec, -1.0, 1.0)
        op(x, ACosPhase1SimGeneric(cfg.coeffs(Sqrt), x))
      }
      case ACosPhase2 => {
        val x  = within(spec, -1.0, 1.0)
        val p1 = ACosPhase1SimGeneric(cfg.coeffs(Sqrt), x)
        op(p1, ACosPhase2SimGeneric(cfg.coeffs(ACosPhase2), p1))
      }
      case ATan2Phase1 => {
        val (x, y, z) = atan2Phase1(cfg, spec)
        Op(x.value.toBigInt, y.value.toBigInt, z.value.toBigInt)
      }
      case ATan2Phase2 => {
        val (_, _, p1) = atan2Phase1(cfg, spec)
        op(p1, ATan2Phase2SimGeneric(cfg.coeffs(ATan2Phase2), p1))
      }
      case Exp => {
        val x = within(spec, -87.0, 88.0)
        op(x, ExpSimGeneric(cfg.coeffs(Exp), x))
      }
      case Log => {
        val x = within(spec, 0.0, 128.0)
        val (cN, cP, cG) = logCoeffs.get
        op(x, LogSimGeneric(cN, cP, cG, x))
      }
      case other => {
        throw new IllegalArgumentException(
          s"no random operation for ${FuncKind.getString(other)}")
      }
    })).toMap
  }

  // --------------------------------------------------------------------------
  // The latency of a function as the class doc lets a user compute it: the
  // registers of its own path plus those of the parts it shares.

  private def toTable(cfg: MathFuncConfig, stages: PipelineStages,
    perFunctionLatency: Boolean, fn: FuncKind): Int = {
    val waits = cfg.funcs.exists(g => FuncKind.needPreMult(g)) &&
      (!perFunctionLatency || FuncKind.needPreMult(fn))
    val preMul = if (waits) {
      stages.latency("preMul")
    } else {
      0
    }
    stages.at("in") + preMul + stages.latency("pre." + FuncKind.pointName(fn))
  }

  private def toOutput(cfg: MathFuncConfig, stages: PipelineStages,
    perFunctionLatency: Boolean, fn: FuncKind): Int = {
    val waits = cfg.funcs.exists(g => FuncKind.needPostMult(g)) &&
      (!perFunctionLatency || FuncKind.needPostMult(fn))
    val postMul = if (waits) {
      stages.at("mulArgs") + stages.latency("postMul")
    } else {
      0
    }
    toTable(cfg, stages, perFunctionLatency, fn) +
      stages.at("table") + stages.at("coef") + stages.latency("poly") +
      postMul + stages.latency("post." + FuncKind.pointName(fn)) +
      stages.at("out")
  }

  // --------------------------------------------------------------------------

  private def runtest(name: String, cfg: MathFuncConfig, spec: RealSpec,
    stage: PipelineStageConfig, mulImpl: MulImpl,
    perFunctionLatency: Boolean, showLatencies: Boolean = false) = {

    it(name) {
      val ops = operations(cfg, spec)
      simulate(new MathFunctions(cfg, spec, stage, mulImpl, mulImpl, mulImpl,
        perFunctionLatency, tagW)) { c =>

        val dTab = cfg.funcs.map(fn => {
          fn -> toTable(cfg, c.stages, perFunctionLatency, fn)
        }).toMap
        val dOut = cfg.funcs.map(fn => {
          fn -> toOutput(cfg, c.stages, perFunctionLatency, fn)
        }).toMap

        cfg.funcs.foreach(fn => {
          assert(dOut(fn) == c.getStage(fn),
            f"${FuncKind.getString(fn)}%s: getStage says ${c.getStage(fn)}%d, " +
            f"the stage points say ${dOut(fn)}%d")
        })
        assert(c.getStage == dOut.values.max,
          f"getStage says ${c.getStage}%d, the longest function takes " +
          f"${dOut.values.max}%d")
        assert(c.io.issueOk.isDefined == perFunctionLatency,
          "io.issueOk is there exactly when the functions have their own " +
          "latencies")

        if (showLatencies) {
          println(f"$name%s: getStage = ${c.getStage}%d")
          cfg.funcs.foreach(fn => {
            println(f"  ${FuncKind.getString(fn)}%-12s ${c.getStage(fn)}%3d")
          })
        }

        // the operations still on their way through the module, and the
        // results they owe
        val inFlight = ArrayBuffer.empty[(FuncKind, Int)]
        val due      = ArrayBuffer.empty[(Int, BigInt, BigInt, FuncKind)]
        var issued   = 0

        /** whether `fn` may be issued this cycle: no operation in flight
         *  reaches the table input or the output selection when it would.
         */
        def free(now: Int, fn: FuncKind): Boolean = {
          !perFunctionLatency || !inFlight.exists { case (g, at) => {
            at + dTab(g) == now + dTab(fn) || at + dOut(g) == now + dOut(fn)
          }}
        }

        for (now <- 0 until cycles + dOut.values.max + 1) {
          inFlight --= inFlight.filter { case (g, at) => at + dOut(g) < now }

          // what the module answers now is decided by the operations issued
          // in earlier cycles, so ask before issuing one more
          val expectOk = cfg.funcs.map(g => g -> free(now, g)).toMap

          val fn = cfg.funcs(r.nextInt(cfg.funcs.length))
          val go = now < cycles && r.nextDouble() < 0.8 && expectOk(fn)

          if (go) {
            val op  = ops(fn)()
            val tag = BigInt(issued % (1 << tagW))
            c.io.in.valid.poke(true.B)
            c.io.in.bits.sel.poke(cfg.signal(fn))
            c.io.in.bits.x.poke(op.x.U(spec.W.W))
            c.io.in.bits.y.foreach(y => y.poke(op.y.U(spec.W.W)))
            c.io.in.bits.tag.get.poke(tag.U(tagW.W))
            inFlight += ((fn, now))
            due      += ((now + dOut(fn), tag, op.z, fn))
            issued   += 1
          } else {
            c.io.in.valid.poke(false.B)
            c.io.in.bits.sel.poke(cfg.signalNone())
            c.io.in.bits.x.poke(0.U(spec.W.W))
            c.io.in.bits.y.foreach(y => y.poke(0.U(spec.W.W)))
            c.io.in.bits.tag.get.poke(0.U(tagW.W))
          }

          c.io.issueOk.foreach(ok => {
            assert(ok(0).peek().litToBoolean,
              f"io.issueOk(signalNone) is low at cycle $now%d")
            cfg.funcs.zipWithIndex.foreach { case (g, i) => {
              assert(ok(i + 1).peek().litToBoolean == expectOk(g),
                f"io.issueOk(${FuncKind.getString(g)}%s) is " +
                f"${ok(i + 1).peek().litToBoolean}%b at cycle $now%d, " +
                f"expected ${expectOk(g)}%b")
            }}
          })

          val expected = due.filter { case (at, _, _, _) => at == now }
          if (c.io.out.valid.peek().litToBoolean) {
            assert(expected.length == 1,
              f"io.out.valid is high at cycle $now%d with " +
              f"${expected.length}%d results due")
            val (_, tag, z, g) = expected.head
            val gotTag = c.io.out.bits.tag.get.peek().litValue
            val gotZ   = c.io.out.bits.z.peek().litValue
            assert(gotTag == tag,
              f"cycle $now%d: tag ${gotTag}%d, expected ${tag}%d " +
              f"(${FuncKind.getString(g)}%s)")
            assert(gotZ == z,
              f"cycle $now%d, ${FuncKind.getString(g)}%s (tag ${tag}%d): " +
              f"z ${gotZ.toString(16)}%s, expected ${z.toString(16)}%s")
            due --= expected
          } else {
            assert(expected.isEmpty,
              f"io.out.valid is low at cycle $now%d with " +
              f"${expected.length}%d results due")
          }

          c.clock.step(1)
        }

        assert(due.isEmpty, f"${due.length}%d results never came out")
        assert(issued > cycles / 10,
          f"only $issued%d operations were issued in $cycles%d cycles")
      }
    }
  }

  // --------------------------------------------------------------------------

  val fp32 = RealSpec.Float32Spec
  val bf16 = RealSpec.BFloat16Spec

  val standardFP32 = MathFuncConfig.standardFP32
  val simpleFP32   = MathFuncConfig.simpleFP32
  val standardBF16 = MathFuncConfig.standardBF16

  /** sin/cos reduced only over the range their arguments stay in. */
  val twoPiFP32 = standardFP32.restrictSinCos(Some(SinCosDomain.twoPi(fp32)))

  /** the register placement measured to reach 500 MHz. */
  val atFiveHundredMHz = PipelineStageConfig.specify(Seq(
    "coef", "mulArgs", "poly.prod1", "poly.prod2", "poly.step2",
    "post.acosPhase1.out", "post.acosPhase1.round", "post.acosPhase2.align",
    "post.acosPhase2.neg", "post.acosPhase2.out", "post.atan2Phase1.man",
    "post.atan2Phase1.out", "post.atan2Phase2.align", "post.atan2Phase2.diff",
    "post.atan2Phase2.out", "post.exp.corr", "post.exp.out", "post.exp.round",
    "post.invSqrt.out", "post.invSqrt.round", "post.log.ex", "post.log.out",
    "post.reciprocal.out", "post.reciprocal.round", "post.sinCos.ex",
    "post.sinCos.out", "post.sqrt.out", "post.sqrt.round",
    "pre.acosPhase1.norm", "pre.acosPhase1.round", "pre.acosPhase1.sub",
    "pre.acosPhase2.align", "pre.acosPhase2.out", "pre.atan2Phase1.cmp",
    "pre.atan2Phase1.ex", "pre.atan2Phase1.out", "pre.atan2Phase2.out",
    "pre.atan2Phase2.shift", "pre.exp.neg", "pre.exp.round", "pre.exp.shift",
    "pre.invSqrt.out", "pre.log.out", "pre.reciprocal.out", "pre.sinCos.fold",
    "pre.sinCos.lead", "pre.sinCos.norm", "pre.sinCos.shift", "pre.sqrt.out",
    "table",
    "poly.step1", "poly.out", "poly.prod1.pp", "poly.prod2.pp", "preMul.pp",
    "preMul.out", "postMul.prod.pp", "postMul.prod", "postMul.out",
    "preMul.chain1").map(n => n -> 1).toMap)

  runtest("standardFP32, six registers", standardFP32, fp32,
    PipelineStageConfig.default(6), MulImpl.plain, true)

  runtest("standardFP32, a register at every boundary", standardFP32, fp32,
    MathTestUtil.allPoints(standardFP32), MulImpl.plain, true)

  runtest("standardFP32 at 500 MHz, sin/cos over (-2pi, 2pi)", twoPiFP32, fp32,
    atFiveHundredMHz, MulImpl.split(17), true, showLatencies = true)

  runtest("simpleFP32, three registers", simpleFP32, fp32,
    PipelineStageConfig.default(3), MulImpl.plain, true)

  runtest("standardBF16, four registers", standardBF16, bf16,
    PipelineStageConfig.default(4), MulImpl.plain, true)

  runtest("standardFP32, six registers, one latency for every function",
    standardFP32, fp32, PipelineStageConfig.default(6), MulImpl.plain, false)
}
