package rial.tests

import chisel3._
import chisel3.simulator.scalatest.ChiselSim

import circt.stage.ChiselStage

import org.scalatest.funspec.AnyFunSpec

import spire.math.SafeLong
import spire.implicits._

import rial.math._
import rial.math.FuncKind._
import rial.arith._
import rial.table.PolynomialCoefficients
import rial.util.PipelineStageConfig
import rial.util.ScalaUtil._

import scala.collection.mutable.Queue
import scala.util.Random

/** A [[rial.math.MathFuncConfig]] whose tables do not share one address
 *  width: what each table is built at, what the software models give through
 *  it, and what the hardware built from it computes.
 */
class TableAddrWidthTest extends AnyFunSpec with ChiselSim {

  val r = new Random(20260930)

  val spec = RealSpec.Float32Spec

  /** the address width the table of each function is given here; every
   *  function left out keeps the default.
   */
  val widths = Map(
    Sqrt        -> 5,
    InvSqrt     -> 6,
    Reciprocal  -> 7,
    Exp         -> 6,
    Sin         -> 6,
    Log         -> 6,
    ACosPhase2  -> 7,
    ATan2Phase2 -> 6,
  )

  val defaultW = 8

  /** how many sub-intervals of its own each table holds per address, i.e. how
   *  many regions it gives the argument.
   */
  val regions = Map(
    Sqrt        -> 2, // one region per parity of the exponent
    InvSqrt     -> 2, // the same
    ACosPhase1  -> 2, // reads the table of Sqrt
    Log         -> 4, // normal, small negative, small positive and padding
    Reciprocal  -> 1,
    Sin         -> 1,
    Cos         -> 1,
    Exp         -> 1,
    ACosPhase2  -> 1,
    ATan2Phase1 -> 1,
    ATan2Phase2 -> 1,
  )

  private def config(addrW: TableAddrWidth,
    funcs: Seq[FuncKind] = MathFuncConfig.standard): MathFuncConfig = {
    new MathFuncConfig(
      funcs   = funcs,
      spec    = spec,
      nOrder  = 2,
      addrW   = addrW,
      inBitW  = spec.manW - defaultW,
      outBitW = spec.manW + 3,
    )
  }

  /** the `standard` set, each table at the width the sweep found for it. */
  val cfg = config(PerFuncAddrWidth(defaultW, widths))

  private def generateRealWithin(p: Double, spec: RealSpec, r: Random) = {
    val x = new RealGeneric(spec, (r.nextDouble() - 0.5) * p)
    new RealGeneric(spec,
      (x.value & (maskSL(spec.exW + 1) << spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  private def generateRealFull(spec: RealSpec, r: Random) = {
    new RealGeneric(spec, SafeLong(BigInt(spec.W, r)))
  }

  // ==========================================================================
  // The tables the configuration builds
  // ==========================================================================

  describe("a per-function table address width") {

    it("builds the table of each function at the width that function was given") {
      for (fn <- cfg.funcs) {
        val own = TableAddrWidth.owner(fn)
        assert(cfg.bitWOf(fn).addrW == widths.getOrElse(own, defaultW),
          f"the table of ${FuncKind.getString(fn)} is not at the width of " +
          f"${FuncKind.getString(own)}")
      }
    }

    it("gives a function that reads another one's table the owner's width") {
      assert(cfg.bitWOf(ACosPhase1 ).addrW == cfg.bitWOf(Sqrt      ).addrW)
      assert(cfg.bitWOf(ATan2Phase1).addrW == cfg.bitWOf(Reciprocal).addrW)
      assert(cfg.bitWOf(Cos        ).addrW == cfg.bitWOf(Sin       ).addrW)
    }

    it("holds as many rows in each table as its own width asks for") {
      for (fn <- cfg.funcs) {
        assert(cfg.coeffs(fn).coeffs.length ==
          regions(fn) * (1 << cfg.bitWOf(fn).addrW),
          f"the table of ${FuncKind.getString(fn)} does not hold the rows " +
          f"an address of ${cfg.bitWOf(fn).addrW}%d bits reads")
      }
    }

    it("reports the widest address of them all as the merged one") {
      assert(cfg.bitW.addrW == cfg.funcs.map(fn => cfg.bitWOf(fn).addrW).max)
      assert(cfg.bitW.addrW == 7)
    }

    it("refuses a width given to a function that reads another one's table") {
      assertThrows[IllegalArgumentException] {
        PerFuncAddrWidth(defaultW, Map(Cos -> 6))
      }
      assertThrows[IllegalArgumentException] {
        PerFuncAddrWidth(defaultW, Map(ACosPhase1 -> 5))
      }
      assertThrows[IllegalArgumentException] {
        PerFuncAddrWidth(defaultW, Map(ATan2Phase1 -> 7))
      }
    }

    it("takes a plain number as one width for every table") {
      val uniform = config(defaultW)
      for (fn <- uniform.funcs) {
        assert(uniform.bitWOf(fn).addrW == defaultW)
      }
    }
  }

  // ==========================================================================
  // The software models through it
  // ==========================================================================

  /** The table of `fn`, built by a configuration of `funcs` alone and read at
   *  the widths the shared polynomial evaluator of `cfg` reads it at: the
   *  table `cfg` should hold for `fn`, with none of the other functions of
   *  `cfg` around it.
   */
  private def alone(fn: FuncKind, addrW: TableAddrWidth,
    funcs: Seq[FuncKind]): PolynomialCoefficients = {
    config(addrW, funcs).tables(fn).generateTable(cfg.bitWOf(fn))
  }

  describe("the software model of a per-function address width") {

    val sqrtAlone  = alone(Sqrt,       5, Seq(Sqrt))
    val recipAlone = alone(Reciprocal, 7, Seq(Reciprocal))

    val acosWidths = PerFuncAddrWidth(defaultW, Map(Sqrt -> 5, ACosPhase2 -> 7))
    val acosFuncs  = FuncKind.normalize(Seq(ACosPhase2))
    val acosAlone1 = alone(Sqrt,       acosWidths, acosFuncs)
    val acosAlone2 = alone(ACosPhase2, acosWidths, acosFuncs)

    val n = 2000

    it("gives sqrt(x) what a configuration of sqrt alone gives") {
      for (i <- 1 to n) {
        val x = generateRealFull(spec, r)
        assert(SqrtSimGeneric(cfg.coeffs(Sqrt), x).value ==
               SqrtSimGeneric(sqrtAlone, x).value,
          f"sqrt(${x.toDouble}) differs from the one sqrt alone gives")
      }
    }

    it("gives 1/x what a configuration of the reciprocal alone gives") {
      for (i <- 1 to n) {
        val x = generateRealFull(spec, r)
        assert(ReciprocalSimGeneric(cfg.coeffs(Reciprocal), x).value ==
               ReciprocalSimGeneric(recipAlone, x).value,
          f"1/${x.toDouble} differs from the one the reciprocal alone gives")
      }
    }

    it("gives acos(x) what a configuration of acos alone gives") {
      for (i <- 1 to n) {
        val x  = generateRealWithin(2.0, spec, r)
        val z  = ACosPhase2SimGeneric(cfg.coeffs(ACosPhase2),
          ACosPhase1SimGeneric(cfg.coeffs(Sqrt), x))
        val z0 = ACosPhase2SimGeneric(acosAlone2,
          ACosPhase1SimGeneric(acosAlone1, x))
        assert(z.value == z0.value,
          f"acos(${x.toDouble}) differs from the one acos alone gives")
      }
    }
  }

  // ==========================================================================
  // The hardware built from it
  // ==========================================================================

  /** one function of `cfg` through the hardware, against its software model.
   *
   *  @param fn        the function to issue
   *  @param stage     where the pipeline registers go
   *  @param name      what the case is called
   *  @param generator the arguments to issue
   *  @param reference the software model of `fn` through `cfg`
   */
  private def runtest(fn: FuncKind, stage: PipelineStageConfig, name: String,
    generator: (RealSpec, Random) => RealGeneric,
    reference: RealGeneric => RealGeneric): Unit = {

    val total = stage.resolve(MathFunctions.stagePoints(cfg)).total

    it(f"${FuncKind.pointName(fn)}(x) ${name}, " +
       f"pipereg ${MathTestUtil.label(cfg, stage)}") {
      val n = 1000
      simulate(new MathFunctions(cfg, stage)) { c =>
        val q = new Queue[(BigInt, BigInt)]
        for (i <- 1 to n + total) {
          val x = generator(spec, r)
          q += ((x.value.toBigInt, reference(x).value.toBigInt))

          c.io.in.valid.poke(true.B)
          c.io.in.bits.sel.poke(cfg.signal(fn))
          c.io.in.bits.x.poke(x.value.toBigInt.U(spec.W.W))

          val zi = c.io.out.bits.z.peek().litValue.toBigInt
          if (i > total) {
            assert(c.io.out.valid.peek().litToBoolean,
              f"io.out.valid is not set at cycle $i%d (latency $total%d)")
            val (xid, z0) = q.dequeue()
            assert(zi == z0,
              f"${FuncKind.getString(fn)}(0x${xid}%08x): " +
              f"0x${zi}%08x != 0x${z0}%08x")
          }
          c.clock.step(1)
        }
      }
    }
  }

  describe("the hardware of a per-function address width") {

    val noPipe    = PipelineStageConfig.none
    val allPoints = MathTestUtil.allPoints(cfg)

    it("elaborates every function of the configuration") {
      ChiselStage.emitCHIRRTL(new MathFunctions(cfg, allPoints))
    }

    val sqrtRef  = SqrtSimGeneric      (cfg.coeffs(Sqrt      ), _: RealGeneric)
    val recipRef = ReciprocalSimGeneric(cfg.coeffs(Reciprocal), _: RealGeneric)

    runtest(Sqrt, noPipe, "over the whole format", generateRealFull(_, _), sqrtRef)
    runtest(Sqrt, allPoints, "over the whole format", generateRealFull(_, _), sqrtRef)

    runtest(Reciprocal, noPipe, "over the whole format",
      generateRealFull(_, _), recipRef)
    runtest(Reciprocal, allPoints, "over the whole format",
      generateRealFull(_, _), recipRef)
  }
}
