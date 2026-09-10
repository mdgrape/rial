package rial.tests

import org.scalatest._

import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._
import rial.arith._
import rial.util.PipelineStageConfig

import java.io.PrintWriter
import scala.collection.mutable.ArrayBuffer

//
// Not an assertion test: drives the HW (MathFunctions) over a sweep of inputs
// and writes the I/O to space-separated text files for visual inspection.
//
//   func_io.dat  : `x sqrt invsqrt reciprocal sin cos acos exp log sigmoid softplus smg`
//   smg_io.dat   : `x smg`               (fine central sweep)
//   atan2_io.dat : `ratio atan2_xpos atan2_xneg`
//
// The functions that the `standard` config covers (sqrt..log, atan2) are
// dumped from a `standard` MathFunctions build; the all-only functions
// (sigmoid / softplus / smg) come from an `all` build.  Point `stdCfg` at
// `allFP32` if you'd rather dump everything from a single all-functions module.
//
// Run with:  sbt "testOnly rial.tests.FuncDumpTest"
//
class FuncDumpTest extends AnyFunSpec
    with ChiselSim with EmitVcd {

  val spec   = RealSpec.Float32Spec
  val stage  = PipelineStageConfig.none
  val total  = stage.total
  val outDir = sys.props.getOrElse("user.dir", ".")

  // HW configs: standard build for the standard funcs, all build for the rest.
  val stdCfg = MathFuncConfig.standardFP32
  val allCfg = MathFuncConfig.allFP32

  private def toBits(x: Double): BigInt = new RealGeneric(spec, x).value.toBigInt
  private def toDouble(bits: BigInt): Double = new RealGeneric(spec, spire.math.SafeLong(bits)).toDouble

  // Drive `inputs` (bit patterns) for function `fn` of config `cfg` through the
  // module `c`, returning the output bit pattern for each input.
  private def evalStream(c: MathFunctions, cfg: MathFuncConfig, fn: FuncKind,
                         inputs: Seq[BigInt], ys: Seq[BigInt]): Seq[BigInt] = {
    val n   = inputs.length
    val out = ArrayBuffer[BigInt]()
    for (i <- 1 to n + total) {
      val idx = i - 1
      if (idx < n) {
        c.io.in.valid.poke(true.B)
        c.io.in.bits.sel.poke(cfg.signal(fn))
        c.io.in.bits.x.poke(inputs(idx).U(spec.W.W))
        if (c.io.in.bits.y.isDefined) c.io.in.bits.y.get.poke(ys(idx).U(spec.W.W))
      } else {
        c.io.in.valid.poke(false.B)
        c.io.in.bits.sel.poke(cfg.signalNone())
        c.io.in.bits.x.poke(0.U(spec.W.W))
        if (c.io.in.bits.y.isDefined) c.io.in.bits.y.get.poke(0.U(spec.W.W))
      }
      val zi = c.io.out.bits.z.peek().litValue
      if (i > total) out += zi
      c.clock.step(1)
    }
    out.toSeq
  }

  private def evalUnary(c: MathFunctions, cfg: MathFuncConfig, fn: FuncKind,
                        inputs: Seq[BigInt]): Seq[BigInt] =
    evalStream(c, cfg, fn, inputs, Seq.fill(inputs.length)(BigInt(0)))

  it("write func_io.dat, smg_io.dat and atan2_io.dat") {
    // ---- input sweeps ----
    val xLo = -10.0; val xHi = 10.0; val nStep = 2000
    val xs    = (0 to nStep).map(i => xLo + (xHi - xLo) * i / nStep)
    val xBits = xs.map(toBits)

    // fine central sweep for SMG (transition is within |x| < ~0.03)
    val sLo = -0.06; val sHi = 0.06; val sStep = 2400
    val sxs   = (0 to sStep).map(i => sLo + (sHi - sLo) * i / sStep)
    val sBits = sxs.map(toBits)

    // atan2 ratio sweep: x>0 -> (y, x)=(r, 1); x<0 -> (y, x)=(-r, -1)
    val rLo = -8.0; val rHi = 8.0; val rStep = 1600
    val rs   = (0 to rStep).map(k => rLo + (rHi - rLo) * k / rStep)
    val posX = rs.map(_ => toBits( 1.0)); val posY = rs.map(r => toBits( r))
    val negX = rs.map(_ => toBits(-1.0)); val negY = rs.map(r => toBits(-r))

    // ---- collected results (filled inside the test blocks below) ----
    var sqrtO, invsqO, recO, sinO, cosO, acosO, expO, logO = Seq.empty[Double]
    var sigO,  spO,    smgO, smgFineO                       = Seq.empty[Double]
    var atanPosO, atanNegO                                  = Seq.empty[Double]

    // ====================================================================
    // standard build: sqrt / invsqrt / reciprocal / sin / cos / acos /
    //                 exp / log / atan2
    // ====================================================================
    simulate(new MathFunctions(stdCfg, spec, stage)) { c =>
      sqrtO  = evalUnary(c, stdCfg, Sqrt,       xBits).map(toDouble)
      invsqO = evalUnary(c, stdCfg, InvSqrt,    xBits).map(toDouble)
      recO   = evalUnary(c, stdCfg, Reciprocal, xBits).map(toDouble)
      sinO   = evalUnary(c, stdCfg, Sin,        xBits).map(toDouble)
      cosO   = evalUnary(c, stdCfg, Cos,        xBits).map(toDouble)
      acosO  = evalUnary(c, stdCfg, ACosPhase2,
                 evalUnary(c, stdCfg, ACosPhase1, xBits)).map(toDouble)
      expO   = evalUnary(c, stdCfg, Exp,        xBits).map(toDouble)
      logO   = evalUnary(c, stdCfg, Log,        xBits).map(toDouble)

      atanPosO = evalUnary(c, stdCfg, ATan2Phase2,
                   evalStream(c, stdCfg, ATan2Phase1, posX, posY)).map(toDouble)
      atanNegO = evalUnary(c, stdCfg, ATan2Phase2,
                   evalStream(c, stdCfg, ATan2Phase1, negX, negY)).map(toDouble)
    }

    // ====================================================================
    // all build: sigmoid / softplus / smg (not in `standard`)
    // ====================================================================
    simulate(new MathFunctions(allCfg, spec, stage)) { c =>
      sigO     = evalUnary(c, allCfg, Sigmoid,              xBits).map(toDouble)
      spO      = evalUnary(c, allCfg, SoftPlus,             xBits).map(toDouble)
      smgO     = evalUnary(c, allCfg, ScaleMixtureGaussian, xBits).map(toDouble)
      smgFineO = evalUnary(c, allCfg, ScaleMixtureGaussian, sBits).map(toDouble)
    }

    // ====================================================================
    // write files
    // ====================================================================
    val sw = new PrintWriter(s"$outDir/func_io.dat")
    sw.println("x sqrt invsqrt reciprocal sin cos acos exp log sigmoid softplus smg")
    for (i <- xs.indices) {
      sw.println(f"${xs(i)}%.8g ${sqrtO(i)}%.8g ${invsqO(i)}%.8g ${recO(i)}%.8g " +
                 f"${sinO(i)}%.8g ${cosO(i)}%.8g ${acosO(i)}%.8g ${expO(i)}%.8g " +
                 f"${logO(i)}%.8g ${sigO(i)}%.8g ${spO(i)}%.8g ${smgO(i)}%.8g")
    }
    sw.close()
    println(s"wrote $outDir/func_io.dat (${xs.length} rows)")

    val mw = new PrintWriter(s"$outDir/smg_io.dat")
    mw.println("x smg")
    for (i <- sxs.indices) mw.println(f"${sxs(i)}%.8g ${smgFineO(i)}%.8g")
    mw.close()
    println(s"wrote $outDir/smg_io.dat (${sxs.length} rows)")

    val aw = new PrintWriter(s"$outDir/atan2_io.dat")
    aw.println("ratio atan2_xpos atan2_xneg")
    for (i <- rs.indices) aw.println(f"${rs(i)}%.8g ${atanPosO(i)}%.8g ${atanNegO(i)}%.8g")
    aw.close()
    println(s"wrote $outDir/atan2_io.dat (${rs.length} rows)")
  }
}
