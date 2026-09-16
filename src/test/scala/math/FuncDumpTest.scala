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
import scala.math._

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
// The multiplier-free set (`standardFixed`) gets one file per function, all
// of them from a single build:
//
//   fixed_sqrt.dat  fixed_invsqrt.dat  fixed_reciprocal.dat  fixed_pow2.dat
//   fixed_log2.dat  fixed_sinpi.dat    fixed_cospi.dat
//   fixed_atan.dat  fixed_acos.dat
//
// Their columns are `x z ref ulp`: the argument, the result of the HW as a
// double, the double-precision reference, and the signed error of the result
// in fp32 units of the last place of that reference (0 where the two are the
// same special value). A file holds one sweep per block of lines, the blocks
// separated by a blank line: a sweep over the range that matters, then, where
// it helps, a logarithmic sweep of the distance from the argument the
// function is hardest at.
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
  val fixCfg = MathFuncConfig.standardFixedFP32

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
    simulate(new MathFunctions(stdCfg, stage)) { c =>
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
    simulate(new MathFunctions(allCfg, stage)) { c =>
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

  // ==========================================================================
  // the multiplier-free set
  // ==========================================================================

  /** `n + 1` points spread evenly over `[lo, hi]`. */
  private def linSweep(lo: Double, hi: Double, n: Int): Seq[Double] = {
    (0 to n).map(i => lo + (hi - lo) * i / n)
  }

  /** `n + 1` points spread evenly in the logarithm over `[lo, hi]`. */
  private def logSweep(lo: Double, hi: Double, n: Int): Seq[Double] = {
    (0 to n).map(i => lo * pow(hi / lo, i.toDouble / n))
  }

  /** The distances `logSweep(lo, hi, n)` laid out on both sides of `centre`,
   *  from the far end below it to the far end above it. */
  private def aroundSweep(centre: Double, lo: Double, hi: Double, n: Int)
    : Seq[Double] = {
    val d = logSweep(lo, hi, n)
    d.reverse.map(v => centre - v) ++ d.map(v => centre + v)
  }

  /** The bit pattern as a monotone integer, so that the distance of two
   *  results counts the representable values between them across zero too. */
  private def ordinal(b: BigInt): BigInt = {
    val sign = BigInt(1) << (spec.W - 1)
    if (b >= sign) { -(b - sign) } else { b }
  }

  /** The signed error of the result `z` against `ref`, in units of the last
   *  place of `ref` in `spec`; 0 where both are the same special value. */
  private def ulpOf(z: BigInt, ref: Double): Long = {
    val zd = toDouble(z)
    if (zd.isNaN && ref.isNaN) {
      0L
    } else if (zd.isInfinite && zd == ref) {
      0L
    } else {
      (ordinal(z) - ordinal(toBits(ref))).toLong
    }
  }

  // sin(pi x) and cos(pi x) of a double, kept away from the cancellation the
  // product pi x suffers where the result vanishes: the argument is folded
  // into [-1/2, 1/2] first, exactly, and the quarter of the period nearest a
  // zero of the result is read off the other function.
  private def refSinPi(x: Double): Double = {
    val r = x - 2.0 * rint(x / 2.0)
    val t = if (r > 0.5) { 1.0 - r } else if (r < -0.5) { -1.0 - r } else { r }
    if (x == 0.0) { x } else { sin(Pi * t) }
  }
  private def refCosPi(x: Double): Double = {
    val r = abs(x - 2.0 * rint(x / 2.0))
    val t = if (r > 0.5) { 1.0 - r } else { r }
    val s = if (r > 0.5) { -1.0 } else { 1.0 }
    if (t >= 0.25) { s * sin(Pi * (0.5 - t)) } else { s * cos(Pi * t) }
  }

  /** Drives `fn` over `xs` and formats one line of `x z ref ulp` per point. */
  private def sweepOf(c: MathFunctions, fn: FuncKind, xs: Seq[Double],
                      ref: Double => Double): Seq[String] = {
    val xBits = xs.map(toBits)
    val zBits = evalUnary(c, fixCfg, fn, xBits)
    xBits.zip(zBits).map { case (xi, zi) => {
      val x0 = toDouble(xi)
      val r0 = ref(x0)
      f"${x0}%.9g ${toDouble(zi)}%.9g ${r0}%.9g ${ulpOf(zi, r0)}%d"
    }}
  }

  /** Writes the sweeps of one function, a blank line between each two. */
  private def writeDump(name: String, sweeps: Seq[Seq[String]]): Unit = {
    val w = new PrintWriter(s"$outDir/$name")
    w.println("x z ref ulp")
    for ((lines, i) <- sweeps.zipWithIndex) {
      if (i != 0) { w.println("") }
      lines.foreach(w.println)
    }
    w.close()
    println(s"wrote $outDir/$name (${sweeps.map(_.length).sum} rows, " +
      s"${sweeps.length} sweeps)")
  }

  it("write the fixed_*.dat of the multiplier-free function set") {
    simulate(new MathFunctions(fixCfg, stage)) { c =>
      // the whole range of the argument, in the logarithm
      val decades = logSweep(pow(2.0, -20), pow(2.0, 20), 3000)
      writeDump("fixed_sqrt.dat",
        Seq(sweepOf(c, Sqrt,       decades, x => sqrt(x))))
      writeDump("fixed_invsqrt.dat",
        Seq(sweepOf(c, InvSqrt,    decades, x => 1.0 / sqrt(x))))
      writeDump("fixed_reciprocal.dat",
        Seq(sweepOf(c, Reciprocal, decades, x => 1.0 / x)))

      // every exponent the result can have, and the two ends where it cannot
      writeDump("fixed_pow2.dat",
        Seq(sweepOf(c, Pow2, linSweep(-130.0, 130.0, 3000), x => pow(2.0, x))))

      // the decades of x, the binade around one, and the approach to one
      val log2Ref = (x: Double) => log(x) / log(2.0)
      writeDump("fixed_log2.dat", Seq(
        sweepOf(c, Log2Fixed, logSweep(pow(2.0, -8), pow(2.0, 8), 2000), log2Ref),
        sweepOf(c, Log2Fixed, linSweep(0.5, 2.0, 2000), log2Ref),
        sweepOf(c, Log2Fixed,
          aroundSweep(1.0, pow(2.0, -24), pow(2.0, -1), 1000), log2Ref)))

      // a whole period, the approach to zero, and, for the cosine, the
      // approach to the half where its own result vanishes
      val period  = linSweep(-2.0, 2.0, 3000)
      val small   = aroundSweep(0.0, pow(2.0, -24), 0.5, 1000)
      writeDump("fixed_sinpi.dat", Seq(
        sweepOf(c, SinPiFixed, period, refSinPi),
        sweepOf(c, SinPiFixed, small,  refSinPi)))
      writeDump("fixed_cospi.dat", Seq(
        sweepOf(c, CosPiFixed, period, refCosPi),
        sweepOf(c, CosPiFixed, small,  refCosPi),
        sweepOf(c, CosPiFixed,
          aroundSweep(0.5, pow(2.0, -24), 0.5, 1000), refCosPi)))

      // the whole domain, and the approach to the zero of the result
      writeDump("fixed_atan.dat", Seq(
        sweepOf(c, ATanFixed, linSweep(-1.0, 1.0, 3000), x => atan(x)),
        sweepOf(c, ATanFixed,
          aroundSweep(0.0, pow(2.0, -24), 1.0, 1000), x => atan(x))))

      // the whole domain, and the approach to either end of it
      val nearEnd = logSweep(pow(2.0, -24), 1.0, 1000)
      writeDump("fixed_acos.dat", Seq(
        sweepOf(c, ACosFixed, linSweep(-1.0, 1.0, 3000), x => acos(x)),
        sweepOf(c, ACosFixed,
          nearEnd.reverse.map(d => -(1.0 - d)) ++ nearEnd.map(d => 1.0 - d),
          x => acos(x))))
    }
  }
}
