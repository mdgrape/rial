package rial.tests

import org.scalatest.funspec.AnyFunSpec

import chisel3._
import chisel3.simulator.scalatest.ChiselSim
import chisel3.simulator.scalatest.Cli.EmitVcd
import chisel3.simulator.scalatest.HasCliOptions.CliOption

import org.apache.commons.math3.special.Erf

import rial.arith.{RealSpec, RealGeneric}
import rial.rng._

import scala.math._

// Tail-focused chi^2 test for GenRandomGaussian (HW simulation standalone,
// checked against the analytic gaussian -- no software Sim reference).
//
// The tail of a gaussian appears with vanishing probability, so the uniform
// random sweep in GenGaussianTest never reaches the deep intervals:
// P(xex = k) = 2^-(k+1). But the generator is a deterministic function of the
// input bits, and *which* tail interval is used is selected by the position of
// the lowest set bit of the xex0 field. So we can force the input into any
// interval directly and test it with full statistical power, instead of waiting
// for the rare event. This is plain importance sampling / conditioning.
//
// Input bit layout (rndW = 64, manW = 23):
//   bit 63        : sign
//   bit 62        : page (0 = central 50%, 1 = outer/tail 50%)
//   bit [61:23]   : xex0   (xex = lowest set bit position; selects the interval)
//   bit [22:0]    : xman   (mantissa; selects x within the interval)
//
// For page = 1 with xex = k the represented uniform value x is uniform on
// [2^-(k+2), 2^-(k+1)] and the output is |z| = sqrt2 * erfcInv(x). Therefore,
// conditioned on the interval, |z| follows the gaussian *conditional* density
// restricted to that interval's z-range, with conditional CDF
//   F(t) = (Phi(t) - Phi(a)) / (Phi(b) - Phi(a)),   a <= t <= b,
// which we check by Pearson's chi^2.
class GenGaussianTailTest extends AnyFunSpec with ChiselSim with EmitVcd {

  addOption(CliOption.int("n",    "chi^2 samples per tail interval (default 5000)"))
  addOption(CliOption.int("kmin", "smallest tail interval xex to test (default 2)"))
  addOption(CliOption.int("kmax", "largest  tail interval xex to test (default xex0Width-2)"))

  val rndW   = 64
  val spec   = RealSpec.Float32Spec
  val manW   = spec.manW                  // 23
  val nOrder = 3
  // xex0 occupies bits [rndW-3 : manW]; its width sets how deep the tail goes.
  val xex0Width = rndW - 2 - manW         // 39

  private val sqrt2 = sqrt(2.0)
  private def phi(t: Double): Double = 0.5 * (1.0 + Erf.erf(t / sqrt2))

  // SplittableRandom, not java.util.Random: the 48-bit LCG correlates the two
  // 32-bit halves of a 64-bit draw (see GenGaussianTest).
  private val sr      = new java.util.SplittableRandom(0xCAFEBABEL)
  private val manMask = (BigInt(1) << manW) - 1

  // Build an input forced into tail interval xex = k (page = 1).
  //  - bit k of xex0 is the lowest set bit  -> selects the interval,
  //  - bits of xex0 above k are random      -> unused by the datapath, kept
  //                                            random so the input is realistic,
  //  - xman is a random mantissa            -> uniform x within the interval,
  //  - sign is random                       -> exercises both halves.
  private def tailInput(k: Int): BigInt = {
    val sign = BigInt(sr.nextLong() >>> 63)
    val page = BigInt(1)
    val xman = BigInt(sr.nextLong()) & manMask
    val xex0 =
      if (k >= xex0Width) BigInt(0) // all-zero xex0 -> xex = xex0Width (deepest)
      else {
        val highW = xex0Width - k - 1
        val high  = if (highW > 0) (BigInt(sr.nextLong()) & ((BigInt(1) << highW) - 1)) else BigInt(0)
        (BigInt(1) << k) | (high << (k + 1))
      }
    (sign << (rndW - 1)) | (page << (rndW - 2)) | (xex0 << manW) | xman
  }

  // z-range of tail interval xex = k. addr = k+1, so x in [2^-(k+2), 2^-(k+1)],
  // and erfcInv is decreasing, so |z| in [a, b].
  private def zRange(k: Int): (Double, Double) = {
    val a = sqrt2 * Erf.erfcInv(pow(2.0, -(k + 1))) // x_hi -> z_lo
    val b = sqrt2 * Erf.erfcInv(pow(2.0, -(k + 2))) // x_lo -> z_hi
    (a, b)
  }

  it("tail intervals follow the gaussian conditional distribution (chi^2)") {
    val n    = getOption[Int]("n").getOrElse(5000)
    val kMin = getOption[Int]("kmin").getOrElse(2)
    // all "normal" loop intervals; the deepest capped one (xex == xex0Width) is
    // covered separately below, since its x is not a single 2^-k octave.
    val kMax = getOption[Int]("kmax").getOrElse(xex0Width - 1)

    val nbins       = 16
    val minExpected = 5.0
    // reduced chi^2 per interval; generous since we take the worst over ~37
    // intervals (multiple comparisons push the max up).
    val threshold   = 3.0

    simulate(new GenRandomGaussian(rndW, spec, nOrder, 0, 0)) { c =>
      c.io.rnd.valid.poke(true.B)
      var worstChi2 = 0.0
      var worstK    = -1
      val failures  = scala.collection.mutable.ArrayBuffer.empty[String]

      for (k <- kMin to kMax) {
        val (a, b) = zRange(k)
        val denom  = phi(b) - phi(a)
        val edges  = (0 to nbins).map(i => a + i * (b - a) / nbins)
        val counts = Array.fill(nbins)(0)
        // Tolerance absorbs the 1-ULP / polynomial-approximation overshoot that
        // a sample landing exactly on an interval boundary produces; anything
        // beyond it means the interval's scale/range is genuinely wrong (a
        // mis-scaled interval lands a full octave, ~b, away -- far past tol).
        val tol = 1e-3 * (b - a)
        var outOfRange = 0

        for (_ <- 0 until n) {
          val rnd = tailInput(k)
          c.io.rnd.bits.poke(rnd.U(rndW.W))
          c.clock.step(1) // push input through the 1-cycle input Pipe
          val zi = c.io.z.bits.peek().litValue
          val az = abs(new RealGeneric(spec, zi).toDouble)
          if (az < a - tol || az >= b + tol) outOfRange += 1
          else {
            // clamp boundary samples (within tol) into the edge bins
            val bin = ((az - a) / (b - a) * nbins).toInt
            counts(min(max(bin, 0), nbins - 1)) += 1
          }
        }

        val inRange  = n - outOfRange
        var chi2     = 0.0
        var usedBins = 0
        for (i <- 0 until nbins) {
          val exp = inRange * (phi(edges(i + 1)) - phi(edges(i))) / denom
          if (exp >= minExpected) {
            val d = counts(i) - exp
            chi2 += d * d / exp
            usedBins += 1
          }
        }
        val ndeg    = max(usedBins - 1, 1)
        val redChi2 = chi2 / ndeg

        if (redChi2 > worstChi2) { worstChi2 = redChi2; worstK = k }
        // a correct interval has only a few boundary ULPs out of range; a wrong
        // scale/range puts a large fraction out, so fail above 0.5%.
        if (outOfRange > 0.005 * n)
          failures += f"xex=$k%2d |z| in [$a%7.4f, $b%7.4f): $outOfRange%d/$n out of range"
        if (redChi2 >= threshold)
          failures += f"xex=$k%2d |z| in [$a%7.4f, $b%7.4f): reduced chi^2 = $redChi2%.3f (dof=$ndeg)"
      }

      println(f"tested xex = $kMin..$kMax, $n samples each; worst reduced chi^2 = $worstChi2%.3f at xex=$worstK")
      assert(failures.isEmpty, "\n  " + failures.mkString("\n  "))
    }
  }

  // The deepest "capped" bucket, reached only when xex0 == 0 (xex = xex0Width).
  // This is the rarest input (prob 2^-xex0Width) and the part whose range was
  // wrong: it must map to the deepest tail, x in [2^-(xex0Width+manW+1),
  // 2^-(xex0Width+1)], i.e. |z| ~ 7..9 sigma, NOT |z| ~ 0. We work in x-space,
  // x_obs = erfc(|z|/sqrt2): only erfc is needed (erfcInv overflows for such
  // tiny x), it is accurate at large arguments, and the old bug (|z| ~ 0 ->
  // x_obs ~ 1) lands far outside the expected [xLo, xHi].
  //
  // We assert range + coverage only. The within-bucket shape is NOT uniform
  // (reported as the uniform-x chi^2 below): a single cubic spans ~23 octaves
  // here, so it is a coarse fit -- acceptable for a 2^-39 event, but see the
  // notes if a tighter deepest-tail is wanted.
  it("deepest capped bucket (xex0 == 0) lands in the deep tail") {
    val n     = getOption[Int]("n").getOrElse(5000)
    val nbins = 16
    val xHi   = pow(2.0, -(xex0Width + 1))            // 2^-40, near |z| ~ 7.14
    val xLo   = pow(2.0, -(xex0Width + manW + 1))     // 2^-63, near the cap ~9.1

    simulate(new GenRandomGaussian(rndW, spec, nOrder, 0, 0)) { c =>
      c.io.rnd.valid.poke(true.B)
      val counts = Array.fill(nbins)(0)
      var outOfRange = 0
      var zmin = Double.MaxValue; var zmax = Double.MinValue; var zsum = 0.0
      val tol = 1e-3 * xHi
      for (_ <- 0 until n) {
        val rnd = tailInput(xex0Width) // xex0 == 0 -> deepest bucket
        c.io.rnd.bits.poke(rnd.U(rndW.W))
        c.clock.step(1)
        val az    = abs(new RealGeneric(spec, c.io.z.bits.peek().litValue).toDouble)
        zmin = min(zmin, az); zmax = max(zmax, az); zsum += az
        val xObs  = Erf.erfc(az / sqrt2)
        if (xObs < xLo - tol || xObs >= xHi + tol) outOfRange += 1
        else {
          val bin = ((xObs - xLo) / (xHi - xLo) * nbins).toInt
          counts(min(max(bin, 0), nbins - 1)) += 1
        }
      }
      val inRange = n - outOfRange
      val exp     = inRange.toDouble / nbins
      val chi2    = counts.map(cnt => { val d = cnt - exp; d * d / exp }).sum
      val redChi2 = chi2 / (nbins - 1)
      val zLo     = sqrt2 * Erf.erfcInv(xHi) // ~7.14 sigma; lower (least deep) edge

      println(f"deepest bucket over $n samples: |z| min=$zmin%.4f mean=${zsum / n}%.4f max=$zmax%.4f")
      println(f"  expected deep-tail edge |z| >= $zLo%.4f; $inRange%d/$n in range x in [$xLo%.3e, $xHi%.3e]")
      println(f"  uniform-x bins (low x = deepest |z|): " + counts.mkString(" "))
      println(f"  uniform-x reduced chi^2 = $redChi2%.3f")

      // The actual regression for the fixed bug: the deepest bucket must land in
      // the deep tail (|z| ~ 7..9), not at |z| ~ 0 (which gives x_obs ~ 1).
      assert(outOfRange <= 0.005 * n,
        f"$outOfRange%d/$n samples outside the deep-tail range (a |z| ~ 0 bug shows as x_obs ~ 1)")
      // and it must actually span the tail, not collapse to a single value.
      assert(zmin < zLo + 0.05 && zmax > zLo + 0.5,
        f"deepest bucket does not span the tail: |z| in [$zmin%.3f, $zmax%.3f]")
    }
  }

  // Regression guard for the negative-clamp spike: at every octave's shallow
  // endpoint the cubic dips just below 0, and a negative result must saturate
  // to 0 (|z| = the octave's small-z edge), NOT wrap to ~1.0 and emit ~2x the
  // scale. The shallow endpoint is hit whenever the within-octave mantissa is 0:
  //   main octave xex=k : xex0 has its lowest set bit at k, xman = 0
  //   deep octave clz=c : xex0 == 0, xman = 1 << (manW-1-c)  (leading 1, no tail)
  // A spike (|z| ~ 2*edge) makes x_obs = erfc(|z|/sqrt2) collapse far below the
  // octave's shallow x; require it to stay close instead.
  it("no 2x-magnitude spike at octave shallow endpoints (xman mantissa = 0)") {
    simulate(new GenRandomGaussian(rndW, spec, nOrder, 0, 0)) { c =>
      c.io.rnd.valid.poke(true.B)
      val bad = scala.collection.mutable.ArrayBuffer.empty[String]

      def check(label: String, rnd: BigInt, shallowX: Double): Unit = {
        c.io.rnd.bits.poke(rnd.U(rndW.W))
        c.clock.step(1)
        val az   = abs(new RealGeneric(spec, c.io.z.bits.peek().litValue).toDouble)
        val xObs = Erf.erfc(az / sqrt2)
        // a 2x spike pushes |z| to ~2*edge, so x_obs ~ shallowX^2 << shallowX.
        if (xObs < shallowX * 0.5)
          bad += f"$label: |z|=$az%.4f x_obs=$xObs%.3e expected ~$shallowX%.3e (2x spike?)"
      }

      // main octaves: xex = k, shallow x = 2^-(k+1)
      for (k <- 0 until xex0Width) {
        val rnd = (BigInt(1) << (rndW - 2)) | (BigInt(1) << (manW + k)) // page1, xex0 bit k, xman=0
        check(f"main xex=$k", rnd, pow(2.0, -(k + 1)))
      }
      // deep octaves: xex0 == 0, leading 1 of xman at position manW-1-c,
      // shallow x = 2^-(xex0Width + 1 + c)
      for (c <- 0 until manW) {
        val rnd = (BigInt(1) << (rndW - 2)) | (BigInt(1) << (manW - 1 - c)) // page1, xex0=0
        check(f"deep clz=$c", rnd, pow(2.0, -(xex0Width + 1 + c)))
      }

      assert(bad.isEmpty, "\n  " + bad.mkString("\n  "))
    }
  }
}
