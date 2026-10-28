package rial.tests

import org.scalatest.funspec.AnyFunSpec

import chisel3.simulator.scalatest.HasCliOptions.CliOption
import chisel3.simulator.scalatest.ChiselSim

import org.apache.commons.math3.special.Erf

import rial.fixed.FixedSpec
import rial.rng._

import scala.math._

// Statistics of GenRandomGaussianFixedSim on uniform random words: the moments,
// the symmetry, the fractions beyond 1, 2 and 3 sigma, a chi^2 over [-5, 5),
// and the largest magnitude. The model is bit-exact to the hardware (see
// GenRandomGaussianFixedSimMatchTest), so this is the distribution the module
// produces.
class GenRandomGaussianFixedTest extends AnyFunSpec with ChiselSim {

  addOption(CliOption.int("n", "number of random words per configuration (default 1000000)"))
  def n: Int = getOption[Int]("n").getOrElse(1000000)

  private val sqrt2 = sqrt(2.0)
  private def cdf(x: Double): Double = 0.5 * (1.0 + Erf.erf(x / sqrt2))

  private def runTest(rndW: Int, zSpec: FixedSpec): Unit = {
    it(f"Q${zSpec.intW}.${zSpec.fracW} at rndW=$rndW follows the gaussian") {
      // SplittableRandom, not java.util.Random: the 48-bit LCG correlates the
      // two halves of a 64-bit draw (see GenGaussianTest)
      val sr   = new java.util.SplittableRandom(20261027L + rndW)
      val mask = (BigInt(1) << rndW) - 1
      val eps  = pow(2.0, -zSpec.fracW)

      val reach   = GenRandomGaussianFixed.reach(rndW, zSpec)
      val zs      = Array.tabulate(n)(_ => {
        val rnd = BigInt(sr.nextLong()) & mask
        GenRandomGaussianFixedSim.toDouble(
          GenRandomGaussianFixedSim.genGaussian(rnd, rndW, zSpec), zSpec)
      })

      // moments; the tail beyond the reach is missing, which is far below
      // the tolerance at this n
      val mean = zs.sum / n
      val vari = zs.map(z => (z - mean) * (z - mean)).sum / (n - 1)
      info(f"mean $mean%+.5f, variance $vari%.5f")
      assert(abs(mean) < 5.0 / sqrt(n))
      assert(abs(vari - 1.0) < 5.0 * sqrt(2.0 / n))

      // symmetry: as many positive as negative results
      val nPos = zs.count(_ > 0.0)
      val nNeg = zs.count(_ < 0.0)
      info(f"positive $nPos, negative $nNeg")
      assert(abs(nPos - nNeg) < 5.0 * sqrt(nPos + nNeg))

      // the fractions beyond 1, 2 and 3 sigma
      for (k <- 1 to 3) {
        val p     = 2.0 * (1.0 - cdf(k.toDouble))
        val count = zs.count(z => abs(z) > k)
        info(f"|z| > $k: ${count.toDouble / n}%.6f, expected $p%.6f")
        assert(abs(count - n * p) < 5.0 * sqrt(n * p * (1.0 - p)))
      }

      // chi^2 over 16 bins of [-5, 5), the bins with fewer than 5 expected
      // samples left out; the threshold is exceeded with a probability of ~3%.
      // The reference is the gaussian clipped at the reach, where the result
      // is held.
      val cdfClipped = (x: Double) => {
        if (x < -reach) { 0.0 } else if (x >= reach) { 1.0 } else { cdf(x) }
      }
      val nBins = 16
      val width = 10.0 / nBins
      val counts = Array.fill(nBins)(0L)
      zs.foreach(z => {
        val b = floor((z + 5.0) / width).toInt
        if (0 <= b && b < nBins) { counts(b) += 1 }
      })
      val terms = (0 until nBins).flatMap(b => {
        val lo  = -5.0 + b * width
        val ref = (cdfClipped(lo + width) - cdfClipped(lo)) * n
        if (ref >= 5.0) { Some((counts(b) - ref) * (counts(b) - ref) / ref) } else { None }
      })
      val chi2 = terms.sum / (terms.length - 1)
      info(f"chi^2 / dof = $chi2%.3f over ${terms.length} bins")
      assert(chi2 < 1.8)

      // no result beyond the reach, and the reach itself as often as the
      // tail beyond it rounds onto it
      val zMaxSeen = zs.map(abs).max
      val nReach   = zs.count(z => abs(z) >= reach)
      val pReach   = 2.0 * (1.0 - cdf(reach - eps / 2))
      info(f"largest |z| $zMaxSeen%.6f, reach $reach%.6f, " +
        f"${nReach} at the reach, expected ${n * pReach}%.1f")
      assert(zMaxSeen <= reach)
      assert(abs(nReach - n * pReach) < 5.0 * sqrt(n * pReach) + 1.0)
    }
  }

  runTest(64, FixedSpec(16, 12))
  runTest(32, FixedSpec(16, 13))
}
