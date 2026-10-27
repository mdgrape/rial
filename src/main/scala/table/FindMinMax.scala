package rial.table

import org.hipparchus.analysis.solvers.LaguerreSolver

object FindMinMax {

  /** The minimum and the maximum of a polynomial on `[xMin, xMax]`.
   *
   * {{{
   * // x^3 - x on [-1, 1]
   * val (lo, hi) = FindMinMax.polynomial(Seq(0.0, -1.0, 0.0, 1.0), -1.0, 1.0)
   * }}}
   *
   * @param coeffs the coefficients, the constant term first
   */
  def polynomial(
    coeffs: Seq[Double],
    xMin: Double,
    xMax: Double
  ): (Double, Double) = {
    require(coeffs.nonEmpty)
    require(xMin <= xMax)

    val deriv = (1 until coeffs.length).map(i => { coeffs(i) * i })

    // the real part of a complex root is not an extremum, but it lies inside
    // the range of the polynomial like any other point of the interval does,
    // and it keeps a real root a solver returns as nearly real
    val inside = roots(deriv).filter(x => { xMin < x && x < xMax })

    val ys = (Seq(xMin, xMax) ++ inside).map(x => { eval(coeffs, x) })
    (ys.min, ys.max)
  }

  private def eval(coeffs: Seq[Double], x: Double): Double = {
    coeffs.init.foldRight(coeffs.last)((c, sum) => { sum * x + c })
  }

  // the real parts of the roots, the constant term first
  private def roots(coeffs: Seq[Double]): Seq[Double] = {
    val cs = coeffs.reverse.dropWhile(_ == 0.0).reverse
    if (cs.length <= 1) {
      Seq.empty
    } else if (cs.length == 2) {
      Seq(-cs(0) / cs(1))
    } else if (cs.length == 3) {
      val (c, b, a) = (cs(0), cs(1), cs(2))
      val disc = b * b - 4.0 * a * c
      if (disc < 0.0) {
        Seq(-b / (2.0 * a))
      } else {
        // the root of larger magnitude first, then the other one from it,
        // so that neither cancels
        val q = -0.5 * (b + (if (b < 0.0) { -1.0 } else { 1.0 }) * math.sqrt(disc))
        if (q == 0.0) {
          Seq(0.0)
        } else {
          Seq(q / a, c / q)
        }
      }
    } else {
      new LaguerreSolver().solveAllComplex(cs.toArray, 0.0).toSeq.map(_.getReal)
    }
  }
}
