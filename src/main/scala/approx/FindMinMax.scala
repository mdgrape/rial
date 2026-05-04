package rial.approx

import org.hipparchus.optim.MaxEval
import org.hipparchus.optim.nonlinear.scalar.GoalType
import org.hipparchus.optim.univariate.{BrentOptimizer, SearchInterval, UnivariateObjectiveFunction}

object FindMinMax {

  def findMin(
    f: Double => Double,
    xMin: Double,
    xMax: Double
  ): Double = {
    require(isUnimodal(f, xMin, xMax))
    val optimizer = new BrentOptimizer(1e-10, 1e-14)
    val minResult = optimizer.optimize(
      new MaxEval(1000),
      new UnivariateObjectiveFunction((x: Double) => f(x)),
      GoalType.MINIMIZE,
      new SearchInterval(xMin, xMax)
    )
    minResult.getValue
  }

  def findMax(
    f: Double => Double,
    xMin: Double,
    xMax: Double
  ): Double = {
    require(isUnimodal(f, xMin, xMax))
    val optimizer = new BrentOptimizer(1e-10, 1e-14)
    val maxResult = optimizer.optimize(
      new MaxEval(1000),
      new UnivariateObjectiveFunction((x: Double) => f(x)),
      GoalType.MAXIMIZE,
      new SearchInterval(xMin, xMax)
    )
    maxResult.getValue
  }

  // check it approximately (since it is a part of polynomial approx, user
  // never pass a function with high-frequency component. it is enough (or
  // the whole process becomes nonsense)
  def isUnimodal(
    f: Double => Double, xMin: Double, xMax: Double,
    n: Int = 1000
  ): Boolean = {
    require(n >= 2)

    val xs = (0 to n).map(i => { xMin + (xMax - xMin) * i / n })
    val ys = xs.map(x => f(x))
    val diffs = (1 until ys.length).map(i => { ys(i) - ys(i-1) })
    val signChanges = (1 until diffs.length).count(i => {
      diffs(i-1) * diffs(i) < 0.0
    })
    signChanges <= 1
  }
}
