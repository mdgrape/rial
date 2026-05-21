package rial.table

import scala.math._
import scala.collection.mutable.{ArrayBuffer, ArraySeq}

class FuncTable(
  // order must be the same for all the intervals (polynomial part is shared)
  val nOrder: Int,
) {
  val intervals = new ArrayBuffer[FuncTableInterval]

  // sometimes, we need to split the domain of the target function
  def addInterval(
    f: Double => Double,
    xMin: Double,
    xMax: Double
  ) = {
    val intervalW = xMax - xMin
    // [f_chebyshev(-1), f_chebyshev(1)] => [f_user(xMin), f_user(xMax)]
    intervals += new FuncTableInterval( (x: Double) => {
      f((x + 1.0) * 0.5 * intervalW + xMin)
    }, nOrder, xMin, xMax)
  }

  def addEqualInterval(
    f: Double => Double,
    nIntervals: Int,
    xMin: Double,
    xMax: Double
  ) = {
    val intervalW = (xMax - xMin) / nIntervals
    for(i <- 0 until nIntervals) {
      val xMinLocal = xMin + intervalW *  i
      val xMaxLocal = xMin + intervalW * (i+1)
      addInterval(f, xMinLocal, xMaxLocal)
    }
  }

  def eval(x: Double, hornerStep: Int = nOrder): Double = {
    intervals.find(iv => {
      iv.inRange(x)
    }).map( iv => {
      iv.eval(iv.mapToDomain(x), hornerStep)
    }).getOrElse(Double.NaN)
  }
}

private[rial] object FuncTable {

  // generate evenly-splitted function tables
  def apply(
    f: Double => Double,
    nOrder: Int,
    nIntervals: Int,
    xMin: Double,
    xMax: Double,
  ): FuncTable = {
    val t = new FuncTable(nOrder)
    t.addEqualInterval(f, nIntervals, xMin, xMax)
    t
  }
}

private[rial] class FuncTableInterval(
  f : Double => Double, // x in [-1, 1]
  val nOrder: Int,
  val xMin: Double,
  val xMax: Double,
) {

  def inRange(x: Double): Boolean = {
    xMin <= x && x < xMax
  }

  val width = xMax - xMin
  val rwidth = 1.0 / width
  // convert [xMin, xMax] -> [-1, 1] (domain of chebyshev)
  def mapToDomain(x: Double): Double = {
    (x - xMin) * rwidth * 2.0 - 1.0
  }

  val coeffs: Seq[Double] = if (nOrder == 0) {

    // For nOrder == 0 we skip the Chebyshev fit and just take f at the left
    // endpoint (t = -1, i.e. x = xMin).
    // Since, in case of nOrder==0, the hardware uses its whole mantissa as the
    // key of the table and that corresponds to the xMin(= -1.0 for chebyshev).
    Seq(f(-1.0))

  } else {
    // do interpolation
    val chebyshev = new Chebyshev(nOrder)
    val nodes = chebyshev.nodes()
    val funcs = nodes.map(x => f(x))

    // chebyshev weights
    val ks = ArraySeq.fill(nOrder+1)(0.0)
    for(i <- 0 to nOrder) {
      val orderTerm = ( if(i == 0) { 1 } else { 2 } )
      ks(i) = (0 to nOrder).map(j => {
        funcs(j) * chebyshev(i, nodes(j))
      }).sum / (nOrder+1) * orderTerm
    }

    // convert weights to polynomial coeff
    (0 to nOrder).map(i => {
      (0 to nOrder).map(j => {
        ks(j) * chebyshev.coeffs(j)(i)
      }).sum
    })
  }

  // x must be in [-1, 1].
  def eval(x: Double, hornerStep: Int = nOrder): Double = {
    // to determine calcBit, we use `hornerStep` cutoff with horner's rule
    //         ((c2) * x + c1) * x + c0
    // step 0.  ^^^^
    // step 1. ^^^^^^^^^^^^^^
    // step 2. ^^^^^^^^^^^^^^^^^^^^^^^^
    var sum = coeffs(nOrder)
    for(i <- 1 to hornerStep) {
      sum = sum * x + coeffs(nOrder-i)
    }
    sum
  }
}

