package rial.approx

import scala.math._
import scala.collection.mutable.{ArrayBuffer, ArraySeq}

class FuncTable(
  f: Double => Double,
  val nOrder: Int,
) {
  val intervals = new ArrayBuffer[FuncTableInterval]

  // sometimes, we need to split the domain of the target function
  def addInterval(
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
    nIntervals: Int,
    xMin: Double,
    xMax: Double
  ) = {
    val intervalW = (xMax - xMin) / nIntervals
    for(i <- 0 until nIntervals) {
      val xMinLocal = xMin + intervalW *  i
      val xMaxLocal = xMin + intervalW * (i+1)
      addInterval(xMinLocal, xMaxLocal)
    }
  }

  def eval(x: Double, n: Int = nOrder): Double = {
    intervals.find(iv => {
      iv.inRange(x)
    }).map( iv => {
      iv.eval(iv.mapToDomain(x), n)
    }).getOrElse(Double.NaN)
  }
}

private[rial] object FuncTable {
  def apply(
    f: Double => Double,
    nOrder: Int,
    nIntervals: Int,
    xMin: Double,
    xMax: Double,
  ): FuncTable = {
    val t = new FuncTable(f, nOrder)
    t.addEqualInterval(nIntervals, xMin, xMax)
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

  private val chebyshev = new Chebyshev(nOrder)
  private val nodes = chebyshev.nodes()
  private val funcs = nodes.map(x => f(x))

  // chebyshev weights
  private val ks = ArraySeq.fill(nOrder+1)(0.0)
  for(i <- 0 to nOrder) {

    val orderTerm = ( if(i == 0) { 1 } else { 2 } )

    ks(i) = (0 to nOrder).map(j => {
      funcs(j) * chebyshev(i, nodes(j))
    }).sum / (nOrder+1) * orderTerm
  }

  // polynomial coefficient!
  val coeffs = (0 to nOrder).map(i => {
    (0 to nOrder).map(j => {
      ks(j) * chebyshev.coeffs(j)(i)
    }).sum
  })

  // x must be in [-1, 1]
  def eval(x: Double, order: Int = nOrder): Double = {
    (0 to order).map(i => { pow(x, i) * coeffs(i) }).sum
  }
}
