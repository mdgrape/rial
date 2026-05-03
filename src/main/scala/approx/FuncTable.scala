package rial.approx

import scala.math._
import scala.collection.mutable.ArraySeq

private[rial] class FuncTable(
  f: Double => Double,
  val cfg: ApproxFuncConfig,
) {
  val intervals = (0 until cfg.nInterval).map(i => {
    val xMin = cfg.xMin + cfg.intervalW * i
    val xMax = cfg.xMin + cfg.intervalW * (i+1)

    // [-1, 1] => [f(xMin), f(xMax)]
    new FuncTableInterval( (x: Double) => {
      f((x + 1.0) * 0.5 * cfg.intervalW + xMin)
    }, cfg.nOrder)
  })

  def eval(x: Double, n: Int = cfg.nOrder): Double = {
    val idx  = floor((x - cfg.xMin) * cfg.rIntervalW)
    val xMin = cfg.xMin + cfg.intervalW * idx

    // [xMin, xMax] => [-1, 1]
    intervals(idx.toInt).eval( (x - xMin) * cfg.rIntervalW * 2.0 - 1.0 )
  }
}

private[rial] class FuncTableInterval(
  f : Double => Double, // already centered to [-1, 1]
  val nOrder: Int,
) {

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

  def eval(x: Double, order: Int = nOrder): Double = {
    (0 to order).map(i => { pow(x, i) * coeffs(i) }).sum
  }
}
