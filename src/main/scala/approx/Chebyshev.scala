package rial.approx

import scala.math._
import scala.collection.mutable.ArraySeq

// n-th order chebyshev polynomial itself
private[rial] class Chebyshev(
  val nOrder: Int,
) {
  val coeffs = ArraySeq.fill(nOrder+1, nOrder+1)(0)

  coeffs(0)(0) = 1
  coeffs(1)(1) = 1

  for(order <- 2 to nOrder) {
    coeffs(order)(0) = -1 * coeffs(order-2)(0)
    for(i <- 1 to nOrder) {
      coeffs(order)(i) = 2 * coeffs(order-1)(i-1) - coeffs(order-2)(i)
    }
  }

  def apply(order: Int, x: Double): Double = {
    (0 to nOrder).map(i => {
      pow(x, i.toDouble) * coeffs(order)(i).toDouble
    }).sum
  }

  def nodes(): Seq[Double] = {
    (0 to nOrder).map(i => {
      cos( (i + 0.5) / (nOrder+1) * Pi )
    })
  }
}
