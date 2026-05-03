package rial.approx

import scala.math._

case class ApproxFuncConfig(
  val nOrder:    Int,
  val nInterval: Int,
  val xMin:      Double,
  val xMax:      Double,
  val fMin:      Double,
  val fMax:      Double,
  val inBitW:    Int, // input    x  width
  val outBitW:   Int, // output f(x) width
) {

  val width = (xMax - xMin)
  val intervalW = width / nInterval
  val rIntervalW = 1.0 / intervalW
}

