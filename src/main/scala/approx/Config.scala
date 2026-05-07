package rial.approx

case class ApproxFuncConfig(
  val nOrder:    Int,
  val nInterval: Int,
  val xMin:      Double,
  val xMax:      Double,
  // x is a signed int. -(1<<(xBitW-1)) means -1.0, (1<<(xBitW-1)) means 1.0
  val inBitW:    Int, // width of x in each interval polynomial.
  val outBitW:   Int, // (1 << outBitW) means 1.0
) {
  val width = (xMax - xMin)
  val intervalW = width / nInterval
  val rIntervalW = 1.0 / intervalW
}

