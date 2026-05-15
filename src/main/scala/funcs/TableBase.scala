package rial.funcs

import rial.approx.{PolynomialBitWidth, PolynomialCoefficients}

private[rial] trait TableBase {
  def getBitWidth(): PolynomialBitWidth
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients
}
