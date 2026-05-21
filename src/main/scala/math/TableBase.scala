package rial.math

import rial.table.{PolynomialBitWidth, PolynomialCoefficients}

private[rial] trait TableBase {
  def getBitWidth(): PolynomialBitWidth
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients
}
