//% @file GenGaussianFixedSim.scala
//
// Bit-exact software model of GenRandomGaussianFixed (see GenGaussianFixed.scala).
// Given the rndW-bit random word the hardware consumes, it returns the same
// two's complement bits. Pipeline registers do not change the result, so the
// model is a pure function of the word and the table parameters.
//
package rial.rng

import rial.fixed.FixedSpec
import rial.table._

object GenRandomGaussianFixedSim {

  /** The output bits of the generator for the random word `rnd`.
   *
   * @return the two's complement bits of the result, `zSpec.W` wide.
   */
  def genGaussian(rnd: BigInt, rndW: Int, zSpec: FixedSpec, nOrder: Int = 3,
                  splitW: Int = 3, guardBits: Int = 3,
                  posW: Option[Int] = None): BigInt = {
    val positionW = GenRandomGaussianFixed.positionWidth(zSpec, posW)
    val nOctaves  = GenRandomGaussianFixed.octaves(rndW, zSpec, posW)
    val table     = GenRandomGaussianFixed.table(zSpec, nOrder, nOctaves,
      positionW, splitW, guardBits)

    val sgn  = (rnd >> (rndW - 1)) & 1
    val page = (rnd >> (rndW - 2)) & 1
    val scan = (rnd >> positionW) & ((BigInt(1) << nOctaves) - 1)
    val pos  = rnd & ((BigInt(1) << positionW) - 1)

    val row =
      if (page == 0)      { 0 }
      else if (scan != 0) { 1 + scan.lowestSetBit }
      else                { table.nRows - 1 }

    val inBitW = table.inBitW
    val addr   = (row << splitW) | (pos >> inBitW).toInt
    val dx     = pos & ((BigInt(1) << inBitW) - 1)

    val polyOut = PolynomialEval.clamp(PolynomialEval.eval(addr, dx, table.coeffs),
      table.coeffs.bitW, clampNegToZero = true)

    magnitudeToBits(GenRandomGaussianFixed.roundGuard(polyOut, zSpec, guardBits), sgn == 1, zSpec)
  }

  private[rial] def magnitudeToBits(mag: BigInt, negative: Boolean, zSpec: FixedSpec): BigInt = {
    val z = if (negative) { -mag } else { mag }
    z & ((BigInt(1) << zSpec.W) - 1)
  }

  /** The signed value of the two's complement bits `bits`, in units of
   *  `2^-fracW`. */
  def toSigned(bits: BigInt, zSpec: FixedSpec): BigInt = {
    if (bits.testBit(zSpec.W - 1)) { bits - (BigInt(1) << zSpec.W) } else { bits }
  }

  /** The value of the two's complement bits `bits`. */
  def toDouble(bits: BigInt, zSpec: FixedSpec): Double = {
    toSigned(bits, zSpec).toDouble * math.pow(2.0, -zSpec.fracW)
  }
}
