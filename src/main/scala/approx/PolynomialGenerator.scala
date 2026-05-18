package rial.approx

import scala.math._
import java.lang.Math.scalb
import chisel3.util.log2Ceil

import rial.util.ScalaUtil._

case class PolynomialBitWidth(
  val nOrder:   Int,
  val addrW:    Int,
  val inBitW:   Int, // width of x in each interval polynomial.
  val outBitW:  Int, // (1 << outBitW) means 1.0 (so it must be f(x) < 1.0)
  val coefBits: Seq[Int],
  val calcBits: Seq[Int],
) {
  assert(coefBits.length == nOrder+1)
  assert(calcBits.length == nOrder+1)

  val inPrec  = inBitW-1 // input is [-1.0, 1.0). (1 << inPrec ) == 1.0
  val outPrec = outBitW  // output is [0.0, 1.0). (1 << outPrec) == 1.0

  def merge(other: PolynomialBitWidth): PolynomialBitWidth = {
    require(nOrder == other.nOrder)
    require(addrW == other.addrW)
    require(inBitW == other.inBitW)
    require(outBitW == other.outBitW)
    PolynomialBitWidth(
      nOrder,
      addrW,
      inBitW,
      outBitW,
      coefBits.zip(other.coefBits).map { case (lhs, rhs) => max(lhs, rhs) },
      calcBits.zip(other.calcBits).map { case (lhs, rhs) => max(lhs, rhs) },
    )
  }
}

case class PolynomialCoefficients(
  val bitW:   PolynomialBitWidth,
  val coeffs: Seq[ Seq[Long /* ; nOrder+1 */] /* ; nInterval */ ]
) {
}

class PolynomialGenerator(
  private val funcTable: FuncTable,
  val inBitW:  Int, // signed, so (1 << (inBitW-1)) = 1.0 / -1.0
  val outBitW: Int, // 1 << outBitW == 1.0 (so, it must be f(x) < 1.0)
) {

  val nOrder = funcTable.nOrder
  val addrW  = log2Ceil(funcTable.intervals.length)

  // ----------------------------------------------------------------------
  // minimum required (possible max) width of coeffs for each order
  // to realize outBitW precision

  private val minCoeffs = (0 to nOrder).map(i => funcTable.intervals.map(_.coeffs(i)).min)
  private val maxCoeffs = (0 to nOrder).map(i => funcTable.intervals.map(_.coeffs(i)).max)

  private val maxCoeffMagnitude = (0 to nOrder).map(i => {
    // println(f"minCoeffs(${i}) = ${minCoeffs(i)}%10.6f, maxCoeffs(${i}) = ${maxCoeffs(i)}%10.6f")
    scala.math.max( abs(minCoeffs(i)), abs(maxCoeffs(i)) )
  })

  private val minimumRequiredCoeffBits = maxCoeffMagnitude.map(x => {
    val coeffInt  = round(scalb(x, outBitW))
    // println(f"maxCoeffMag = ${x}%10.6f, log2(coeffInt)=${log2Ceil(coeffInt).toInt}")
    log2Ceil(coeffInt).toInt + 1 // +1 for SInt
  })

  // ----------------------------------------------------------------------
  // minimum required (possible max) width of calc result for each order
  // to realize outBitW precision

  private val minCalc = (0 to nOrder).map(order => {
    funcTable.intervals.map(iv => {
      FindMinMax.findMin( (x: Double) => iv.eval(x, order), 0.0, 1.0 )._1
    }).min
  })

  private val maxCalc = (0 to nOrder).map(order => {
    funcTable.intervals.map(iv => {
      FindMinMax.findMax( (x: Double) => iv.eval(x, order), 0.0, 1.0 )._1
    }).max
  })

  private val maxCalcMagnitude = (0 to nOrder).map(i => {
    scala.math.max( abs(minCalc(i)), abs(maxCalc(i)) )
  })

  private val minimumRequiredCalcBits = maxCalcMagnitude.map(x => {
    val coeffInt = round(scalb(x, outBitW))
    log2Ceil(coeffInt).toInt + 1 // +1 for SInt
  })

  // --------------------------------------------------------------------------

  def getBitWidth(): PolynomialBitWidth = {
    PolynomialBitWidth(
      nOrder,
      addrW,
      inBitW,
      outBitW,
      minimumRequiredCoeffBits,
      minimumRequiredCalcBits,
    )
  }

  def generateTable(
    bitW: PolynomialBitWidth
  ): PolynomialCoefficients = {
    PolynomialCoefficients(
      bitW,
      funcTable.intervals.map(iv => {
        assert(iv.coeffs.length == bitW.coefBits.length)
        iv.coeffs.zip(bitW.coefBits).map { case (coeff, bitW) => {
          val bits = round(scalb(coeff, outBitW))
          assert( -(1L << (coeffs.bitW.coefBits(i)-1)) - 1 < bits, f"coeff bit(${bits}) must be > than 0b100..0(${-(1L << (coeffs.bitW.coefBits(i)-1))})" )
          assert( bits < (1L << (coeffs.bitW.coefBits(i)-1)) + 1,  f"coeff bit(${bits}) must be < than 0b011..1(${ (1L << (coeffs.bitW.coefBits(i)-1))})" )

          val minVal = -(1L << (bitW-1))
          val maxVal = (1L << (bitW-1)) - 1
          // Clip to the valid range for signed integer
          bits.max(minVal).min(maxVal)
        }}
      })
    )
  }
}

object PolynomialGenerator {
  def packTable(
    coeffs: PolynomialCoefficients
  ): Seq[BigInt] = {
    coeffs.coeffs.map(cs => {
      var lsb    = 0
      var result = BigInt(0)
      assert(cs.length == coeffs.bitW.coefBits.length)
      for(i <- 0 until cs.length) {
        assert( abs(cs(i)) <= (1L << (coeffs.bitW.coefBits(i)-1)) )
        val c = BigInt(cs(i)) & maskBI(coeffs.bitW.coefBits(i))
        result |= (c << lsb)
        lsb = lsb + coeffs.bitW.coefBits(i)
      }
      result
    })
  }
}

