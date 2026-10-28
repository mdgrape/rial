package rial.table

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

  /** The widths that cover both: the wider coefficient and accumulator of the
   *  two at every order, and the wider of the two addresses. The two have to
   *  agree on the order, the input and the output.
   */
  def merge(other: PolynomialBitWidth): PolynomialBitWidth = {
    require(nOrder == other.nOrder)
    require(inBitW == other.inBitW)
    require(outBitW == other.outBitW)
    PolynomialBitWidth(
      nOrder,
      max(addrW, other.addrW),
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

/** The integer coefficients of the polynomial of every interval of `funcTable`,
 * and the widths they and the intermediate results of [[PolynomialEval]] need.
 *
 * [[generateTable]] refuses a table whose polynomial could reach a value the
 * evaluator does not handle correctly: a result below zero leaves the output as
 * all ones, i.e. as ~1.0, unless the evaluator clamps it to zero. It reads a
 * row at every position of `inBitW` bits, so a table whose argument carries
 * fewer bits than that is held to positions it never reaches as well.
 *
 * @param funcTable         the function and the intervals to approximate
 * @param inBitW            width of the position inside an interval; signed,
 *                          so `1 << (inBitW-1)` is 1.0
 * @param outBitW           `1 << outBitW` is 1.0, so the function must stay
 *                          below 1.0
 * @param negativeIsClamped whether the evaluator that reads this table clamps a
 *                          negative result to zero, i.e. is built with
 *                          `PolynomialEval(clampNegToZero = true)`; a table of
 *                          such an evaluator may go below zero, so nothing is
 *                          refused
 */
class PolynomialGenerator(
  private val funcTable: FuncTable,
  val inBitW:  Int, // signed, so (1 << (inBitW-1)) = 1.0 / -1.0
  val outBitW: Int, // 1 << outBitW == 1.0 (so, it must be f(x) < 1.0)
  val negativeIsClamped: Boolean = false,
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
    log2Ceil(max(coeffInt, 1L)).toInt + 1 // +1 for SInt; a zero takes one bit
  })

  // ----------------------------------------------------------------------
  // minimum required (possible max) width of the horner accum at each stage
  //
  // the accumulator of an order is the polynomial of the coefficients from
  // that order up.
  //
  //          ((c2) * x + c1) * x + c0
  // order 0: ^^^^^^^^^^^^^^^^^^^^^^^^ c0 + c1 x + c2 x^2
  // order 1:  ^^^^^^^^^^^^^           c1 + c2 x
  // order 2:   ^^                     c2
  private val calcRangeOf = (0 to nOrder).map(order => {
    funcTable.intervals.map(iv => {
      FindMinMax.polynomial(iv.coeffs.drop(order), -1.0, 1.0)
    }).toSeq
  })
  private val minCalcOf = calcRangeOf.map(_.map(_._1))
  private val minCalc = minCalcOf.map(_.min)
  private val maxCalc = calcRangeOf.map(_.map(_._2).max)

  private val maxCalcMagnitude = (0 to nOrder).map(i => {
    scala.math.max( abs(minCalc(i)), abs(maxCalc(i)) )
  })

  private val minimumRequiredCalcBits = maxCalcMagnitude.map(x => {
    val coeffInt = round(scalb(x, outBitW))
    log2Ceil(max(coeffInt, 1L)).toInt + 1 // +1 for SInt; a zero takes one bit
  })

  println(f"outBitW=${outBitW}, " +
    f"coefBits=${minimumRequiredCoeffBits.mkString("[",",","]")}, " +
    f"calcBits=${minimumRequiredCalcBits .mkString("[",",","]")}")

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

  // The output of the evaluator is unsigned, so a row that reaches a negative
  // value leaves it as all ones, i.e. as ~1.0, where the function is ~0.0.
  //
  // Rounding the coefficients and the truncation of a Horner step move a row by
  // at most this many LSBs of the output: 1.5 per coefficient (half an LSB of
  // rounding, or a whole one where a coefficient that does not fit is clamped),
  // and 2 per step, which truncates both the product and the position it
  // multiplies. Errors of a step are multiplied by |x| <= 1 afterwards, so they
  // add up. One LSB on top is a margin for the rounding of the minimum.
  private val quantError = 1.5 * (nOrder + 1) + 2.0 * nOrder + 1.0

  /** Fails unless every row of `table` stays at or above zero. Rows whose
   *  double polynomial keeps [[quantError]] of headroom cannot go below it;
   *  the rest are evaluated at every position of the interval. */
  private def requireNonNegative(table: PolynomialCoefficients): Unit = {
    if (negativeIsClamped) { return }

    // a table of order zero is a plain lookup: it reads no position at all
    val nPos = if (nOrder == 0) { 1 } else { 1 << inBitW }

    funcTable.intervals.zipWithIndex.foreach { case (iv, i) =>
      if (scalb(minCalcOf(0)(i), outBitW) - quantError < 0.0) {
        val minVal = (0 until nPos).view.map(dx => {
          PolynomialEval.eval(i, BigInt(dx), table)
        }).min

        if (minVal < 0) {
          throw new IllegalArgumentException(
            f"interval ${i}, x in [${iv.xMin}, ${iv.xMax}), reaches " +
            f"${minVal} below zero")
        }
      }
    }
  }

  def generateTable(
    bitW: PolynomialBitWidth
  ): PolynomialCoefficients = {
    val table = PolynomialCoefficients(
      bitW,
      funcTable.intervals.map(iv => {
        assert(iv.coeffs.length == bitW.coefBits.length)
        iv.coeffs.zip(bitW.coefBits).map { case (coeff, coefBitW) => {
          val bits = round(scalb(coeff, outBitW))

          val (minVal, maxVal) = if(nOrder == 0) {
            // there is no calculation, so coefficient itself is the output.
            // since output is mantissa, it must be in [0, 1) (excl hidden bit).
            val minV = 0L
            val maxV = (1L << (coefBitW-1)) - 1
            (minV, maxV)

          } else {
            val minV = -(1L << (coefBitW-1))
            val maxV =  (1L << (coefBitW-1)) - 1
            (minV, maxV)
          }

          // overflow by rounding (e.g. 0.999 -> 1.000) is okay.
          // overflow by the function itself (e.g. f(x)=1.001) is not okay.
          assert( bits >= minVal - 1, f"coeff bit(${bits}) must be > than 0b100..0(${minVal})" )
          assert( bits <= maxVal + 1, f"coeff bit(${bits}) must be < than 0b011..1(${maxVal})" )

          bits.max(minVal).min(maxVal)
        }}
      }).toSeq
    )
    requireNonNegative(table)
    table
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

