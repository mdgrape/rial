package rial.approx

import scala.math._
import scala.collection.mutable.ArraySeq
import java.lang.Math.scalb

/**
 * @param nOrder   polynomial order
 * @param coefBits bit widths required to store coefficients
 * @param calcBits bit widths required to calculate polynomial
 */
case class PolynomialBitWidth(
  val nOrder: Int,
  val coefBits: Seq[Int],
  val calcBits: Seq[Int],
) {
  assert(coefBits.length == nOrder+1)
  assert(calcBits.length == nOrder+1)

  def merge(other: PolynomialBitWidth): PolynomialBitWidth = {
    require(nOrder == other.nOrder)
    PolynomialBitWidth(
      nOrder,
      coefBits.zip(other.coefBits).map { case (lhs, rhs) => max(lhs, rhs) },
      calcBits.zip(other.calcBits).map { case (lhs, rhs) => max(lhs, rhs) },
    )
  }
}

class PolynomialGenerator(
  f: Double => Double,
  val cfg: ApproxFuncConfig,
) {
  require(cfg.nOrder > 0, "if nOrder == 0, then use f(x) value itself.")

  private val funcTable = new FuncTable(f, cfg)

  private def log2(x: Double): Double = { log(x) / log(2.0) }

  // ----------------------------------------------------------------------
  // minimum required (possible max) width of coeffs for each order
  // to realize cfg.outBitW precision

  private val minCoeffs = (0 to cfg.nOrder).map(i => funcTable.intervals.map(_.coeffs(i)).min)
  private val maxCoeffs = (0 to cfg.nOrder).map(i => funcTable.intervals.map(_.coeffs(i)).max)

  private val maxCoeffMagnitude = (0 to cfg.nOrder).map(i => {
    // println(f"minCoeffs(${i}) = ${minCoeffs(i)}%10.6f, maxCoeffs(${i}) = ${maxCoeffs(i)}%10.6f")
    scala.math.max( abs(minCoeffs(i)), abs(maxCoeffs(i)) )
  })

  private val minimumRequiredCoeffBits = maxCoeffMagnitude.map(x => {
    val coeffInt  = round(scalb(x, cfg.outBitW))
    // println(f"maxCoeffMag = ${x}%10.6f, log2(coeffInt)=${ceil(log2(coeffInt)).toInt}")
    ceil(log2(coeffInt)).toInt + 1 // +1 for SInt
  })

  // ----------------------------------------------------------------------
  // minimum required (possible max) width of calc result for each order
  // to realize cfg.outBitW precision

  private val minCalc = (0 to cfg.nOrder).map(order => {
    funcTable.intervals.map(iv => {
      FindMinMax.findMin( (x: Double) => iv.eval(x, order), 0.0, 1.0 )
    }).min
  })

  private val maxCalc = (0 to cfg.nOrder).map(order => {
    funcTable.intervals.map(iv => {
      FindMinMax.findMax( (x: Double) => iv.eval(x, order), 0.0, 1.0 )
    }).max
  })

  private val maxCalcMagnitude = (0 to cfg.nOrder).map(i => {
    scala.math.max( abs(minCalc(i)), abs(maxCalc(i)) )
  })

  private val minimumRequiredCalcBits = maxCalcMagnitude.map(x => {
    val coeffInt  = round(scalb(x, cfg.outBitW))
    ceil(log2(coeffInt)).toInt + 1 // +1 for SInt
  })

  // --------------------------------------------------------------------------

  def getBitWidth(): PolynomialBitWidth = {
    PolynomialBitWidth(
      cfg.nOrder,
      minimumRequiredCoeffBits,
      minimumRequiredCalcBits,
    )
  }

  def generateTable(
    bitWidth: PolynomialBitWidth
  ): Seq[ Seq[Long /* ; cfg.nOrder+1 */] /* ; cfg.nInterval */ ] = {
    funcTable.intervals.map(iv => {
      assert(iv.coeffs.length == bitWidth.coefBits.length)
      iv.coeffs.zip(bitWidth.coefBits).map { case (coeff, bitW) => {
        // println(f"coeff=${coeff}, outBitW=${cfg.outBitW}, bitW=${bitW}")
        val bits = round(scalb(coeff, cfg.outBitW))
        assert(-(1L << (bitW-1)) <= bits && bits < (1L << (bitW-1)))
        bits
      }}
    })
  }

  def packTable(
    bitWidth: PolynomialBitWidth,
    coeffs: Seq[ Seq[Long /* ; cfg.nOrder+1 */] /* ; cfg.nInterval */ ]
  ): Seq[BigInt] = {
    coeffs.map(cs => {
      var lsb    = 0
      var result = BigInt(0)
      assert(cs.length() == bitWidth.coeffBits.length())
      for(i <- 0 until cs.length()) {
        assert( abs(cs(i)) <= (1L << (bitWidth.coeffBits(i)-1)) )
        result += (BigInt(cs(i)) << lsb)
        lsb += bitWidth.coeffBits(i)
      }
      result
    })
  }

  private def evalImpl(
    x: Long,
    bitWidth: PolynomialBitWidth,
    tableElem: Seq[Long],
  ): BigInt = {

    var calcW = bitWidth.calcBits(tableElem.length-1)
    var sum   = BigInt(tableElem (tableElem.length-1))

    // println(f"i=${tableElem.length-1}: inBitW=${cfg.inBitW}, coefW=${bitWidth.coefBits(tableElem.length-1)}, calcW=${calcW}, x=${x}, sum=${sum}, coef=${tableElem(tableElem.length-1)}")
    for(i <- (tableElem.length-2) to 0 by -1) {
      val xW   = min(calcW, cfg.inBitW) // minimum required bitW
      val prod = BigInt(x >> (cfg.inBitW - xW)) * sum

      //                                 -16                            15
      // [-1, 1] -> signed fixedpoint [10000,..11111, 00000, 00001,..01111]
      //
      // since x is signed and in [-1,1), the position of decimal point is
      // slightly different from its width.
      // if x is unsigned and in [0,1), decimal point is the same as its width.
      // but if x is signed and in [-1,1), decimal point locates width-1.
      sum   = (prod >> (xW-1)) + tableElem(i)
      calcW = bitWidth.calcBits(i)
      // println(f"i=${i}: inBitW=${cfg.inBitW}, coefW=${bitWidth.coefBits(i)}, calcW=${calcW}, x=${x}, prod=${prod>>xW}, sum=${sum}, coef=${tableElem(i)}")
    }
    sum
  }

  def eval(
    x: Double,
    bitWidth: PolynomialBitWidth,
    table: Seq[Seq[Long]],
  ): Double = {
    val idx = floor((x - cfg.xMin) * cfg.rIntervalW).toInt
    val xiv = (x - cfg.xMin - idx * cfg.intervalW) * cfg.rIntervalW * 2.0 - 1.0
    assert(-1.0 <= xiv && xiv <= 1.0)

    val xl0 = round(scalb(xiv, cfg.inBitW-1))
    val xlmin = -(1L << (cfg.inBitW - 1))
    val xlmax =  (1L << (cfg.inBitW - 1)) - 1L
    val xl = max(min(xl0, xlmax), xlmin)

    // println(f"xiv=${xiv}%8g, xl0=${xl0}, clamp(${xlmin}, ${xlmax})=${xl}, xiv-xl=${ abs(xiv - scalb(xl, -cfg.inBitW+1)) }%8g")

    val z = evalImpl(xl, bitWidth, table(idx))
    scalb(z.toDouble, -cfg.outBitW)
  }
}
