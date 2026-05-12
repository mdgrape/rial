package rial.approx

import scala.math._
import scala.collection.mutable.ArraySeq
import java.lang.Math.scalb

import chisel3.util.log2Ceil

case class PolynomialBitWidth(
  val nOrder:    Int,
  val nInterval: Int,
  val inBitW:    Int, // width of x in each interval polynomial.
  val outPrec:   Int, // (1 << outPrec) means 1.0
  val coefBits:  Seq[Int],
  val calcBits:  Seq[Int],
) {
  assert(coefBits.length == nOrder+1)
  assert(calcBits.length == nOrder+1)

  val addrW = log2Ceil(nInterval)
  val inPrec = inBitW-1 // input is [-1.0, 1.0). (1 << inPrec) == 1.0

  def merge(other: PolynomialBitWidth): PolynomialBitWidth = {
    require(nOrder == other.nOrder)
    PolynomialBitWidth(
      nOrder,
      nInterval,
      inBitW,
      outPrec,
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
  f: Double => Double,
  val cfg: ApproxFuncConfig,
) {
  require(cfg.nOrder > 0, "if nOrder == 0, then use f(x) value itself.")

  private val funcTable = new FuncTable(f, cfg)

  // ----------------------------------------------------------------------
  // range checks

  private val (xAtMinF, minF) = FindMinMax.findMin(f, 0.0, 1.0)
  private val (xAtMaxF, maxF) = FindMinMax.findMax(f, 0.0, 1.0)
  if(minF < 0.5) {
    println(f"WARNING: f(x=${xAtMinF}) = ${minF} < 0.5 may lose precision.")
  }
  assert(maxF <= 1.0 f"WARNING: f(x=${xAtMaxF}) = ${maxF} > 1.0 may cause overflow.")

  // ----------------------------------------------------------------------
  // minimum required (possible max) width of coeffs for each order
  // to realize cfg.outPrec precision

  private val minCoeffs = (0 to cfg.nOrder).map(i => funcTable.intervals.map(_.coeffs(i)).min)
  private val maxCoeffs = (0 to cfg.nOrder).map(i => funcTable.intervals.map(_.coeffs(i)).max)

  private val maxCoeffMagnitude = (0 to cfg.nOrder).map(i => {
    // println(f"minCoeffs(${i}) = ${minCoeffs(i)}%10.6f, maxCoeffs(${i}) = ${maxCoeffs(i)}%10.6f")
    scala.math.max( abs(minCoeffs(i)), abs(maxCoeffs(i)) )
  })

  private val minimumRequiredCoeffBits = maxCoeffMagnitude.map(x => {
    val coeffInt  = round(scalb(x, cfg.outPrec))
    // println(f"maxCoeffMag = ${x}%10.6f, log2(coeffInt)=${log2Ceil(coeffInt).toInt}")
    log2Ceil(coeffInt).toInt + 1 // +1 for SInt
  })

  // ----------------------------------------------------------------------
  // minimum required (possible max) width of calc result for each order
  // to realize cfg.outPrec precision

  private val minCalc = (0 to cfg.nOrder).map(order => {
    funcTable.intervals.map(iv => {
      FindMinMax.findMin( (x: Double) => iv.eval(x, order), 0.0, 1.0 )._1
    }).min
  })

  private val maxCalc = (0 to cfg.nOrder).map(order => {
    funcTable.intervals.map(iv => {
      FindMinMax.findMax( (x: Double) => iv.eval(x, order), 0.0, 1.0 )._1
    }).max
  })

  private val maxCalcMagnitude = (0 to cfg.nOrder).map(i => {
    scala.math.max( abs(minCalc(i)), abs(maxCalc(i)) )
  })

  private val minimumRequiredCalcBits = maxCalcMagnitude.map(x => {
    val coeffInt  = round(scalb(x, cfg.outPrec))
    log2Ceil(coeffInt).toInt + 1 // +1 for SInt
  })

  // --------------------------------------------------------------------------

  def getBitWidth(): PolynomialBitWidth = {
    PolynomialBitWidth(
      cfg.nOrder,
      cfg.nInterval,
      cfg.inBitW,
      cfg.outPrec,
      minimumRequiredCoeffBits,
      minimumRequiredCalcBits,
    )
  }

  def generateTable(
    bitWidth: PolynomialBitWidth
  ): PolynomialCoefficients = {
    PolynomialCoefficients(bitWidth,
      funcTable.intervals.map(iv => {
        assert(iv.coeffs.length == bitWidth.coefBits.length)
        iv.coeffs.zip(bitWidth.coefBits).map { case (coeff, bitW) => {
          // println(f"coeff=${coeff}, outPrec=${cfg.outPrec}, bitW=${bitW}")
          val bits = round(scalb(coeff, cfg.outPrec))
          assert(-(1L << (bitW-1)) <= bits && bits < (1L << (bitW-1)))
          bits
        }}
      })
    )
  }

  def packTable(
    coeffs: PolynomialCoefficients
  ): Seq[BigInt] = {
    coeffs.coeffs.map(cs => {
      var lsb    = 0
      var result = BigInt(0)
      assert(cs.length() == coeffs.bitW.coeffBits.length())
      for(i <- 0 until cs.length()) {
        assert( abs(cs(i)) <= (1L << (coeffs.bitW.coeffBits(i)-1)) )
        result += (BigInt(cs(i)) << lsb)
        lsb += coeffs.bitWidth.coeffBits(i)
      }
      result
    })
  }
}


