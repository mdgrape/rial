//% @file ScaleMixtureGaussianSim.scala
//
// Simulator for d/dx log p(x), where p is a 2-Gaussian mixture with stddev
// `sgmA` (large) and `sgmB` (small, sgmB/sgmA << 1).
//
// Math:
//   1/sgm'^2 = (sgmA^2 - sgmB^2) / (sgmA^2 sgmB^2)
//   g(x)     = exp(-x^2 / (2 sgm'^2))
//   f(x)     = -x/sgmA^2 [ (sgmA^2/sgm'^2) / (sgmB/(sgmA g(x)) + 1) + 1 ]
//                          ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//                          this term is approximated by the polynomial.
//
// The table approximates  (sgmA^2/sgm'^2) / (sgmB/(sgmA g(x))+1)  scaled by
// 1/tableMaxValue so the result fits in [0, 1).  The domain is [0, domainMax)
// (mapped from [0, 1) at the polynomial level), beyond which g(x) -> 0 and the
// term becomes negligible.
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math._

import spire.math.SafeLong
import spire.implicits._

import rial.approx._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding

class ScaleMixtureGaussianTable(
  val sgmA:    Double,
  val sgmB:    Double,
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  assert(sgmB / sgmA < 1.0, "sgmB / sgmA << 1")

  // The output of this table is f(x) / tableMaxValue, where
  //   f(x) = (sgmA^2/sgm'^2) / (sgmB/(sgmA g(x)) + 1),
  // scaled so the range fits in [0, 1).
  private val sgmBoverA   = sgmB / sgmA
  private val sgmA2       = sgmA * sgmA
  private val sgmB2       = sgmB * sgmB
  private val sgmPrime2   = (sgmA2 * sgmB2) / (sgmA2 - sgmB2)
  private val sgmA2overP2 = sgmA2 / sgmPrime2

  // Domain: input x in [0, domainMax), where domainMax = 2^domainDigit so that
  // f(x) at x = domainMax is small (the table value rounds to 0).
  val domainDigit  = ScaleMixtureGaussianTable.tableDomainDigit(outBitW, sgmA, sgmB)
  val maxBitDigit  = ScaleMixtureGaussianTable.tableMaxBitDigit(sgmA, sgmB)
  val domainMax    = pow(2.0, domainDigit)

  private val func = (x01: Double) => {
    val x = x01 * domainMax
    val g = exp(-x * x / (2 * sgmPrime2))
    val v = sgmA2overP2 / (sgmBoverA / g + 1.0)
    v * pow(2.0, -maxBitDigit) // scale into [0, 1)
  }

  private val nIntervals = 1 << addrW
  private val ft  = FuncTable(func, nOrder, nIntervals, 0.0, 1.0)
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW)

  def getBitWidth(): PolynomialBitWidth = { gen.getBitWidth() }
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    gen.generateTable(bitW)
  }
}

object ScaleMixtureGaussianTable {

  /** log2 of x's value above which the table value is essentially 0 (rounds
   *  away within `precW` bits of precision). Match the old math/ port.
   */
  def tableDomainDigit(precW: Int, sgmA: Double, sgmB: Double): Int = {
    val log2 = (x: Double) => log(x) / log(2.0)

    val sgmA2     = sgmA * sgmA
    val sgmB2     = sgmB * sgmB
    val sgmPrime2 = (sgmA2 * sgmB2) / (sgmA2 - sgmB2)

    val boundary = sqrt(2.0 * sgmPrime2 * log(
      sgmA / sgmB * (pow(2.0, precW) * sgmA2 / sgmPrime2 - 1.0)
    ))

    ceil(log2(boundary)).toInt
  }

  /** The maximum value of the un-scaled table function (taken at x=0). */
  def tableMaxValue(sgmA: Double, sgmB: Double): Double = {
    (sgmA / sgmB) * ((sgmA - sgmB) / sgmB)
  }

  /** ceil(log2(tableMaxValue)). */
  def tableMaxBitDigit(sgmA: Double, sgmB: Double): Int = {
    val log2 = (x: Double) => log(x) / log(2.0)
    val d = ceil(log2(tableMaxValue(sgmA, sgmB))).toInt
    assert(d > 0)
    d
  }
}

/** An interface for software simulator.
 */
object ScaleMixtureGaussianSimGeneric {

  def apply(
    t:    PolynomialCoefficients,
    x:    RealGeneric,
    sgmA: Double = exp(-1.0),
    sgmB: Double = exp(-6.0),
  ): RealGeneric = {

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    // smg is the *derivative* of log p; the result has the opposite sign of x.
    val zSgn = if (x.sgn == 0) 1 else 0

    if (x.isNaN)      { return RealGeneric.nan (spec)       }
    if (x.isInfinite) { return RealGeneric.inf (spec, zSgn) }
    if (x.isZero)     { return RealGeneric.zero(spec)       }

    val addrW     = t.bitW.addrW
    val outPrec   = t.bitW.outPrec
    val bp        = outPrec
    val extraBits = outPrec - manW
    val nOrder    = t.bitW.nOrder

    // ---- scale x into the table domain ------------------------------------
    val xTableMaxDigit = ScaleMixtureGaussianTable.tableDomainDigit(outPrec, sgmA, sgmB)
    val xTableShift0   = xTableMaxDigit - (x.ex - exBias)
    val xTableShift    = abs(xTableShift0)
    val xTableMan      = x.manW1 >> xTableShift

    // ---- polynomial table query ------------------------------------------
    val xW   = manW - addrW
    val addr = slice(xW, addrW, xTableMan).toInt
    val xiv  = slice(0,  xW,    xTableMan).toBigInt
    val zTable0 = PolynomialEval.eval(addr, xiv, t)

    // if |x| >= domainMax, g(x) -> 0 and the table term vanishes.
    val zTable: Long = if (xTableShift0 > 0) zTable0.toLong else 0L

    // ---- z0 = (table scaled back) + 1 -------------------------------------
    val zTableScaleDigit = ScaleMixtureGaussianTable.tableMaxBitDigit(sgmA, sgmB)
    val zTableScaled     = SafeLong(zTable) << zTableScaleDigit
    val zman0            = zTableScaled + (SafeLong(1) << bp)
    val zman0W           = binaryWidthSL(zman0)
    val zman0Shift       = zman0W - (manW + 1 + extraBits)
    val zman0Rounded     = if (zman0Shift > 0) {
      (zman0 >> zman0Shift) + bit(zman0Shift - 1, zman0)
    } else {
      zman0
    }
    val zman0MoreThan2AfterRound = bit(manW + 1 + extraBits + 1, zman0Rounded)

    val z0ManW1 = if (zman0MoreThan2AfterRound == 1) { zman0Rounded >> 1 } else { zman0Rounded }
    val z0ExNoBias = zman0W - 1 - bp + zman0MoreThan2AfterRound

    // ---- x / sgmA^2 (emulated pre-mult) -----------------------------------
    val rsgmA2 = new RealGeneric(spec, 1.0 / (sgmA * sgmA))
    val xOverSgmA2Prod          = x.manW1 * rsgmA2.manW1
    val xOverSgmA2ProdMoreThan2 = bit((1 + manW) + (1 + manW) - 1, xOverSgmA2Prod)
    val xOverSgmA2ProdRoundedW1 = Rounding.roundToEven(
      manW + xOverSgmA2ProdMoreThan2, xOverSgmA2Prod)
    val xOverSgmA2ProdMoreThan2AfterRounding = bit(manW + 1, xOverSgmA2ProdRoundedW1)
    val xOverSgmA2ProdEx = x.ex + rsgmA2.ex - rsgmA2.spec.exBias +
      xOverSgmA2ProdMoreThan2 + xOverSgmA2ProdMoreThan2AfterRounding

    // ---- final multiply: (x/sgmA^2) * z0 ----------------------------------
    val zProd                    = xOverSgmA2ProdRoundedW1 * z0ManW1
    val zProdMoreThan2           = bit((1 + manW) + (1 + manW + extraBits) - 1, zProd)
    val zProdRoundedW1           = Rounding.roundToEven(
      (manW + extraBits) + zProdMoreThan2, zProd)
    val zProdMoreThan2AfterRound = bit(manW + 1, zProdRoundedW1)
    val zProdEx = xOverSgmA2ProdEx + z0ExNoBias +
      zProdMoreThan2 + zProdMoreThan2AfterRound

    val zProdMan = slice(0, manW, zProdRoundedW1)

    if (zProdEx >= maskI(exW)) {
      new RealGeneric(spec, zSgn, maskI(exW), 0) // overflow to Inf
    } else {
      new RealGeneric(spec, zSgn, zProdEx.toInt, zProdMan)
    }
  }
}
