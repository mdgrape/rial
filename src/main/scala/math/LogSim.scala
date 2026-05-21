//% @file LogSim.scala
//
// Simulator for log(x).
//
// Three regions, three tables:
// - x < 0.5 or 2 <= x: log(x) = ex * log(2) + log(1.man): uses  log2(1+t)    table
// -     0.5 <= x < 1 : log(x) = -log(x)/(1-x) * -(1-x)  : uses -log(y)/(1-y) table
// -     1.0 <= x < 2 : log(x) =  log(x)/(x-1) *  (x-1)  : uses  log(y)/(y-1) table
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math._

import spire.math.SafeLong
import spire.math.Real
import spire.implicits._

import rial.approx._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding._

class LogTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  private val log2 = (a: Double) => log(a) / log(2.0)

  // path 0 (x < 0.5 or x >= 2): log2(1+t) for t in [0, 1)
  private val funcN = (t: Double) => { log2(1.0 + t) }

  // path 1 (0.5 <= x < 1.0): -log((1+t)/2) / (1 - (1+t)/2) - 1 for t in [0, 1)
  private val funcG = (t: Double) => {
    if (t == 1.0) {
      0.0
    } else {
      val x = (1.0 + t) * 0.5
      -log(x) / (1.0 - x) - 1.0
    }
  }
  // path 2 (1.0 <= x < 2.0): log(1+t)/t  (=> 1 at t=0) for t in [0, 1)
  private val funcP = (t: Double) => {
    if (t == 0.0) {
      1.0 
    } else {
      log(1.0 + t) / t
    }
  }

  private val nIntervals = 1 << addrW

  private val ftN = FuncTable(funcN, nOrder, nIntervals, 0.0, 1.0)
  private val ftG = FuncTable(funcG, nOrder, nIntervals, 0.0, 1.0)
  private val ftP = FuncTable(funcP, nOrder, nIntervals, 0.0, 1.0)

  private val genN = new PolynomialGenerator(ftN, inBitW, outBitW)
  private val genG = new PolynomialGenerator(ftG, inBitW, outBitW)
  private val genP = new PolynomialGenerator(ftP, inBitW, outBitW)

  def getBitWidth(): PolynomialBitWidth =
    genN.getBitWidth().merge(genG.getBitWidth()).merge(genP.getBitWidth())

  def genNormalCoeffs  (bitW: PolynomialBitWidth): PolynomialCoefficients = {genN.generateTable(bitW)}
  def genSmallNegCoeffs(bitW: PolynomialBitWidth): PolynomialCoefficients = {genG.generateTable(bitW)}
  def genSmallPosCoeffs(bitW: PolynomialBitWidth): PolynomialCoefficients = {genP.generateTable(bitW)}

  /** Concatenated form: [normal] ++ [smallNeg] ++ [smallPos] ++ [zeros].
   *  Consumed by the HW (LogTableCoeff) which slices it back into three
   *  independent ROMs. */
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    val cN = genNormalCoeffs(bitW)
    val cG = genSmallNegCoeffs(bitW)
    val cP = genSmallPosCoeffs(bitW)
    val zeros = Seq.fill(nIntervals)(Seq.fill(nOrder + 1)(0L))
    PolynomialCoefficients(bitW, cN.coeffs ++ cG.coeffs ++ cP.coeffs ++ zeros)
  }
}

object LogSimGeneric {

  /** Compute log(x).
   *
   *  @param tNormal   polynomial coefficients for the log2(1+t) table
   *  @param tSmallPos polynomial coefficients for the log(1+t)/t   table
   *  @param tSmallNeg polynomial coefficients for the -log((1+t)/2)/(1-(1+t)/2) - 1 table
   *  @param x         input
   */
  def apply(
    tNormal:   PolynomialCoefficients,
    tSmallPos: PolynomialCoefficients,
    tSmallNeg: PolynomialCoefficients,
    x:         RealGeneric,
  ): RealGeneric = {

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val xnan  = x.isNaN
    val xinf  = x.isInfinite
    val xzero = x.isZero
    val xneg  = x.sgn == 1

    if (xnan)         { return RealGeneric.nan(spec)         }
    if (xinf && !xneg) { return RealGeneric.inf(spec, sgn = 0) }
    if (xzero)        { return RealGeneric.inf(spec, sgn = 1) }
    if (xneg)         { return RealGeneric.nan(spec)         }
    if (x.man == SafeLong(0) && x.ex == exBias) {
      return RealGeneric.zero(spec)
    }

    val zsgn = if (x.ex < exBias) 1 else 0

    // All tables produce a fracW-wide unsigned value scaled by 2^outPrec.
    val addrW   = tNormal.bitW.addrW
    val inBitW  = tNormal.bitW.inBitW
    val outPrec = tNormal.bitW.outPrec
    val fracW   = outPrec
    assert(tSmallPos.bitW.addrW == addrW && tSmallNeg.bitW.addrW == addrW)
    assert(tSmallPos.bitW.outPrec == outPrec && tSmallNeg.bitW.outPrec == outPrec)

    val xexNobias = x.ex - exBias
    val zint0     = if (xexNobias >= 0) SafeLong(xexNobias) else SafeLong(-xexNobias - 1)

    val xW = manW - addrW

    val (ymanW1, zex0, cmanW1) = if (x.ex == exBias - 1) {
      // -----------------------------------------------------------------------
      // 0.5 <= x < 1.0  ->  -log(x)/(1-x) * -(1-x)

      val addr = slice(xW, addrW, x.man).toInt
      val xiv  = slice(0,  xW,    x.man).toBigInt

      val res = PolynomialEval.eval(addr, xiv, tSmallNeg)
      val ymanW1 = if (res < 0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        SafeLong(1) << fracW
      } else if (res >= (SafeLong(1) << fracW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        maskSL(fracW + 1)
      } else {
        (SafeLong(1) << fracW) + SafeLong(res)
      }
      // result will be in [ln(2), 1.0).  ln(2) > 0.5  =>  yex == -1.
      val yex0 = -1
      assert(bit(fracW, ymanW1) == 1, f"ymanW1 = ${ymanW1.toLong.toBinaryString}")

      val c0        = (SafeLong(1) << (manW + 1)) - x.manW1
      val c0Shift   = manW + 1 - binaryWidthSL(c0)
      val c0Shifted = c0 << c0Shift
      assert((SafeLong(1) << manW) <= c0Shifted && c0Shifted < (SafeLong(1) << (manW + 1)))

      val zex0 = exBias + yex0 - c0Shift
      (ymanW1, zex0.toInt, c0Shifted)

    } else if (x.ex == exBias) {
      // -----------------------------------------------------------------------
      // 1.0 <= x < 2.0  ->  log(x)/(x-1) * (x-1)

      val addr = slice(xW, addrW, x.man).toInt
      val xiv  = slice(0,  xW,    x.man).toBigInt

      val res = PolynomialEval.eval(addr, xiv, tSmallPos)
      val ymanW1 = if (res < 0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        SafeLong(0)
      } else if (res >= (SafeLong(1) << fracW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        maskSL(fracW + 1)
      } else {
        SafeLong(res) << 1 // normalize: bring fracW-th bit to 1
      }
      val yex0 = 0
      assert(bit(fracW, ymanW1) == 1, f"ymanW1 = ${ymanW1.toLong.toBinaryString}")

      val c0        = x.man
      val c0Shift   = manW + 1 - binaryWidthSL(c0)
      val c0Shifted = c0 << c0Shift
      assert((SafeLong(1) << manW) <= c0Shifted && c0Shifted < (SafeLong(1) << (manW + 1)))

      val zex0 = exBias + yex0 - c0Shift - 1 // -1 because ymanW1 was shifted
      (ymanW1, zex0.toInt, c0Shifted)

    } else {
      // -----------------------------------------------------------------------
      // polynomial path (x < 0.5 or x >= 2.0): log(x) = ex*log(2) + log(1.man)

      val addr = slice(xW, addrW, x.man).toInt
      val xiv  = slice(0,  xW,    x.man).toBigInt

      val zfrac0Pos0 = PolynomialEval.eval(addr, xiv, tNormal)
      val zfrac0Pos = if (zfrac0Pos0 < 0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        SafeLong(0)
      } else if (zfrac0Pos0 >= (SafeLong(1) << fracW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        maskSL(fracW)
      } else {
        SafeLong(zfrac0Pos0)
      }

      val zfrac0 = if (xexNobias >= 0) {
        zfrac0Pos
      } else if (zfrac0Pos == SafeLong(0)) {
        maskSL(fracW)
      } else {
        (SafeLong(1) << fracW) - zfrac0Pos
      }

      val zfrac  = zfrac0 & maskSL(fracW)
      val zfull0 = (zint0 << fracW) + zfrac

      assert(0 <= zfrac && zfrac < (SafeLong(1) << fracW))
      assert(0 <= zfull0)
      val zfullW  = binaryWidthSL(zfull0)
      val zShiftW = exW + fracW - zfullW
      assert(0 <= zShiftW)

      val zShifted = zfull0 << zShiftW
      assert(bit(exW + fracW - 1, zShifted) == 1)

      val yman0  = slice(exW - 1, fracW, zShifted) // -1 for the hidden bit
      val ymanW1 = (SafeLong(1) << fracW) + yman0
      val zex0   = exBias + (exW - 1) - zShiftW - 1 // ln2.ex == -1
      val cmanW1 = Real.log(Real.two)(manW + 1)
      assert((SafeLong(1) << manW) < cmanW1 && cmanW1 < (SafeLong(1) << (manW + 1)))

      (ymanW1, zex0.toInt, cmanW1)
    }

    // -----------------------------------------------------------------------
    // multiply

    val zmanProd          = ymanW1 * cmanW1
    val zmanProdMoreThan2 = bit((fracW + 1) + (manW + 1) - 1, zmanProd).toInt
    val zmanRoundBits     = fracW + zmanProdMoreThan2
    val zmanRound         = roundBySpec(RoundSpec.roundToEven, zmanRoundBits, zmanProd)
    val zmanRoundMoreThan2 = bit(manW + 1, zmanRound).toInt
    val lnman = zmanRound >> zmanRoundMoreThan2
    val lnex  = zex0 + zmanProdMoreThan2 + zmanRoundMoreThan2

    new RealGeneric(spec, zsgn, lnex, lnman)
  }
}
