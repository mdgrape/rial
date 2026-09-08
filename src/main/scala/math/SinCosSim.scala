//% @file SinCosSim.scala
//
// Simulators for sin/cos function (HW pre/post multipliers are emulated
// at full precision; the precision of 1/pi used by the emulated pre-mult
// is exposed as an explicit parameter so the hardware can match it.)
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._

import spire.math.SafeLong
import spire.math.Real
import spire.implicits._

import rial.table._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.RoundSpec
import rial.arith.Rounding._

// The table approximates  z(y0) = sin(Pi * y0 / 2) / y0    for y0 in (0, 1)
//                          z(0) = Pi/2                     (by limit)
// z is in [1, Pi/2), so the stored value is `z - 1` in [0, Pi/2 - 1).
// (y0 at runtime corresponds to 2y, where y in [0, 1/2) is the argument-
//  reduced |x|/pi after quadrant folding.)
class SinCosTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  private val func = (y0: Double) => {
    val z = if (y0 == 0.0) Pi / 2.0 else sin(Pi * y0 / 2.0) / y0
    z - 1.0
  }

  private val nIntervals = 1 << addrW
  private val ft  = FuncTable(func, nOrder, nIntervals, 0.0, 1.0)
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW)

  def getBitWidth(): PolynomialBitWidth = {
    gen.getBitWidth()
  }
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    gen.generateTable(bitW)
  }
}

object SinCosTable {
  def fp32(): SinCosTable = {
    val addrW     = 8
    val extraBits = 3
    new SinCosTable(
      nOrder  = 2,
      addrW   = addrW,
      inBitW  = (RealSpec.Float32Spec.manW - addrW),
      outBitW = (RealSpec.Float32Spec.manW + extraBits),
    )
  }
  def bf16(): SinCosTable = {
    val addrW     = RealSpec.BFloat16Spec.manW
    val extraBits = 0
    new SinCosTable(
      nOrder  = 0,
      addrW   = addrW,
      inBitW  = (RealSpec.BFloat16Spec.manW - addrW),
      outBitW = (RealSpec.BFloat16Spec.manW + extraBits),
    )
  }

  private lazy val fp32Table = fp32().generateTable(fp32().getBitWidth())

  /** sin simulation with the default fp32 configuration */
  def sinDefaultFP32(x: Float): Float = {
    SinCosSimGeneric(true,  fp32Table,
      RealGeneric.fromFloat(RealSpec.Float32Spec, x)).toFloat
  }
  /** cos simulation with the default fp32 configuration */
  def cosDefaultFP32(x: Float): Float = {
    SinCosSimGeneric(false, fp32Table,
      RealGeneric.fromFloat(RealSpec.Float32Spec, x)).toFloat
  }
}

/** Software simulator for the hardware sincos pipeline.
 *
 *  The hardware performs an `x * (1/pi)` multiplication in its pre-process
 *  stage; this simulator emulates that multiplication at the same precision
 *  as the hardware. The width of the 1/pi pad (`oneOverPiPad`) is exposed
 *  explicitly so that the eventual hardware port can use the same value.
 */
object SinCosSimGeneric {

  /** Default precision-pad bits for the 1/pi constant.
   *
   *  The 1/pi constant width is `manW + 2 + oneOverPiPad`. The default
   *  `max(manW+1, 12)` matches the original math/ port, where 12 was the
   *  minimum width needed to keep BFloat16 results within 2 LSB of an
   *  FP64-reference computation.
   */
  def defaultOneOverPiPad(spec: RealSpec): Int = max(spec.manW + 1, 12)

  /** Compute sin(x) or cos(x).
   *
   *  @param isSin        true for sin(x), false for cos(x)
   *  @param t            polynomial coefficients (from SinCosTable)
   *  @param x            input
   *  @param oneOverPiPad precision pad for the 1/pi constant; must match
   *                      the hardware. `None` uses `defaultOneOverPiPad`.
   */
  def apply(
    isSin: Boolean,
    t:     PolynomialCoefficients,
    x:     RealGeneric,
    oneOverPiPad: Option[Int] = None,
  ): RealGeneric = {

    val spec   = x.spec
    val manW   = spec.manW
    val exBias = spec.exBias

    if (x.isNaN)      { return RealGeneric.nan(spec) }
    if (x.isInfinite) { return RealGeneric.nan(spec) }

    val pad = oneOverPiPad.getOrElse(defaultOneOverPiPad(spec))

    // -----------------------------------------------------------------------
    // emulated pre-mult: x * (1/pi)
    //
    // 1/2 > 1/pi > 1/4, so with the precision `manW + 2 + pad`, the highest
    // set bit of the integer representation of 1/pi lands at position
    // `manW + pad`.

    val oneOverPi = (Real.one / Real.pi)(manW + 2 + pad)
    assert(bit(manW + pad, oneOverPi) == 1L)

    // No rounding here: we will later subtract 0.5 or 1.5 from this value,
    // so the low bits matter.
    val xOverPiProd          = x.manW1 * oneOverPi
    val xOverPiProdMoreThan2 = bit((1 + manW) + (1 + manW + pad) - 1, xOverPiProd)
    val xOverPiEx            = if (x.ex == 0) 0 else x.ex - 2 + xOverPiProdMoreThan2

    // extend the fraction one more bit (or not, depending on the overflow)
    val xOverPi      = xOverPiProd << (1 - xOverPiProdMoreThan2)
    val xOverPiFracW = (1 + manW) + (1 + manW + pad) - 1

    assert((SafeLong(1) << xOverPiFracW) < xOverPi &&
                                           xOverPi < (SafeLong(1) << (xOverPiFracW + 1)))

    // -----------------------------------------------------------------------
    // fold full-range x/pi into [0, 2)
    //
    // x/pi underflows when its exponent is at or below the smallest normal
    // one; the reduced argument is then zero, so sin(x) is +-0 and cos(x)
    // is +-1.

    val xOverPiIsZero = xOverPiEx <= 0

    val xOverPiAligned: SafeLong = if (xOverPiIsZero) {
      SafeLong(0)
    } else if (xOverPiEx >= exBias) {
      // drop bits representing ≥ 2
      slice(0, 1 + xOverPiFracW, xOverPi << (xOverPiEx - exBias))
    } else {
      xOverPi >> (exBias - xOverPiEx)
    }

    val xOverPiAligned2MSBs          = slice(1 + xOverPiFracW - 2, 2, xOverPiAligned)
    val xOverPiAlignedMoreThan3over2 = xOverPiAligned2MSBs == 3
    val xOverPiAlignedMoreThan1      = xOverPiAligned2MSBs == 2
    val xOverPiAlignedMoreThan1over2 = xOverPiAligned2MSBs == 1

    // -----------------------------------------------------------------------
    // sign of result, derived from the 2 MSBs of (|x|/pi mod 2)
    //
    // sin(x)                 cos(x)
    // 2MSBs 00 01 10 11      2MSBs 00 01 10 11
    //      |  :  :  :  :          |  :  :  :  :
    //      | .-. :  :  :        .-|-.:  :  :.-:
    // _____|'_:_': _:__:     __.__|__:__:__:__:
    //   _.'|  :  :'._.':          |  :'.:.':  :
    //      |  :  :  :  :          |  :  :  :  :
    //            '-----'             '-----'
    //            negative            negative

    val zSgn: Int = if (isSin) {
      val s = if (xOverPiAlignedMoreThan1 || xOverPiAlignedMoreThan3over2) 1 else 0
      if (x.sgn == 1) 1 - s else s
    } else {
      if (xOverPiAlignedMoreThan1 || xOverPiAlignedMoreThan1over2) 1 else 0
    }

    // -----------------------------------------------------------------------
    // fold further into y in [0, 1/2)
    //
    // sin: pi/2-x or pi-x or x-pi or 2pi-x
    // cos: pi/2-x or x-pi/2 or 3pi/2-x or x-3pi/2

    val (yex0, yman0): (Int, SafeLong) = if (isSin) {
      if (xOverPiAlignedMoreThan3over2 || xOverPiAlignedMoreThan1over2) {
        // pi/2 < |x|%(2pi) < pi    or    3/2 pi < |x|%(2pi) < 2pi
        (exBias,
         (SafeLong(1) << (1 + xOverPiFracW - 1))
           - slice(0, 1 + xOverPiFracW - 1, xOverPiAligned))
      } else if (xOverPiAlignedMoreThan1 || xOverPiIsZero ||
                 (exBias <= xOverPiEx && xOverPiAligned2MSBs == 0)) {
        (exBias, slice(0, 1 + xOverPiFracW - 1, xOverPiAligned))
      } else {
        // 0 < |x| < pi/2  -- use the un-aligned x/pi to keep precision
        assert(xOverPiEx < exBias, f"xOverPiEx = ${xOverPiEx}, exBias = ${exBias}")
        (xOverPiEx, xOverPi)
      }
    } else {
      if (xOverPiAlignedMoreThan3over2 || xOverPiAlignedMoreThan1over2) {
        (exBias, slice(0, 1 + xOverPiFracW - 2, xOverPiAligned))
      } else {
        (exBias,
         (SafeLong(1) << (1 + xOverPiFracW - 2))
           - slice(0, 1 + xOverPiFracW - 2, xOverPiAligned))
      }
    }

    // round/normalize yman0 to a manW-bit mantissa.
    val (yex, yman) = if (yman0 == SafeLong(0)) {
      (0, SafeLong(0))
    } else {
      val yman0W         = binaryWidthSL(yman0)
      val yman0Shift     = 1 + xOverPiFracW - yman0W
      val yman0Shifted   = yman0 << yman0Shift
      val yman0RoundBit  = xOverPiFracW - manW
      val yman0Rounded   = (yman0Shifted >> yman0RoundBit) +
                              bit(yman0RoundBit - 1, yman0Shifted)
      val yman0MoreThan2 = bit(manW + 1, yman0Rounded)
      assert(yman0MoreThan2 == 1 || bit(manW, yman0Rounded) == 1)
      val yexRounded     = yex0 - yman0Shift + yman0MoreThan2
      if (yexRounded == 0) (0, SafeLong(0))
      else                 (yexRounded.toInt, slice(0, manW, yman0Rounded))
    }

    assert(yex <= exBias - 1)
    assert(yex != exBias - 1 || yman == SafeLong(0))
    assert(yman < (SafeLong(1) << manW))

    // -----------------------------------------------------------------------
    // special small-y cases

    // y == 0.5 exactly  ->  result is ±1
    if (yex == exBias - 1) {
      assert(yman == SafeLong(0))
      return new RealGeneric(spec, zSgn, exBias, SafeLong(0))
    }
    // y == 0  ->  result is ±0
    if (yex == 0) {
      assert(yman == SafeLong(0), f"yex == 0, but yman == ${yman}")
      return new RealGeneric(spec, zSgn, 0, SafeLong(0))
    }

    // -----------------------------------------------------------------------
    // table lookup of  z(2y) = sin(Pi*y) / y  in [1, Pi/2)
    //
    // The address used by the table is the bit pattern of `2y` (a manW-bit
    // unsigned in [0, 2^manW)). Note y in [0, 1/2) here, so 2y in [0, 1).

    val ymanW1   = yman + (SafeLong(1) << manW)
    val yAligned = ymanW1 >> (exBias - 1 - yex)
    assert(bit(manW, yAligned) != 1)

    val nOrder  = t.bitW.nOrder
    val addrW   = t.bitW.addrW
    val outPrec = t.bitW.outPrec

    assert(addrW <= manW || nOrder == 0,
      f"addrW(${addrW}) > manW(${manW}), but nOrder is not zero.")
    assert(addrW == manW || nOrder != 0,
      f"addrW(${addrW}) must == manW(${manW}) when nOrder == 0.")

    val xW = manW - addrW
    val (addr, xiv) = if (nOrder == 0) {
      (yAligned.toInt, BigInt(0))
    } else {
      val a = slice(xW, addrW, yAligned).toInt
      val v = slice(0, xW, yAligned).toBigInt
      (a, v)
    }

    val res0 = PolynomialEval.eval(addr, xiv, t)

    val res = if (res0 < 0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x = ${x.toDouble}, sin(x) = ${sin(x.toDouble)}")
      SafeLong(0)
    } else if (res0 >= (SafeLong(1) << outPrec)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x = ${x.toDouble}, sin(x) = ${sin(x.toDouble)}")
      maskSL(outPrec)
    } else {
      SafeLong(res0)
    }

    // wmanW1 represents z = sin(Pi*y)/y in [1, Pi/2) as a fixed-point with
    // the hidden 1 at position outPrec.
    val wmanW1 = res + (SafeLong(1) << outPrec)

    // -----------------------------------------------------------------------
    // emulated post-mult: sin(x) = (sin(Pi*y)/y) * y = wmanW1 * ymanW1

    val zProd                = wmanW1 * ymanW1
    val zMoreThan2           = bit((manW + 1) + (outPrec + 1) - 1, zProd)
    val zRoundBit            = outPrec + zMoreThan2
    val zRounded             = roundBySpec(RoundSpec.roundToEven, zRoundBit, zProd)
    val zMoreThan2AfterRound = bit(manW + 1, zRounded)
    val zExInc               = zMoreThan2 + zMoreThan2AfterRound
    val zManW1               = if (zMoreThan2AfterRound == 1) SafeLong(1) << manW
                               else                            zRounded

    val zMan = slice(0, manW, zManW1)
    val wex  = exBias - 1
    val zEx  = wex - exBias + yex - exBias + zExInc + exBias + 2 // 2 for table

    new RealGeneric(spec, zSgn, zEx, zMan)
  }
}
