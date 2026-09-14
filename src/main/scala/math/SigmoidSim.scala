//% @file SigmoidSim.scala
//
// Simulator for sigmoid(x) = 1 / (1 + exp(-x)).
//
// Algorithm:
//   1. saturate large |x|: sigmoid(|x|>=rangeMax) ~ 1 (or 0 if x < 0).
//   2. very small |x| (xScaled <= 1): result is 0.5.
//   3. otherwise look up the polynomial table at addr=xScaled.upperBits;
//      the table approximates  2 - 2/(1+exp(-|x|)) = 2*(1 - sigmoid(|x|))
//      for |x| in [0, rangeMax).
//   4. fold the table value into the final answer:
//        for x > 0:  sigmoid(x) = 1 - tableValue / 2
//        for x < 0:  sigmoid(x) =     tableValue / 2
//      then normalize the (calcW+1)-bit numerator into a (1+manW)-bit mantissa.
//
// `rangeMaxLog2` is chosen so that  2 - 2/(1+exp(-rangeMax)) < 2^-outBitW,
// i.e., the table function falls below one LSB beyond the table range -- so
// the sigmoid result is indistinguishable from saturation in the output FP.
//
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._

import spire.math.SafeLong
import spire.implicits._

import rial.table._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric

class SigmoidTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  // [0, 1) -> |x| in [0, rangeMax) -> [1, 0) (approximately).
  val rangeMaxLog2 = SigmoidTable.rangeMaxLog2(outBitW)
  val rangeMax     = 1 << rangeMaxLog2

  private val func = (x01: Double) => {
    val x = x01 * rangeMax
    val z = 2.0 - (2.0 / (1.0 + exp(-x)))
    assert(z <= 1.0, f"x = ${x}, z = ${z}")
    z
  }

  private val nIntervals = 1 << addrW
  private val ft  = FuncTable(func, nOrder, nIntervals, 0.0, 1.0)
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW)

  def getBitWidth(): PolynomialBitWidth = { gen.getBitWidth() }
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    gen.generateTable(bitW)
  }
}

object SigmoidTable {

  /** Smallest power-of-two threshold above which sigmoid saturates within
   *  `precW` bits of precision.  Mirrors the original math/ port.
   */
  def rangeMaxLog2(precW: Int): Int = {
    val log2  = (x: Double) => log(x) / log(2.0)
    val delta = pow(2.0, -precW)
    // we want `2 - 2/(1+exp(-rangeMax)) < delta` so that the table value at
    // the upper end is smaller than one LSB of the output precision.
    val threshold = -log(2.0 / (2.0 - delta) - 1.0)
    ceil(log2(threshold)).toInt
  }

  def fp32(): SigmoidTable = {
    val addrW     = 8
    val extraBits = 3
    new SigmoidTable(
      nOrder  = 2,
      addrW   = addrW,
      inBitW  = (RealSpec.Float32Spec.manW - addrW),
      outBitW = (RealSpec.Float32Spec.manW + extraBits),
    )
  }
  def bf16(): SigmoidTable = {
    val addrW     = RealSpec.BFloat16Spec.manW
    val extraBits = 0
    new SigmoidTable(
      nOrder  = 0,
      addrW   = addrW,
      inBitW  = (RealSpec.BFloat16Spec.manW - addrW),
      outBitW = (RealSpec.BFloat16Spec.manW + extraBits),
    )
  }

  private lazy val fp32Table = fp32().generateTable(fp32().getBitWidth())
  def sigmoidDefaultFP32(x: Float): Float = {
    SigmoidSimGeneric(fp32Table,
      RealGeneric.fromFloat(RealSpec.Float32Spec, x)).toFloat
  }
}

/** An interface for software simulator.
 */
object SigmoidSimGeneric {

  def apply(t: PolynomialCoefficients, x: RealGeneric): RealGeneric = {

    if (x.isNaN)                    { return RealGeneric.nan (x.spec)      }
    if (x.isInfinite && x.sgn == 0) { return new RealGeneric (x.spec, 1.0) }
    if (x.isInfinite && x.sgn == 1) { return RealGeneric.zero(x.spec)      }

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val addrW     = t.bitW.addrW
    val outPrec   = t.bitW.outPrec
    val calcW     = outPrec
    val extraBits = outPrec - manW

    val rangeMaxLog2 = SigmoidTable.rangeMaxLog2(outPrec)
    val rangeMaxEx   = exBias + rangeMaxLog2 - 1

    // ---- saturation: large |x| ----
    if (x.ex > rangeMaxEx) {
      if (x.sgn == 0) { return new RealGeneric(spec, 1.0) }
      else            { return RealGeneric.zero(spec)     }
    }

    // ---- scale |x| into the [0, rangeMax) fixed-point address space ----
    val xScaled = x.manW1 >> (rangeMaxEx - x.ex + 1)

    // ---- very small |x| shortcut: sigmoid ~ 0.5 ----
    if (xScaled == SafeLong(0) || xScaled == SafeLong(1)) {
      return new RealGeneric(spec, 0.5)
    }

    // ---- polynomial table query ----
    val addr = TableAddress.adr(xScaled, manW, addrW)
    val xiv  = TableAddress.dx(xScaled, manW, addrW, t.bitW.inBitW)
    val zScaled0 = PolynomialEval.eval(addr, xiv, t)

    val zScaled = if (zScaled0 < 0) {
      println(f"WARNING (${this.getClass.getName}) : polynomial value negative at x=${x.toDouble}")
      SafeLong(0)
    } else if (zScaled0 >= (SafeLong(1) << calcW)) {
      println(f"WARNING (${this.getClass.getName}) : polynomial range overflow at x=${x.toDouble}")
      maskSL(calcW)
    } else {
      SafeLong(zScaled0)
    }

    // ---- fold and normalize ----
    //   zScaled  =  2 * (1 - sigmoid(|x|))    in [0, 1) fixed-point at `calcW`
    //   x > 0  =>  sigmoid(x)  =  1 - zScaled / 2  =  (2 - zScaled) / 2
    //   x < 0  =>  sigmoid(x)  =  zScaled / 2
    val zShifted = if (x.sgn == 0) { (SafeLong(1) << (calcW + 1)) - zScaled }
                   else            { zScaled }

    val zShift   = (calcW + 2) - binaryWidthSL(zShifted)
    val zmanW1   = zShifted << zShift
    val zmanW1Rounded = (zmanW1 >> (extraBits + 1)) + bit(extraBits, zmanW1)

    assert(zShift >= 0)

    val zMan = if (zmanW1Rounded >= (SafeLong(1) << (manW + 1))) {
      println(f"WARNING (${this.getClass.getName}) : rounding overflow at x=${x.toDouble}, zmanW1 = ${zmanW1Rounded.toLong.toBinaryString}")
      maskSL(manW)
    } else {
      slice(0, manW, zmanW1Rounded)
    }

    val zEx = exBias - zShift
    new RealGeneric(spec, 0, zEx, zMan)
  }
}
