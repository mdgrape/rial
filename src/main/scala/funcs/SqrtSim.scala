//% @file sqrtSim.scala
//
// Simulators for sqrt function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.funcs

import scala.math._
import java.lang.Math.scalb

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.approx._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric

class SqrtTable(
  val nOrder:     Int,
  val nIntervals: Int,
  val inBitW:     Int,
  val outBitW:    Int,
) extends TableBase {

  // to distinguish ex=odd and ex=even cases, we use the LSB of exbit at the
  // top of the table address.
  // XXX: assuming exBias is odd. That means that, in case of exNobias == 0,
  //      the MSB is 1.
  //
  // [0.0, 0.5) -> [2, 4), sqrt([2, 4)) = [1.414, 2)
  // [0.5, 1.0) -> [1, 2), sqrt([1, 2)) = [1, 1.414)
  private val func = (x: Double) => {
    if(x < 0.5) {
      sqrt(x*4.0+2.0)-1.0
    } else {
      sqrt(x*2.0)-1.0
    }
  }

  // since we use different table for even and odd exponent, we double the
  // number of intervals
  private val ft = FuncTable(func, nOrder, nIntervals*2, 0.0, 1.0)
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW)

  // ------------------------------------------------------------------------

  def getBitWidth(): PolynomialBitWidth = {
    gen.getBitWidth()
  }
  def generateTable(
    bitW: PolynomialBitWidth
  ): PolynomialCoefficients = {
    gen.generateTable(bitW)
  }
}

object SqrtTable {

  def fp32(): SqrtTable = {

    val addrW = 8
    val extraBits = 3

    new SqrtTable(
      nOrder     = 2,
      nIntervals = (1 << addrW),
      inBitW     = (RealSpec.Float32Spec.manW - addrW),
      outBitW    = (RealSpec.Float32Spec.manW + extraBits),
    )
  }
  def bf16(): SqrtTable = {

    val addrW = RealSpec.BFloat16Spec.manW
    val extraBits = 0

    new SqrtTable(
      nOrder     =  0,
      nIntervals = (1 << addrW),
      inBitW     = (RealSpec.BFloat16Spec.manW - addrW),
      outBitW    = (RealSpec.BFloat16Spec.manW + extraBits),
    )
  }

  val sqrtDefaultFP32Table = fp32().generateTable(fp32().getBitWidth())

  /** sqrt simulation with the default configuration (2, 8, 3)
   */
  def sqrtDefaultFP32(x: Float): Float = {
    SqrtSimGeneric(
      sqrtDefaultFP32Table,
      RealGeneric.fromFloat(RealSpec.Float32Spec, x)
    ).toFloat
  }
}

/** An interface for software simulator
 */
object SqrtSimGeneric {

  // sqrt( [1, 4) ) -> [1, 2) (without the hidden bit)
  def apply( t: PolynomialCoefficients, x: RealGeneric ): RealGeneric = {

    if (x.isNaN)  {return RealGeneric.nan (x.spec)}
    if (x.isZero) {return RealGeneric.zero(x.spec)}
    if (x.isInfinite && x.sgn == 0) {return RealGeneric.inf (x.spec, sgn=0)}
    if (x.sgn == 1) {return RealGeneric.zero(x.spec)} // not NaN, for usability

//     println(f"x    = ${x.toDouble}(${sgn}|${x.ex}(${ex})|${x.man.toLong.toBinaryString})")
//     println(f"xman = ${man.toLong.toBinaryString}")

    val addrW   = t.bitW.addrW // XXX this is already +1-ed. see sqrtTableGeneration
    val nOrder  = t.bitW.nOrder
    val inPrec  = t.bitW.inPrec  // (1 <<  inPrec) == 1.0
    val outPrec = t.bitW.outPrec // (1 << outPrec) == 1.0

    val exW       = x.spec.exW
    val manW      = x.spec.manW
    val exBias    = x.spec.exBias
    val extraBits = outPrec - manW

    val sgn   = x.sgn
    val exRaw = x.ex
    val ex    = exRaw - exBias
    val man   = slice(0, manW+1, x.value) // x is even -> x.man, odd -> x.manW1

    assert(manW+1 == addrW + t.bitW.inBitW)

    val zSgn = 0
    val zEx  = (ex >> 1) + exBias

    assert(addrW <  manW || nOrder == 0, f"table address width(${addrW}) >= "+
      f"mantissa width(${manW}), but polynomial order is not zero.")
    assert(addrW == manW || nOrder != 0, f"table address width(${addrW}) "+
      f" must be == mantissa width(${manW}), for polynomial order is zero.")

    // val zman0 = if (nOrder==0) {
    //
    //   val addr = man.toInt
    //   PolynomialEval.eval(addr, BigInt(0), t)
    //
    // } else {

    val xW    = manW+1 - addrW
    val addr  = slice(xW, addrW, man).toInt
    val xiv   = slice( 0,    xW, man).toBigInt
    val zman0 = PolynomialEval.eval(addr, xiv, t)

    val zman = if(zman0 < 0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
      SafeLong(0)
    } else if (zman0 >= (SafeLong(1)<<outPrec)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
      maskSL(outPrec)
    } else {
      SafeLong(zman0)
    }

    // Simple rounding
    val zmanRounded = if (extraBits > 0) {
      val rounded = (zman >> extraBits) + bit(extraBits-1, zman)
      if(rounded >= (SafeLong(1) << manW)) {
        maskSL(manW)
      } else {
        rounded
      }
    } else {
      zman
    }
    val z = new RealGeneric(x.spec, zSgn, zEx, zmanRounded)
//     println(f"z    = ${z.toDouble}(${zSgn}|${z.ex}(${z.ex-exBias})|${z.man.toLong.toBinaryString})")
//     println(f"zman = ${zmanRounded.toLong.toBinaryString}")
    z
  }
}
