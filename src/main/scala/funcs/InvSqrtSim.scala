//% @file InvSqrtSim.scala
//
// Simulators for 1/sqrt function
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

class InvSqrtTable(
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
  // [0.0, 0.5) -> [2, 4), 1/sqrt([2, 4)) = [0.5, 0.707)
  // [0.5, 1.0) -> [1, 2), 1/sqrt([1, 2)) = [0.707, 1.0)
  private val func = (x: Double) => {
    if(x < 0.5) {
      (2.0/sqrt(x*4.0+2.0)) - 1.0
    } else {
      (2.0/sqrt(x*2.0)) - 1.0
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

object InvSqrtTable {

  def fp32(): InvSqrtTable = {

    val addrW = 8
    val extraBits = 3

    new InvSqrtTable(
      nOrder     = 2,
      nIntervals = (1 << addrW),
      inBitW     = (RealSpec.Float32Spec.manW - addrW),
      outBitW    = (RealSpec.Float32Spec.manW + extraBits),
    )
  }
  def bf16(): InvSqrtTable = {

    val addrW = RealSpec.BFloat16Spec.manW
    val extraBits = 0

    new InvSqrtTable(
      nOrder     =  0,
      nIntervals = (1 << addrW),
      inBitW     = (RealSpec.BFloat16Spec.manW - addrW),
      outBitW    = (RealSpec.BFloat16Spec.manW + extraBits),
    )
  }

  val invsqrtDefaultFP32Table = fp32().generateTable(fp32().getBitWidth())

  /** invsqrt simulation with the default configuration (2, 8, 3)
   */
  def invsqrtDefaultFP32(x: Float): Float = {
    InvSqrtSimGeneric(
      invsqrtDefaultFP32Table,
      RealGeneric.fromFloat(RealSpec.Float32Spec, x)
    ).toFloat
  }
}

/** An interface for software simulator
 */
object InvSqrtSimGeneric {

  // 1/sqrt( [1, 4) ) -> [0.5, 1) (without the hidden bit)
  def apply( t: PolynomialCoefficients, x: RealGeneric ): RealGeneric = {

    if (x.isNaN)  {return RealGeneric.nan (x.spec)}
    if (x.isZero) {return RealGeneric.inf (x.spec, sgn=0)}
    if (x.isInfinite) {return RealGeneric.zero(x.spec)}
    if (x.sgn == 1) {return RealGeneric.inf(x.spec, sgn=0)} // sqrt(-x) = inf

    val addrW   = t.bitW.addrW // XXX this is already +1-ed. see invsqrtTableGeneration
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
    val zEx0 = exBias - (ex >> 1) - 1

    // invsqrt uses manW+1 address bits (including LSB of exponent for even/odd distinction)
    assert(addrW <= manW+1 || nOrder == 0, f"table address width(${addrW}) > "+
      f"mantissa width+1(${manW+1}), but polynomial order is not zero.")
    assert(addrW == manW+1 || nOrder != 0, f"table address width(${addrW}) "+
      f" must be == mantissa width+1(${manW+1}), for polynomial order is zero.")

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

    // Special case: if x.man == 0 && x.ex is odd, then z.man = 0 and z.ex += 1
    val (zEx, zmanFinal) = if(x.ex % 2 == 1 && x.man == SafeLong(0) && !x.isZero) {
      (zEx0+1, SafeLong(0))
    } else {
      (zEx0, zmanRounded)
    }

    val z = new RealGeneric(x.spec, zSgn, zEx, zmanFinal)
    z
  }
}
