//% @file ReciprocalSim.scala
//
// Simulators for reciprocal function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.math._
import java.lang.Math.scalb

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.table._

import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric

class ReciprocalTable(
  val nOrder:  Int,
  val addrW:   Int,
  val inBitW:  Int,
  val outBitW: Int,
) extends TableBase {

  // 1/x = 2^(-e-1) * 2/1.m
  // We compute 2/(1+x) - 1 where x is the normalized mantissa
  private val func = (x: Double) => {
    2.0/(1.0+x) - 1.0
  }

  private val nIntervals = 1 << addrW
  private val ft = FuncTable(func, nOrder, nIntervals, 0.0, 1.0)
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

object ReciprocalTable {

  def fp32(): ReciprocalTable = {

    val addrW = 8
    val extraBits = 3

    new ReciprocalTable(
      nOrder  = 2,
      addrW   = addrW,
      inBitW  = (RealSpec.Float32Spec.manW - addrW),
      outBitW = (RealSpec.Float32Spec.manW + extraBits),
    )
  }
  def bf16(): ReciprocalTable = {

    val addrW = RealSpec.BFloat16Spec.manW
    val extraBits = 0

    new ReciprocalTable(
      nOrder  = 0,
      addrW   = addrW,
      inBitW  = (RealSpec.BFloat16Spec.manW - addrW),
      outBitW = (RealSpec.BFloat16Spec.manW + extraBits),
    )
  }

  val reciprocalDefaultFP32Table = fp32().generateTable(fp32().getBitWidth())

  /** reciprocal simulation with the default configuration (2, 8, 3)
   */
  def reciprocalDefaultFP32(x: Float): Float = {
    ReciprocalSimGeneric(
      reciprocalDefaultFP32Table,
      RealGeneric.fromFloat(RealSpec.Float32Spec, x)
    ).toFloat
  }
}

/** An interface for software simulator
 */
object ReciprocalSimGeneric {

  // 1/x: compute 2/(1+x) - 1
  def apply( t: PolynomialCoefficients, x: RealGeneric ): RealGeneric = {

    if (x.isNaN)  {return RealGeneric.nan (x.spec)}
    if (x.isZero) {return RealGeneric.inf (x.spec, sgn=0)}
    if (x.isInfinite) {return RealGeneric.zero(x.spec)}

    val addrW   = t.bitW.addrW
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
    val man   = slice(0, manW, x.value)

    val zSgn = sgn

    // 1/x = 2^(-e-1) * 2/1.m
    val xmanNonZero = if(man != SafeLong(0)) {1} else {0}
    val zEx0 = (exBias << 1) - exRaw - xmanNonZero

    // Check overflow/underflow
    if (zEx0 >= maskI(exW+1)) {
      return RealGeneric.inf(x.spec, sgn)
    } else if (zEx0 <= 0) {
      return RealGeneric.zero(x.spec)
    }

    val zEx = zEx0

    assert(addrW <  manW || nOrder == 0, f"table address width(${addrW}) >= "+
      f"mantissa width(${manW}), but polynomial order is not zero.")
    assert(addrW == manW || nOrder != 0, f"table address width(${addrW}) "+
      f" must be == mantissa width(${manW}), for polynomial order is zero.")

    // If man == 0, that means x = 2^N, then z = 2^-N, so zman = 0
    val zman0 = if (man == SafeLong(0)) {
      SafeLong(0)
    } else {
      val xW    = manW - addrW
      val addr  = slice(xW, addrW, man).toInt
      val xiv   = slice( 0,    xW, man).toBigInt
      SafeLong(PolynomialEval.eval(addr, xiv, t))
    }

    val zman = if(zman0 < SafeLong(0)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
      SafeLong(0)
    } else if (zman0 >= (SafeLong(1)<<outPrec)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
      maskSL(outPrec)
    } else {
      zman0
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
    z
  }
}
