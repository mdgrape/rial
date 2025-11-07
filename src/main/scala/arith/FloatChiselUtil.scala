//% @file FloatChiselUtil.scala
//
// Some functions to be used float arithmetics
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.arith

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._

/** utility functions for floating-point handling.
 */
object FloatChiselUtil {

  /** decompose floating point into (sign, exponent, and mantissa) based on RealSpec.
   *
   * {{{
   * val (sgn, ex, man) = FloatChiselUtil.decompose(spec, f)
   * }}}
   *
   * @param spec floating point spec of input UInt.
   * @param x    floating point value.
   * @return a tuple of (sign, exponent, mantissa).
   */
  def decompose( spec: RealSpec, x: UInt ) = {
    if (x.getWidth != spec.W) {
      throw new RuntimeException(f"ERROR (${this.getClass.getName}) : Width mismatch (expected ${spec.W}, actual ${x.getWidth})")
    }
    val sgn = if (spec.disableSign) 0.U else x(spec.W-1)
    val ex  = x(spec.manW+spec.exW-1,spec.manW)
    val man = x(spec.manW-1,0)
    (sgn, ex, man)
  }

  /** decompose floating point into (sign, exponent, and mantissa) based on RealSpec.
   *
   * unlike [[rial.arith.FloatChiselUtil.decompose]], it adds hidden 1 to the decomposed mantissa.
   *
   * {{{
   * val (sgn, ex, manWith1) = FloatChiselUtil.decompose(spec, f)
   * }}}
   *
   * @param spec floating point spec of input UInt.
   * @param x    floating point value.
   * @return a tuple of (sign, exponent, mantissa). mantissa has the hidden bit at the MSB.
   */
  def decomposeWithLeading1( spec: RealSpec, x: UInt ) = {
    if (x.getWidth != spec.W) {
      throw new RuntimeException(f"ERROR (${this.getClass.getName}) : Width mismatch (expected ${spec.W}, actual ${x.getWidth})")
    }
    val sgn = if (spec.disableSign) 0.U else x(spec.W-1)
    val ex  = x(spec.manW+spec.exW-1,spec.manW)
    val man = 1.U(1.W) ## x(spec.manW-1,0)
    (sgn, ex, man)
  }

  /** checks if a floating point value is zero, inf, or nan.
   *
   * {{{
   * val (isZero, isInf, isNan) = FloatChiselUtil.checkValue(spec, f)
   * }}}
   *
   * @param spec floating point spec of input UInt.
   * @param x    floating point value.
   * @return a tuple of (isZero, isInf, isNan).
   */
  def checkValue( spec: RealSpec, x: UInt ) = {
    val ex  = x(spec.manW+spec.exW-1,spec.manW)
    val man = x(spec.manW-1,0)
    val inf  = ex.andR.asBool
    val zero = !(ex.orR.asBool)
    val nan  = (!spec.disableNaN).B && inf && man.orR.asBool
    (zero, inf, nan)
  }

  /** calculate carry bit according to [[rial.arith.RoundSpec]].
   *
   * @param spec   rounding spec.
   * @param lsb    the lsb of the fp.
   * @param round  the rounded bit.
   * @param sticky the sticky bit.
   * @param sgn    the sign of the rounded floating point.
   * @return carry bit.
   *
   * @see [[rial.arith.RoundSpec]]
   */
  def roundIncBySpec( spec: RoundSpec, lsb: Bool, round: Bool, stickey: Bool, sgn: Bool = false.B ) : Bool = {
    // All should be width 1
    val res = spec match {
      case RoundSpec.round          => round
      case RoundSpec.roundToEven    => round && (stickey || lsb)
      case RoundSpec.truncate       => false.B
      case RoundSpec.truncateToZero => sgn
      case RoundSpec.ceil           => round||stickey
      case RoundSpec.ceilToInfinite => sgn^(round||stickey)
      case _ => false.B
    }
    res
  }

}
