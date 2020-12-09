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

object FloatChiselUtil {

  def decompose( spec: RealSpec, x: UInt ) = {
    if (x.getWidth != spec.W) {
      throw new RuntimeException(f"ERROR (${this.getClass.getName}) : Width mismatch (expected ${spec.W}, actual ${x.getWidth})")
    }
    val sgn = if (spec.disableSign) 0.U else x(spec.W-1)
    val ex  = x(spec.manW+spec.exW-1,spec.manW)
    val man = x(spec.manW-1,0)
    (sgn, ex, man)
  }

  def decomposeWithLeading1( spec: RealSpec, x: UInt ) = {
    if (x.getWidth != spec.W) {
      throw new RuntimeException(f"ERROR (${this.getClass.getName}) : Width mismatch (expected ${spec.W}, actual ${x.getWidth})")
    }
    val sgn = if (spec.disableSign) 0.U else x(spec.W-1)
    val ex  = x(spec.manW+spec.exW-1,spec.manW)
    val man = 1.U(1.W) ## x(spec.manW-1,0)
    (sgn, ex, man)
  }

  def checkValue( spec: RealSpec, x: UInt ) = {
    val ex  = x(spec.manW+spec.exW-1,spec.manW)
    val man = x(spec.manW-1,0)
    val inf  = ex.andR.asBool
    val zero = !(ex.orR.asBool)
    val nan  = (!spec.disableNaN).B && inf && man.orR.asBool
    (zero, inf, nan)
  }

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
