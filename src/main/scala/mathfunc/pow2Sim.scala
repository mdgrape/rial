//% @file pow2Sim.scala
//
// Simulators for pow2(x)
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import chisel3.util.log2Up

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.table._
import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.Rounding._
import rial.arith._

object MathFuncPow2Sim {

  def pow2SimGeneric( t : FuncTableInt, x : RealGeneric ): RealGeneric = {

//     println("==================================================")
    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val adrW      = t.adrW
    val fracW     = t.bp
    val extraBits = fracW - manW

    val xsgn   = x.sgn

//     println(f"x = ${x.sgn}|${x.ex}(${x.ex})|${x.man.toLong.toBinaryString}")
//     println(f"Pow2Sim: x = ${x.toDouble}, pow2(x) = ${pow2(x.toDouble)}")

    val xnan  = x.isNaN
    val xinf  = x.isInfinite
    val xzero = x.isZero

    if(xnan) {
      return RealGeneric.nan(x.spec)
    }
    if(xinf && xsgn == 0) {
      return RealGeneric.inf(x.spec, 0)
    }
    if(xinf && xsgn == 1) {
      return new RealGeneric(x.spec, 1.0)
    }

    val xexNobias = x.ex - exBias

    // xint  = if(x>0) {floor(x)} else {ceil(x)}
    // xfrac = x - xint
    //
    // 2^x = 2^xint * 2^xfrac
    //
    // xfrac is in [0, 1), so 2^xfrac is in [1, 2). So
    //
    // ex  = xint
    // man = 2^xfrac - 1
    //
    // xint is in [2^x.ex, 2^(x.ex+1) ).

    val log2 = (a:Double) => {log(a) / log(2.0)}

    // if xex exceeds this, the result overflows
    val xExOvfLimit = math.ceil(log2(maskL(exW)-exBias)).toLong // log2(255-127 = 128) = 7
    val xExUdfLimit = math.ceil(log2(abs(0 - exBias)  )).toLong // log2(|0-127| = 127) > 6
    // check if x is large enough to cause over/under flow
    if (xsgn == 0 && xExOvfLimit <= xexNobias) {
      return RealGeneric.inf(x.spec, 0)
    } else if(xsgn == 1 && xExUdfLimit <= xexNobias) {
      return RealGeneric.zero(x.spec)
    }
    // Note: still the result might over/underflows because log2(...) might not be
    //       an integer value.

    val padding = extraBits

    // If integer part of x is larger than 2^ovflim, then it overflows.
    // So we need integer part bitwidth smaller than the limit.
    val xIntW  = max(xExOvfLimit, xExUdfLimit).toInt
    val xFracW = manW + padding

    //  .---xIntW----. .--- xFracW ----.
    // | integer part | fractional part |

    val xValW = xIntW + xFracW + 1
    val xVal  = x.manW1 << (xIntW + padding)

    assert(xVal < maskL(xValW))
    assert(xexNobias < xIntW)

    val xshift = xIntW - xexNobias
    val xValShifted = (xVal >> xshift) + bit(xshift-1, xVal) // simple rounding
    val xint0  = slice(xFracW, xIntW,  xValShifted).toLong
    val xfrac0 = slice(0,      xFracW, xValShifted).toLong
    val (xint, xfrac) = if(xsgn == 0) {
      (xint0, xfrac0)
    } else {
      if(xfrac0 != 0) {
        (-xint0 - 1L, (1L<<xFracW) - xfrac0)
      } else {
        (-xint0, 0L)
      }
    }
    val dxbp = manW-adrW-1
    val d    = slice(padding, dxbp+1, xfrac) - (SafeLong(1) << dxbp) // DO NOT round here, because we have correction term
    val adr  = slice(padding+dxbp+1, adrW, xfrac)

    val zman = t.interval(adr.toInt).eval(d.toLong, dxbp)
    assert(0 <= zman && zman < (1<<fracW))

    val zmanRound = if(extraBits > 0) {
      (zman >> extraBits) + bit(extraBits-1, zman)
    } else {
      zman
    }

    val corrTermW   = 10 // XXX determined empirically
    val coefficient = new RealGeneric(x.spec, log(2.0))
    val xShiftedOut = slice(0, padding, xfrac)
    val zCorrectionCoef0 = (coefficient.manW1 * xShiftedOut)
    val zCorrectionCoef  = zCorrectionCoef0 >> ((manW + 1) + padding - corrTermW)
    val zCorrectionTerm0 = (((1<<manW) + zmanRound) * zCorrectionCoef)
    val zCorrectionTerm  = zCorrectionTerm0 >> (manW+1 + corrTermW - corrTermW)
    val zCorrectionShifted = bit(corrTermW-1, zCorrectionTerm) + bit(corrTermW-2, zCorrectionTerm)
    assert(0 <= zCorrectionShifted && zCorrectionShifted <= 1)

    val zmanCorrected = zmanRound + zCorrectionShifted
    val zmanMoreThan2AfterCorrection = bit(manW, zmanCorrected)
    val z = slice(0, manW, zmanCorrected)

    // 2^(-x) = 2^(-xInt - xFrac)
    //        = 2^(-xInt - 1 + 1 - xFrac)
    //        = 2^(-xInt - 1 + (1 - xFrac))
    //                         ^^^^^^^^^^ positive

    val zex = xint + exBias + zmanMoreThan2AfterCorrection

    if      (zex>=maskI(exW)) { return RealGeneric.inf (x.spec, 0) }
    else if (zex<=0)          { return RealGeneric.zero(x.spec) }

    new RealGeneric(x.spec, 0, zex.toInt, SafeLong(z))
  }
}
