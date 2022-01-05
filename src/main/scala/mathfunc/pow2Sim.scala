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
    //

    val log2 = (x:Double) => {log(x) / log(2.0)}

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

    // If integer part of x is larger than 2^ovflim, then it overflows.
    // So we need integer part bitwidth smaller than the limit.
    val xIntW  = max(xExOvfLimit, xExUdfLimit).toInt
    val xFracW = manW + extraBits // extrabits is for rounding. after shift it, round it before extracting dx.

    //  .---xIntW----. .--- xFracW ----.
    // | integer part | fractional part |

    val xValW = xIntW + xFracW + 1
    val xVal  = x.manW1 << (xIntW + extraBits)

    assert(xVal < maskL(xValW))
    assert(xexNobias < xIntW)

    val xshift = xIntW - xexNobias
    val xValShifted = xVal >> xshift

    val xint  = slice(xFracW, xIntW,  xValShifted)
    val xfrac = slice(0,      xFracW, xValShifted)

    val zex  = if(xsgn == 0) {
       xint + exBias
    } else { // x is negative!
      -xint + exBias
    }

    if      (zex>=maskI(exW)) { return RealGeneric.inf (x.spec, 0) }
    else if (zex<=0)          { return RealGeneric.zero(x.spec) }

    val dxbp = manW-adrW-1
    val d    = slice(extraBits, dxbp+1, xfrac) - // round!
               (SafeLong(1) << dxbp)
    val adr  = slice(extraBits+dxbp+1, adrW, xfrac)

    val zman = t.interval(adr.toInt).eval(d.toLong, dxbp)

    val zmanRound = if(extraBits > 0) {
      (zman >> extraBits) + bit(extraBits-1, zman)
    } else {
      zman
    }

    val z = if (zman<0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
      0L
    } else if (zmanRound >= (1L<<manW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
      maskL(manW)
    } else {
      zmanRound
    }
    new RealGeneric(x.spec, 0, zex.toInt, SafeLong(z))
  }
}
