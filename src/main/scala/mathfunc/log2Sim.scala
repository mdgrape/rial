//% @file log2Sim.scala
//
// Simulators for log2(x)
// Copyright (C) Toru Niina RIKEN BDR 2022
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

object MathFuncLog2Sim {

  //
  // log2(2^ex * 1.man) = log2(2^ex) + log2(1.man)
  //                    = ex + log2(1.man)
  //
  def log2SimGeneric( t : FuncTableInt, x : RealGeneric ): RealGeneric = {

//     println("==================================================")
    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val adrW      = t.adrW
    val fracW     = t.bp
    val extraBits = fracW - manW

    val log2 = (a:Double) => {log(a) / log(2.0)}

//     println(f"x = ${x.sgn}|${x.ex}(${x.ex})|${x.man.toLong.toBinaryString}")
//     println(f"Log2Sim: x = ${x.toDouble}, log2(x) = ${log2(x.toDouble)}")

    // --------------------------------------------------------------------------
    // check special cases

    // - log2(nan) -> nan
    // - log2(inf) -> inf
    // - log2(0)   -> -inf
    // - log2(1)   -> 0
    // - log2(neg) -> nan

    val xnan  = x.isNaN
    val xinf  = x.isInfinite
    val xzero = x.isZero
    val xneg  = x.sgn == 1

    if(xnan) {
      return RealGeneric.nan(x.spec)
    }
    if(xinf && !xneg) {
      return RealGeneric.inf(x.spec, /*sgn = */0)
    }
    if(xzero) {
      return RealGeneric.inf(x.spec, /*sgn = */1)
    }
    if(xneg) {
      return RealGeneric.nan(x.spec)
    }

    val zsgn = if(x.ex < exBias) {1} else {0}

    // --------------------------------------------------------------------------

    val xexNobias = x.ex.toLong - exBias.toLong
    val zint0     = if(xexNobias >= 0) {xexNobias} else {-xexNobias - 1}

//     println(f"xexNobias = ${xexNobias}")
//     println(f"zint0     = ${zint0}")

    // --------------------------------------------------------------------------
    // taylor
    //
    // log(1+x) = x/ln(2) - x^2/2ln(2) + x^3/3ln(2) + O(x^4)
    //          = x(1 - x/2 + x^3/3) / ln(2)
    //
    // if x < 2^-8, then x^3 < 2^-24
    // use x(1 - x/2)ln(2)

    val xman = x.man.toLong // x - 1
    if((xexNobias == 0 || xexNobias == -1) && xman < (1 << (manW-8))) {

      if(xman == 0) {
        return RealGeneric.zero(x.spec)
      }

      val invln2 = math.round((1.0 / log(2.0)) * (1 << (manW+extraBits))).toLong // > 1
      val xmanbp = xman.toBinaryString.length
      val xln2prod  = invln2 * xman // xmanbp+1+manW+extrabits
      val xln2MoreThan2 = bit(xmanbp + manW + extraBits, xln2prod)
      val xln2shift = (xmanbp-1) + xln2MoreThan2
      val xln2      = xln2prod >> xln2shift

//       println(f"x/ln2 = ${xln2.toLong.toBinaryString}")
      assert(bit(manW+extraBits, xln2) == 1 && (xln2 >> (manW+extraBits+1)) == 0)

      val oneMinusHalfx = (1L << (manW+extraBits)) - (xman << (extraBits-1)) // < 1

//       println(f"oneMinusHalfx = ${oneMinusHalfx.toLong.toBinaryString}")
      assert(bit(manW+extraBits-1, oneMinusHalfx) == 1 && slice(manW+extraBits, 64 - (manW+extraBits), oneMinusHalfx) == 0)

      val resProd    = xln2 * oneMinusHalfx // W = 2 * (manW+extraBits) + 1
//       println(f"resProd    = ${resProd.toLong.toBinaryString}")

      val resProdMoreThan1 = bit(2*(manW+extraBits), resProd)
      val resShift   = (manW + extraBits - 1 + resProdMoreThan1)
      val resShifted = resProd >> resShift
//       println(f"resShifted = ${resShifted.toLong.toBinaryString}")
//       println(f"resProdMT1 = ${resProdMoreThan1}")
      assert(bit(manW+extraBits, resShifted) == 1)
      assert(resShifted >> (manW+extraBits+1) == 0)

      val zmanTaylor = (resShifted >> extraBits) + bit(extraBits-1, resShifted)
      val zexTaylor = -(manW - (xmanbp-1)) + xln2MoreThan2 - 1 + resProdMoreThan1

      return new RealGeneric(x.spec, zsgn, zexTaylor.toInt + exBias, zmanTaylor)
    }

    // --------------------------------------------------------------------------
    // polynomial

    val dxbp = manW-adrW-1
    val d    = slice(0,      dxbp+1, x.man) - (SafeLong(1) << dxbp)
    val adr  = slice(dxbp+1, adrW,   x.man)

    val zfrac0Pos = t.interval(adr.toInt).eval(d.toLong, dxbp)
    val zfrac0 = if(xexNobias >= 0) {zfrac0Pos} else {(1<<fracW) - zfrac0Pos}
//     println(f"s: zfrac0Pos = ${zfrac0Pos.toBinaryString}")
//     println(f"s: zfrac0    = ${zfrac0.toBinaryString   }")

    assert(0L <= zfrac0 && zfrac0 < (1L<<fracW))

    val zfull0 = (zint0 << fracW) + zfrac0.toLong
//     println(f"zfull0  = ${zfull0.toBinaryString }")
    assert(0 <= zfull0)

    val zfullW  = zfull0.toBinaryString.length
    val zShiftW = exW + fracW - zfullW
//     println(f"s: zfullW  = ${zfullW }")
//     println(f"s: zShiftW = ${zShiftW}")
    assert(0 <= zShiftW)

    val zShifted = zfull0 << zShiftW
//     println(f"s: zShifted = ${zShifted.toBinaryString }")
    assert(bit(exW + fracW - 1, zShifted) == 1)

    val zman0       = slice(exW-1, fracW, zShifted) // -1 for the hidden bit
    val zmanRounded = slice(extraBits, manW, zman0) + bit(extraBits-1, zman0)
//     println(f"s: zman0 = ${slice(extraBits, manW, zman0).toBinaryString } + ${bit(extraBits-1, zman0)}")
//     println(f"zman0 = ${zman0.toBinaryString }")
//     println(f"zmanR = ${zmanRounded.toBinaryString }")

    val zex0  = exBias + (exW-1) - zShiftW

    val zman = if (zmanRounded<0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
      0L
    } else if (zmanRounded >= (1L<<manW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
      maskL(manW)
    } else {
      zmanRounded.toLong
    }
    val zex = zex0.toInt

    new RealGeneric(x.spec, zsgn, zex, SafeLong(zman))
  }
}
