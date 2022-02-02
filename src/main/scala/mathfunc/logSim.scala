//% @file logSim.scala
//
// Simulators for log(x)
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

object MathFuncLogSim {

  def logSimGeneric( t : FuncTableInt, tSmallPos : FuncTableInt, tSmallNeg : FuncTableInt, x : RealGeneric ): RealGeneric = {

//     println("==================================================")
    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val adrW      = t.adrW
    val fracW     = t.bp
    val extraBits = fracW - manW

    val log2 = (a:Double) => {log(a) / log(2.0)}

//     println(f"x = ${x.sgn}|${x.ex}(${x.ex})|${x.man.toLong.toBinaryString}")
//     println(f"LogSim: x = ${x.toDouble}, log(x) = ${log(x.toDouble)}")

    // --------------------------------------------------------------------------
    // check special cases

    // - log(nan) -> nan
    // - log(inf) -> inf
    // - log(0)   -> -inf
    // - log(1)   -> 0
    // - log(neg) -> nan

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
    if(x.ex == exBias && x.man == 0) {
      return RealGeneric.zero(x.spec)
    }

    val zsgn = if(x.ex < exBias) {1} else {0}

    // --------------------------------------------------------------------------

    val xexNobias = x.ex.toLong - exBias.toLong
    val zint0     = if(xexNobias >= 0) {xexNobias} else {-xexNobias - 1}

//     println(f"xexNobias = ${xexNobias}")
//     println(f"zint0     = ${zint0}")

    val taylorThreshold = MathFuncLog2Sim.calcTaylorThreshold(x.spec)

    val (log2ex0, log2man0) = if(xexNobias == 0 && x.man.toLong < (1L << (manW-taylorThreshold))) {

      // log(1+x) = x/ln(2) - x^2/2ln(2) + x^3/3ln(2) + O(x^4)
      //          = x(1 - x/2 + x^2/3) / ln(2)

      val xman   = x.man // x - 1
      val xmanbp = xman.toLong.toBinaryString.length

      // 1/ln2 > 1
      val invln2   = math.round((1.0 / log(2.0)) * (1L << fracW)).toLong
      val oneThird = math.round((1.0 / 3.0)      * (1L << fracW)).toLong

      // 1 - x/2 < 1
      val oneMinusHalfx = (1L << fracW) - (xman << (extraBits-1))

      // x^2/3
      val xsq      = xman * xman
      val xsqThird = ((xsq * oneThird) >> (xmanbp + xmanbp)) >> ((manW - xmanbp) * 2)

      // 1 - x/2 + x^2/3 < 1
      val taylorTerm = oneMinusHalfx + xsqThird

      val lntermProd      = xman * taylorTerm
      val lntermMoreThan2 = bit(xmanbp + manW + extraBits, lntermProd)
      val lnterm          = lntermProd >> (xmanbp-1 + lntermMoreThan2)

      assert(bit(manW+extraBits-1, oneMinusHalfx) == 1 && slice(manW+extraBits, 64 - (manW+extraBits), oneMinusHalfx) == 0)

      val resProd    = invln2 * lnterm // W = 2 * (manW+extraBits) + 2

      val resProdMoreThan2 = bit(2*(manW+extraBits)+1, resProd)
      val resShift   = (manW + extraBits + resProdMoreThan2)
      val resShifted = resProd >> resShift
      assert(bit(manW+extraBits, resShifted) == 1)
      assert((resShifted >> (manW+extraBits+1)) == 0)

      val zmanTaylor = (resShifted >> extraBits) + bit(extraBits-1, resShifted)
      val zexTaylor = -(manW - (xmanbp-1)) + lntermMoreThan2 + resProdMoreThan2

      (zexTaylor.toInt, slice(0, manW, zmanTaylor.toLong))

    } else if(xexNobias == -1 && (((1L<<manW) - x.man.toLong) < (1L << (manW-taylorThreshold+1)))) {

      // log(1-x) = -x/ln(2) - x^2/2ln(2) - x^3/3ln(2) - O(x^4)
      //          = -x(1 + x/2 + x^2/3) * (1 / ln(2))

      // the exponent is -1, so xman is "multiplied" by 2. So the threshold and
      // exponent calculation become different
      val xman   = (1L<<manW) - x.man.toLong
      val xmanbp = xman.toBinaryString.length

      val invln2   = math.round((1.0 / log(2.0)) * (1 << (manW+extraBits))).toLong // > 1
      val oneThird = math.round((1.0 / 3.0)      * (1L << fracW)).toLong

      // 1 + x/2 > 1
      val onePlusHalfx = (1L << (manW+extraBits)) + (xman << (extraBits-2))

      // x^2/3
      val xsq      = xman * xman
      val xsqThird = ((xsq * oneThird) >> (xmanbp + xmanbp)) >> ((manW - xmanbp+1) * 2)

      val taylorTerm = onePlusHalfx + xsqThird

      val lntermProd      = xman * taylorTerm
      val lntermMoreThan2 = bit(xmanbp + manW + extraBits, lntermProd)
      val lnterm          = lntermProd >> (xmanbp - 1 + lntermMoreThan2)
      assert(bit(manW+extraBits, lnterm) == 1)
      assert(lnterm >> (manW+extraBits+1) == 0)

      val resProd = invln2 * lnterm // W = 2 * (manW+extraBits) + 2

      val resProdMoreThan2 = bit(2*(manW+extraBits)+1, resProd)
      val resShift   = (manW + extraBits + resProdMoreThan2)
      val resShifted = resProd >> resShift
      assert(bit(manW+extraBits, resShifted) == 1)
      assert(resShifted >> (manW+extraBits+1) == 0)

      val zmanTaylor = (resShifted >> extraBits) + bit(extraBits-1, resShifted)
      val zexTaylor = -(manW - (xmanbp-1)) + lntermMoreThan2 - 1 + resProdMoreThan2

      (zexTaylor.toInt, slice(0, manW, zmanTaylor.toLong))

    } else if(xexNobias == 0) {

      // --------------------------------------------------------------------------
      // polynomial (x is in [1, 2))
      //
      // if x == 0,
      // log2(x) = log2(2^ex * 1.man)
      //         = ex + log2(1.man)
      //         = log2(1.man)
      //
      val xman  = x.man.toLong
      val xex   = -(manW - xman.toBinaryString.length)

      val spAdrW  = tSmallPos.adrW
      val spFracW = tSmallPos.bp
      val spExtraBits = spFracW - manW

      val dxbp = manW - spAdrW - 1
      val d    = slice(0,      dxbp+1, xman) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, spAdrW, xman)
      val z0   = tSmallPos.interval(adr.toInt).eval(d.toLong, dxbp)

      val (zex0, zman0) = if(bit(spFracW-1, z0) == 1) {
        (xex  , z0 << 1) // the result from table < 1, but manW1 should be >1
      } else if (bit(spFracW-2, z0) == 1) {
        (xex-1, z0 << 2)
      } else {
        assert(bit(spFracW-3, z0) == 1)
        (xex-2, z0 << 3)
      }

      val zRound = (zman0 >> spExtraBits) + bit(spExtraBits-1, zman0)

      val zman   = slice(0, manW, zRound)
      val zex    = zex0 + bit(manW+1, zRound)

      (zex.toInt, slice(0, manW, zman.toLong))
    } else if(xexNobias == -1) {

      // polynomial (x is in [0.5, 1))

      val xman  = ((1L<<manW) - x.man.toLong) // 1-x
      val xex   = -(manW - xman.toBinaryString.length)

      val snAdrW  = tSmallNeg.adrW
      val snFracW = tSmallNeg.bp
      val snExtraBits = snFracW - manW

      val dxbp = manW - snAdrW - 1
      val d    = slice(0,      dxbp+1, xman) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, snAdrW, xman)

      val z0 = tSmallNeg.interval(adr.toInt).eval(d.toLong, dxbp)

      val (zex0, zman0) = if(bit(snFracW-1, z0) == 1) {
        (xex-1, z0 << 1) // the result from table < 1, but manW1 should be >1
      } else if (bit(snFracW-2, z0) == 1) {
        (xex-2, z0 << 2)
      } else {
        assert(bit(snFracW-3, z0) == 1)
        (xex-3, z0 << 3)
      }

      val zRound = (zman0 >> snExtraBits) + bit(snExtraBits-1, zman0)

      val zman   = slice(0, manW, zRound)
      val zex    = zex0 + bit(manW+1, zRound)

      (zex.toInt, slice(0, manW, zman.toLong))
    } else {
      // --------------------------------------------------------------------------
      // polynomial (0.0 < x < 0.5, 2.0 < x < inf)

      val dxbp = manW-adrW-1
      val d    = slice(0,      dxbp+1, x.man) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, adrW,   x.man)

      val zfrac0Pos = t.interval(adr.toInt).eval(d.toLong, dxbp)
      val zfrac0 = if(xexNobias >= 0) {zfrac0Pos} else {(1<<fracW) - zfrac0Pos}
      val zfrac  = zfrac0 & maskL(fracW)
      val zfull0 = (zint0 << fracW) + zfrac0.toLong

      assert(0L <= zfrac && zfrac < (1L<<fracW)) // avoid overflow in polynomial
      assert(0 <= zfull0)
      val zfullW  = zfull0.toBinaryString.length
      val zShiftW = exW + fracW - zfullW
      assert(0 <= zShiftW)

      val zShifted = zfull0 << zShiftW
      assert(bit(exW + fracW - 1, zShifted) == 1)

      val zman0       = slice(exW-1, fracW, zShifted) // -1 for the hidden bit
      val zmanRounded = slice(extraBits, manW, zman0) + bit(extraBits-1, zman0)
      val zman = if (zmanRounded<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
        0L
      } else if (zmanRounded >= (1L<<manW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
        maskL(manW)
      } else {
        zmanRounded.toLong
      }

      val zex  = (exW-1) - zShiftW

      (zex.toInt, slice(0, manW, zman.toLong))
    }

    // --------------------------------------------------------------------------
    // convert log2 to ln

    // 1/log2(e) < 1
    val oneOverLog2e = math.round(1.0 / log2(math.E) * (1 << (manW+1))).toLong
    assert((1<<manW) < oneOverLog2e && oneOverLog2e < (1<<(manW+1)))

    val zmanProd = ((1<<manW) + log2man0) * oneOverLog2e
    val zmanProdMoreThan2 = bit((manW+1)*2-1, zmanProd).toInt
    val zmanRound = slice(manW + zmanProdMoreThan2, manW, zmanProd) +
                    bit(manW - 1 + zmanProdMoreThan2, zmanProd)

    val zmanRoundMoreThan2 = bit(manW, zmanRound).toInt
    val zman = slice(0, manW, zmanRound)

    val zex = log2ex0 + exBias + zmanProdMoreThan2 + zmanRoundMoreThan2 - 1 // 1/log2(e) < 1

    return new RealGeneric(x.spec, zsgn, zex, SafeLong(zman))
  }
}
