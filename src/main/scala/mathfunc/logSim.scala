//% @file logSim.scala
//
// Simulators for log2 and log(x)
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

  def calcTaylorThreshold(spec: RealSpec): Int = {
    // log(1+x) = 1/ln(2) * (x - x^2/2 + x^3/3 - x^4/4 + O(x^5))
    //          = x(1 - x/2 + x^2/3 - x^3/4) / ln(2)
    // cond: x^3/4          < 2^-manW
    //   <=> x^3            < 2^-manW+2
    //   <=> x              < 2^(-manW+2)/3
    //   <=> 2^x.ex * 1.man < 2^(-manW+2)/3
    //   <=> 2^x.ex         < 2^((-manW+2)/3) / 2
    //   <=> x.ex           < (-manW+2)/3 - 1
    -floor((2-spec.manW)/3 - 1).toInt
    // log(1-x) = 1/ln(2) * (-x - x^2/2 - x^3/3 - x^4/4 + O(x^5))
    //          = -x(1 + x/2 + x^2/3 + x^3/4) / ln(2)
    // the condition is the same
  }

  def logNormalTableGeneration(spec: RealSpec,
    order:     Int =  2,
    adrW:      Int =  8,
    extraBits: Int =  2,
    calcWidthSetting: Option[Seq[Int]] = None,
    cbitSetting: Option[Seq[Int]] = None
    ): FuncTableInt = {

    val log2 = (a:Double) => {log(a) / log(2.0)}
    val fracW = spec.manW + extraBits

    val tableNormalD = new FuncTableDouble( x => log2(1.0 + x), order)
    tableNormalD.addRange(0.0, 1.0, 1<<adrW)
    return new FuncTableInt(tableNormalD, fracW, calcWidthSetting, cbitSetting)
  }

  // XXX:  these special tables has extremely large c2 width.
  // TODO: reduce the bitwidth by removing taylor part from this table range

  // table for [1.0, 2.0)
  def logSmallPositiveTableGeneration(spec: RealSpec,
      order:     Int =  2,
      adrW:      Int =  8,
      extraBits: Int =  2,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ): FuncTableInt = {

    val log2 = (a:Double) => {log(a) / log(2.0)}
    val fracW = spec.manW + extraBits

    // log(1+x) = 1/ln(2)[x - x^2/2 + x^3/3 + O(x^4)]
    //          = 1.44.. [x - x^2/2 + x^3/3 + O(x^4)]
    //          < 2x
    // x < log2(1+x), 0<=x<=1. (log2 is convex upward, log2(1+0)=0, log2(1+1)=1)
    // so
    // x < log2(1+x) < 2x
    //
    // 2^xex * 1.xman < log2(1+x)            < 2^(xex+1) * 1.xman
    // 2^xex          < log2(1+x)            < 2^(xex+1) * 2
    // 1              < log2(1+x) * 2^-xex   < 2^2
    // 0.25           < log2(1+x) * 2^-xex-2 < 1

    val taylorThreshold = -calcTaylorThreshold(spec)

    val f = (x: Double) => {
      val f64      = RealSpec.Float64Spec
      val xi       = java.lang.Double.doubleToRawLongBits(x)
      val xex      = slice(f64.manW, f64.exW, xi) - f64.exBias.toLong
      val res      = if(xex < taylorThreshold) {
        0.0
      } else {
        val baseline = pow(2.0, (-xex-2).toDouble)
        val z        = log2(1.0 + x)
        z * baseline
      }
      res
    }
    val tableD = new FuncTableDouble( f, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    return new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }

  // table for [0.5, 1.0)
  def logSmallNegativeTableGeneration(spec: RealSpec,
      order:     Int =  2,
      adrW:      Int =  8,
      extraBits: Int =  2,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ): FuncTableInt = {

    val log2 = (a:Double) => {log(a) / log(2.0)}
    val fracW = spec.manW + extraBits

    val taylorThreshold = -calcTaylorThreshold(spec)

    val f = (x: Double) => {
      val f64      = RealSpec.Float64Spec
      val xi       = java.lang.Double.doubleToRawLongBits(x)
      val xex      = slice(f64.manW, f64.exW, xi) - f64.exBias.toLong
      val res      = if(xex < taylorThreshold) {
        0.0
      } else {
        val baseline = pow(2.0, (-xex-1).toDouble)
        val z        = -log2(1.0 - x*0.5)
        z * baseline
      }
      res
    }

    val tableD = new FuncTableDouble( f, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)

    return new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }


  //
  // log(2^ex * 1.man) = log(2^ex) + log(1.man)
  //                    = ex + log(1.man)
  //
  def logSimGeneric( islog2: Boolean,
    t : FuncTableInt, tSmallPos : FuncTableInt, tSmallNeg : FuncTableInt,
    x : RealGeneric ): RealGeneric = {

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
    if(x.man == 0 && x.ex == exBias) {
      return RealGeneric.zero(x.spec)
    }
    // XXX table for x in (0.5, 1] fails if x == 0.5
    if(x.man == 0 && x.ex == exBias-1) {
      if(islog2) {
        return new RealGeneric(x.spec, 1, exBias, 0)
      } else {
        // to reproduce rounding
        val log2ex0  = exBias
        val log2man0 = 0L
        val oneOverLog2e = math.round(log(2.0) * (1 << (fracW+1))).toLong
        assert((1<<fracW) < oneOverLog2e && oneOverLog2e < (1<<(fracW+1)))
        val zmanProd = ((1<<fracW) + log2man0) * oneOverLog2e
        val zmanProdMoreThan2 = bit((fracW+1)*2-1, zmanProd).toInt
        val zmanRound = slice(fracW + extraBits + zmanProdMoreThan2, manW, zmanProd) +
                        bit(fracW + extraBits - 1 + zmanProdMoreThan2, zmanProd)
        val zmanRoundMoreThan2 = bit(manW, zmanRound).toInt
        val zman = slice(0, manW, zmanRound)
        val zex = log2ex0 + zmanProdMoreThan2 + zmanRoundMoreThan2 - 1 // 1/log2(e) < 1

        return new RealGeneric(x.spec, 1, zex, zman)
      }
    }

    val zsgn = if(x.ex < exBias) {1} else {0}

    // --------------------------------------------------------------------------

    val xexNobias = x.ex.toLong - exBias.toLong
    val zint0     = if(xexNobias >= 0) {xexNobias} else {-xexNobias - 1}

//     println(f"xexNobias = ${xexNobias}")
//     println(f"zint0     = ${zint0}")

    val taylorThreshold = calcTaylorThreshold(x.spec)

    val (log2ex0, log2man0) = if(xexNobias == 0 && x.man.toLong < (1L << (manW-taylorThreshold))) {
      // ----------------------------------------------------------------------
      // taylor
      //
      // log(1+x) = x/ln(2) - x^2/2ln(2) + x^3/3ln(2) + O(x^4)
      //          = x(1 - x/2 + x^2/3) / ln(2)
      //

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

      // 1 - x/2 + x^2/3 < 1, x < 2^-8
      val taylorTerm = oneMinusHalfx + xsqThird // W = fracW
      // x < x/ln2 ~ x * 1.44 < 2x
      val convTerm   = invln2 * xman            // W = fracW + xmanbp

      // x < x/ln2 * (1 - x/2 + x^2/3) < 2x
      val resProd      = convTerm * taylorTerm  // W = 2*fracW + xmanbp
      val resMoreThan2 = bit(xmanbp + fracW*2, resProd)
      val resShifted   = (resProd >> (fracW+xmanbp-1 + resMoreThan2)) +
                         bit(fracW+xmanbp-1 + resMoreThan2-1, resProd)
      val resShiftedMoreThan2 = bit(fracW+1, resShifted) // resShifted includes 1 at 2^0.

      val zmanTaylor = slice(0, fracW, resShifted)
      val zexTaylor = -(manW - (xmanbp-1)) + resMoreThan2 + resShiftedMoreThan2

      (zexTaylor.toInt + exBias, zmanTaylor)

    } else if(xexNobias == -1 && (((1L<<manW) - x.man.toLong) < (1L << (manW-taylorThreshold+1)))) {

      // --------------------------------------------------------------------------
      //
      // log(1-x) = -x/ln(2) - x^2/2ln(2) - x^3/3ln(2) - O(x^4)
      //          = -x(1 + x/2 + x^2/3) * (1 / ln(2))
      //
      // the exponent is -1, so xman is "multiplied" by 2. So the threshold and
      // exponent calculation become different
      val xman   = ((1L<<manW) - x.man.toLong).toBigInt
      val xmanbp = xman.toLong.toBinaryString.length

      val invln2   = math.round((1.0 / log(2.0)) * (1L << fracW)).toBigInt // > 1
      val oneThird = math.round((1.0 / 3.0)      * (1L << fracW)).toBigInt

      // 1 + x/2 > 1
      val onePlusHalfx = (1L << fracW).toBigInt + (xman << (extraBits-2))

      // x^2/3
      val xsq      = xman * xman
      val xsqThird = ((xsq * oneThird) >> (xmanbp + xmanbp)) >> ((manW - (xmanbp-1)) * 2)

      // 1 + x/2 + x^2/3 > 1
      val taylorTerm = onePlusHalfx + xsqThird
      val convTerm   = xman * invln2

      val resProd = convTerm * taylorTerm
      val resMoreThan2 = bit(xmanbp + fracW + fracW, resProd)
      val resShifted   = (resProd >> (xmanbp - 1 + fracW + resMoreThan2)) +
                         bit(xmanbp - 1 + fracW + resMoreThan2 - 1, resProd)
                         // resShifted include 1 at 2^0.
      val resMoreThan2AfterRound = bit(fracW+1, resShifted)

      val zmanTaylor = slice(0, fracW, resShifted)
      val zexTaylor = -(manW - (xmanbp-1)) - 1 + resMoreThan2 + resMoreThan2AfterRound

      (zexTaylor.toInt + exBias, SafeLong(zmanTaylor))

    } else if(xexNobias == 0) {

      // --------------------------------------------------------------------------
      // polynomial (x is in [1, 2))
      //
      // if x == 0,
      // log(x) = log(2^ex * 1.man)
      //         = ex + log(1.man)
      //         = log(1.man)
      //
      // log table should return full precision
      //
      val xman  = x.man.toLong
      val xex   = -(manW - xman.toBinaryString.length)

      assert(tSmallPos.adrW == t.adrW)
      assert(tSmallPos.bp   == t.bp)

      val dxbp = manW - adrW - 1
      val d    = slice(0,      dxbp+1, xman) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, adrW, xman)

      val z00  = tSmallPos.interval(adr.toInt).eval(d.toLong, dxbp)
      val z0 = if (z00<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        0L
      } else if (z00 >= (1L<<fracW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        maskL(fracW)
      } else {
        z00.toLong
      }

      val (zex0, zman0) = if(bit(fracW-1, z0) == 1) {
        (xex   + exBias, z0 << 1) // the result from table < 1, but manW1 should be >1
      } else if (bit(fracW-2, z0) == 1) {
        (xex-1 + exBias, z0 << 2)
      } else {
        assert(bit(fracW-3, z0) == 1)
        (xex-2 + exBias, z0 << 3)
      }
      assert(bit(fracW, zman0) == 1)
      assert(zman0 < (1<<(fracW+1))) // zman0.getWidth == 1+fracW

      val zman   = slice(0, fracW, zman0)
      val zex    = zex0 + bit(fracW+1, zman0)

      (zex.toInt, SafeLong(zman))

    } else if(xexNobias == -1) {

      // --------------------------------------------------------------------------
      // polynomial (x is in [0.5, 1))
      //
      val xman  = ((1L<<manW) - x.man.toLong) // 1-x
      val xex   = -(manW - xman.toBinaryString.length)

      assert(tSmallNeg.adrW  == t.adrW)
      assert(tSmallNeg.bp    == t.bp)

      val dxbp = manW - adrW - 1
      val d    = slice(0,      dxbp+1, xman) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, adrW, xman)

      val z00 = tSmallNeg.interval(adr.toInt).eval(d.toLong, dxbp)
      val z0 = if (z00<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        0L
      } else if (z00 >= (1L<<fracW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        maskL(fracW)
      } else {
        z00.toLong
      }

      val (zex0, zman0) = if(bit(fracW-1, z0) == 1) {
        (xex-1 + exBias, z0 << 1) // the result from table < 1, but manW1 should be >1
      } else if (bit(fracW-2, z0) == 1) {
        (xex-2 + exBias, z0 << 2)
      } else {
        assert(bit(fracW-3, z0) == 1)
        (xex-3 + exBias, z0 << 3)
      }
      assert(bit(fracW, zman0) == 1)
      assert(zman0 < (1<<(fracW+1))) // zman0.getWidth == 1+fracW

      val zman   = slice(0, fracW, zman0)
      val zex    = zex0 + bit(fracW+1, zman0)

      (zex.toInt, SafeLong(zman))
    } else {
      // --------------------------------------------------------------------------
      // polynomial (0.0 < x < 0.5, 2.0 < x < inf)

      val dxbp = manW-adrW-1
      val d    = slice(0,      dxbp+1, x.man) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, adrW,   x.man)

      val zfrac0Pos0 = t.interval(adr.toInt).eval(d.toLong, dxbp)
      // overflow check
      val zfrac0Pos = if (zfrac0Pos0<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        0L
      } else if (zfrac0Pos0 >= (1L<<fracW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        maskL(fracW)
      } else {
        zfrac0Pos0.toLong
      }

      val zfrac0 = if(xexNobias >= 0) {zfrac0Pos} else {
        if(zfrac0Pos == 0) {maskL(fracW)} else {(1<<fracW) - zfrac0Pos}
      }
      assert(0L <= zfrac0 && zfrac0 < (1L<<fracW))
      val zfrac  = zfrac0 & maskL(fracW)
      val zfull0 = (zint0 << fracW) + zfrac.toLong

//       println(f"sim: zfrac0Pos0 = ${zfrac0Pos0.toLong.toBinaryString}")
//       println(f"sim: zfrac0Pos  = ${zfrac0Pos .toLong.toBinaryString}")
//       println(f"sim: zint0  = ${zint0.toLong.toBinaryString}")
//       println(f"sim: zfrac0 = ${zfrac0.toLong.toBinaryString}")
//       println(f"sim: zfull0 = ${zfull0.toLong.toBinaryString}")

      assert(0L <= zfrac && zfrac < (1L<<fracW)) // avoid overflow in polynomial
      assert(0 <= zfull0)
      val zfullW  = zfull0.toBinaryString.length
      val zShiftW = exW + fracW - zfullW
//       println(f"sim: zfullW = ${zfullW}")
//       println(f"sim: exW+fracW = ${exW+fracW}")
      assert(0 <= zShiftW)

      val zShifted = zfull0 << zShiftW
      assert(bit(exW + fracW - 1, zShifted) == 1)

      val zman0 = slice(exW-1, fracW, zShifted) // -1 for the hidden bit
      val zex0  = exBias + (exW-1) - zShiftW

      (zex0.toInt, SafeLong(zman0))
    }

//     println(f"sim: log2xEx0  = ${log2ex0 .toLong.toBinaryString}")
//     println(f"sim: log2xMan0 = ${log2man0.toLong.toBinaryString}")

    if(islog2) {
      //
      // now, log2man0 has `fracW` prec. round it to manW
      //
      val zmanRounded = slice(extraBits, manW, log2man0) + bit(extraBits-1, log2man0)
      val zmanMoreThan2 = bit(manW, zmanRounded)
      val zman = slice(0, manW, zmanRounded)
      val zex  = log2ex0 + zmanMoreThan2

      new RealGeneric(x.spec, zsgn, zex, zman)
    } else {
      // --------------------------------------------------------------------------
      // convert log2 to ln

      // 1/log2(e) < 1
      val oneOverLog2e = math.round(log(2.0) * (1 << (fracW+1))).toLong
      assert((1<<fracW) < oneOverLog2e && oneOverLog2e < (1<<(fracW+1)))

  //     println(f"sim: oneOverLog2e = ${oneOverLog2e.toLong.toBinaryString}")

      val zmanProd = ((1<<fracW) + log2man0) * oneOverLog2e
      val zmanProdMoreThan2 = bit((fracW+1)*2-1, zmanProd).toInt
      val zmanRound = slice(fracW + extraBits + zmanProdMoreThan2, manW, zmanProd) +
                      bit(fracW + extraBits - 1 + zmanProdMoreThan2, zmanProd)

      val zmanRoundMoreThan2 = bit(manW, zmanRound).toInt
      val zman = slice(0, manW, zmanRound)
      val zex = log2ex0 + zmanProdMoreThan2 + zmanRoundMoreThan2 - 1 // 1/log2(e) < 1

      return new RealGeneric(x.spec, zsgn, zex, zman)
    }
  }
}
