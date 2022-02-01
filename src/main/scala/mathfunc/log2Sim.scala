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

  def calcTaylorThreshold(spec: RealSpec): Int = {
    // log(1+x) = x/ln(2) - x^2/2ln(2) + x^3/3ln(2) + O(x^4)
    //          = x(1 - x/2 + x^2/3) / ln(2)
    // cond: x^2/3 < 2^-manW
    //   <=> x^2   < 2^-manW+1
    //   <=> x     < 2^(-manW+1)/2
    //   <=> x.ex  < (-manW+1)/2
    -floor((1-spec.manW)/2).toInt
  }

  // table for [1.0, 2.0)
  def log2SmallPositiveTableGeneration(spec: RealSpec,
      order:     Int =  2,
      adrW:      Int = 11,
      extraBits: Int =  4,
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

    val f = (x: Double) => {
      val f64      = RealSpec.Float64Spec
      val xi       = java.lang.Double.doubleToRawLongBits(x)
      val xex      = slice(f64.manW, f64.exW, xi) - f64.exBias.toLong
      val baseline = pow(2.0, -xex-2)
      val z        = log2(1.0 + x)
      z * baseline
    }
    val tableD = new FuncTableDouble( f, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    return new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }

  // table for [0.5, 1.0)
  def log2SmallNegativeTableGeneration(spec: RealSpec,
      order:     Int =  2,
      adrW:      Int =  8,
      extraBits: Int =  2,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ): FuncTableInt = {

    val log2 = (a:Double) => {log(a) / log(2.0)}
    val resW = spec.exW + spec.manW
    val fracW = resW + extraBits

    val f = (x: Double) => {
      val z   = log2((1.0 + x) * 0.5)
      val zFP = (new RealGeneric(spec, z)).value
      val res = slice(0, resW, zFP).toDouble * pow(2.0, -resW)
      res
    }

    val tableD = new FuncTableDouble( f, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)

    return new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }


  //
  // log2(2^ex * 1.man) = log2(2^ex) + log2(1.man)
  //                    = ex + log2(1.man)
  //
  def log2SimGeneric( t : FuncTableInt, tSmallPos : FuncTableInt, tSmallNeg : FuncTableInt, x : RealGeneric ): RealGeneric = {

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
    if(x.ex == exBias && x.man == 0) {
      return RealGeneric.zero(x.spec)
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
    //          = x(1 - x/2 + x^2/3) / ln(2)
    //
    // if x < 2^-12, then x^2 < 2^-24
    // use x(1 - x/2)/ln(2)
    //    -x(1 + x/2)/ln(2)

    val taylorThreshold = calcTaylorThreshold(x.spec)

    if(xexNobias == 0 && x.man.toLong < (1L << (manW-taylorThreshold))) {
      val xman = x.man.toLong // x - 1
      val xmanbp = xman.toBinaryString.length

      val invln2 = math.round((1.0 / log(2.0)) * (1 << (manW+extraBits))).toLong // > 1

      val oneMinusHalfx = (1L << (manW+extraBits)) - (xman << (extraBits-1)) // < 1

      val lntermProd      = xman * oneMinusHalfx
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

      return new RealGeneric(x.spec, zsgn, zexTaylor.toInt + exBias, zmanTaylor)
    }

    // the exponent is -1, so xman is "multiplied" by 2. So the threshold and
    // exponent calculation become different
    if(xexNobias == -1 && (((1L<<manW) - x.man.toLong) < (1L << (manW-taylorThreshold+1)))) {
      val xman   = (1L<<manW) - x.man.toLong
      val xmanbp = xman.toBinaryString.length

      val invln2 = math.round((1.0 / log(2.0)) * (1 << (manW+extraBits))).toLong // > 1

      val onePlusHalfx = (1L << (manW+extraBits)) + (xman << (extraBits-2))

      val lntermProd      = xman * onePlusHalfx
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

      return new RealGeneric(x.spec, zsgn, zexTaylor.toInt + exBias, zmanTaylor)
    }

    // --------------------------------------------------------------------------
    // polynomial (x is in [1, 2))
    //
    // if x == 0,
    // log2(x) = log2(2^ex * 1.man)
    //         = ex + log2(1.man)
    //         = log2(1.man)
    //
    // log2 table should return full precision
    //
    if(xexNobias == 0) {
      val xman  = x.man.toLong
      val xex   = -(manW - xman.toBinaryString.length)

      val spAdrW  = tSmallPos.adrW
      val spFracW = tSmallPos.bp
      val spExtraBits = spFracW - manW

      val dxbp = manW - spAdrW - 1
      val d    = slice(0,      dxbp+1, xman) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, spAdrW, xman)
//       println(f"xex  = ${x.ex.toLong}")
//       println(f"xman = ${xman.toLong.toBinaryString}%32s")
//       println(f"dx   = ${slice(0,      dxbp+1, xman).toLong.toBinaryString}%32s")
//       println(f"d    = ${(if(d < 0) {"-" + (-d).toLong.toBinaryString} else {d.toLong.toBinaryString})}%32s")
//       println(f"adr  = ${adr .toLong.toBinaryString}%32s")

      val z0     = tSmallPos.interval(adr.toInt).eval(d.toLong, dxbp)

//       println(f"           3         2         1")
//       println(f"         21098765432109876543210987654321")
//       println(f"z0     = ${z0.toLong.toBinaryString}%32s")

      val (zex0, zman0) = if(bit(spFracW-1, z0) == 1) {
        (xex   + exBias, z0 << 1) // the result from table < 1, but manW1 should be >1
      } else if (bit(spFracW-2, z0) == 1) {
        (xex-1 + exBias, z0 << 2)
      } else {
        assert(bit(spFracW-3, z0) == 1)
        (xex-2 + exBias, z0 << 3)
      }
//       println(f"zman0  = ${zman0.toLong.toBinaryString}%32s")
//       println(f"dex    = ${zex0 - xex}")

      val zRound = (zman0 >> spExtraBits) + bit(spExtraBits-1, zman0)
//       println(f"zRound = ${zRound.toLong.toBinaryString}%32s")

      val zman   = slice(0, manW, zRound)
      val zex    = zex0 + bit(manW+1, zRound)

      return new RealGeneric(x.spec, zsgn, zex.toInt, zman.toLong)
    }

    // --------------------------------------------------------------------------
    // polynomial (x is in [0.5, 1))
    //
    // if x < 1,
    // log2(x) = log2(2^ex * 1.man)
    //         = ex + log2(1.man)
    //         = -(-ex - log2(1.man))
    //         = -(|ex| - 1 + 1 - log2(1.man))
    //             ^^^^^^^^   ^^^^^^^^^^^^^^^
    //             zint       zfrac
    // if ex == -1,
    // log2(x) = log2(2^ex * 1.man)
    //         = -1 + log2(1.man) // cancellation!
    //
    // log2(0.5 * 1.man)
    //
    // table only takes mantissa. We need to subtract 1 from it, anyway...
    //
//     if(xexNobias == -1) {
//       val xman0  = (1L<<manW) - x.man.toLong // 1-x
//       val xmanbp = xman0.toBinaryString.length
//       val xex    = manW - xmanbp // this does not consider the x.ex
//       val mask   = maskL(xmanbp - 1)
//       val xman   = xman0 & mask
// //       println( "         432109876543210987654321")
// //       println(f"xman   = ${x.man.toLong.toBinaryString}%24s")
// //       println(f"xman0  = ${xman0.toLong.toBinaryString}%24s")
// //       println(f"xex    = ${xex}")
// //       println(f"mask   = ${mask .toLong.toBinaryString}%24s")
// //       println(f"xman   = ${xman .toLong.toBinaryString}%24s")
//
//       val tableNegI = tSmallNeg(xex)
//       val smallAdrW = adrW
//
//       val xmanShifted = xman << xex
//
//       val dxbp = manW-smallAdrW-1-1
//       val d    = slice(0,      dxbp+1, xmanShifted) - (SafeLong(1) << dxbp)
//       val adr  = slice(dxbp+1, smallAdrW, xmanShifted)
//
// //       println(f"d      = ${slice(0,      dxbp+1, xmanShifted).toLong.toBinaryString}%24s")
// //       println(f"adr    = ${adr  .toLong.toBinaryString}%24s")
//
//       val zfrac0 = tableNegI.interval(adr.toInt).eval(d.toLong, dxbp)
//
// //       println(f"x in table = ${1.0 - ((1.0+xman.toDouble/(mask+1)) * pow(2.0, -xex-2))}")
// //       println(f"f in table = ${-log2(1.0 - ((1.0+xman.toDouble/(mask+1)) * pow(2.0, -xex-2)))}")
// //       println(f"zfrac0 = ${zfrac0.toLong.toBinaryString}")
// //       println(f"zfrac  = ${zfrac0.toDouble / (1<<fracW)}")
//
//       val (zmanW10, zex0) = if(bit(fracW-1, zfrac0) == 1) {
//         (zfrac0 << 1, -xex - 1)
//       } else if(bit(fracW-2, zfrac0) == 1) {
//         (zfrac0 << 2, -xex - 2)
//       } else {
//         (zfrac0 << 3, -xex - 3)
//       }
//
//       val zmanW1 = (zmanW10 >> extraBits) + bit(extraBits-1, zmanW10)
//       val zMoreThan2 = bit(manW+1, zmanW1)
//       val zman      = if(zMoreThan2 == 1) {
//         slice(1, manW, zmanW1)
//       } else {
//         slice(0, manW, zmanW1)
//       }
//       val zex = zex0 + zMoreThan2
//
//       return new RealGeneric(x.spec, zsgn, zex.toInt + exBias, zman)
//     }

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
