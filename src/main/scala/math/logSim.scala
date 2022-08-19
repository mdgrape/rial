//% @file logSim.scala
//
// Simulators for log(x)
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import chisel3.util.log2Up

import spire.math.SafeLong
import spire.math.Numeric
import spire.math.Real
import spire.implicits._

import rial.table._
import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.Rounding._
import rial.arith._

//
// case 1: x < 0.5 or 2 <= x. (ex + log2(1.man)) * ln2
// case 2: 0.5 <= x < 1. calc -log(x)/(1-x) * (1-x)
// case 3: 1.0 <= x < 2. calc  log(x)/(x-1) * (x-1)
//
// no taylor expansion needed.
//

object LogSim {

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

  // table for [1.0, 2.0).
  // calc log(x) / (x-1). it will be in [ln(2)~0.6, 1.0)
  def logSmallPositiveTableGeneration(spec: RealSpec,
      order:     Int =  2,
      adrW:      Int =  8,
      extraBits: Int =  2,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ): FuncTableInt = {

    val fracW = spec.manW + extraBits
    val f = (x0: Double) => {
      // log(x)/(x-1) -> 1.0 when x -> 1. without this, it becomes NaN.
      if(x0 == 0.0) {
        1.0
      } else {
        val x = x0 + 1.0 // 0~1 -> 1~2
        val res = log(x) / x0 // log(x) / (x-1)
        assert(0.0 <= res && res < 1)
        res
      }
    }
    val tableD = new FuncTableDouble( f, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    return new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }

  // table for [0.5, 1.0).
  // calc -log(x) / (1-x) - 1. -1 is needed to fit the result range into [0,1).
  def logSmallNegativeTableGeneration(spec: RealSpec,
      order:     Int =  2,
      adrW:      Int =  8,
      extraBits: Int =  2,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ): FuncTableInt = {

    val fracW = spec.manW + extraBits
    val f = (x0: Double) => {
      // -log(x)/(1-x) -> 1.0 when x -> 1. without this, it becomes NaN.
      if(x0 == 1.0) {
        0.0 // -log(x)/(1-x) - 1 -> 0 when x -> 1
      } else {
        val x = (1.0 + x0) * 0.5 // 0~1 -> 0.5~1
        val res = -log(x) / (1.0 - x) - 1.0 // -log(x) / (1-x) - 1.0
        assert(0.0 <= res && res < 1)
        res
      }
    }
    val tableD = new FuncTableDouble( f, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    return new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }

  //
  // log(2^ex * 1.man) = log(2^ex) + log(1.man)
  //                   = ex + log(1.man)
  //
  def logSimGeneric(
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

    val zsgn = if(x.ex < exBias) {1} else {0}

    // --------------------------------------------------------------------------

    val xexNobias = x.ex - exBias
    val zint0     = if(xexNobias >= 0) {SafeLong(xexNobias)} else {SafeLong(-xexNobias - 1)}

//     println(f"xexNobias = ${xexNobias}")
//     println(f"zint0     = ${zint0}")

    // All path calculates z = (table result) * (constant) = y * c.
    // calc y and c.
    // y is normalized as its mantissa width is fracW.
    val (ymanW1, zex0, cmanW1) = if(x.ex == exBias - 1) {
      // ----------------------------------------------------------------------
      // 0.5 <= x < 1.0.
      //
      // y = -log(x) / (1-x),
      // c = (1-x).
      //

      val dxbp = manW - adrW - 1
      val d    = slice(0,      dxbp+1, x.man) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, adrW, x.man)

      val res = tSmallNeg.interval(adr.toInt).eval(d.toLong, dxbp)
      val ymanW1 = if (res<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        (SafeLong(1) << fracW)
      } else if (res >= (SafeLong(1)<<fracW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        maskSL(fracW+1)
      } else {
        (SafeLong(1) << fracW) + SafeLong(res)
      }
      // result will be in [ln(2), 1.0). ln(2) > 0.5, so yex == -1.
      val yex0 = -1
      assert(bit(fracW, ymanW1) == 1, f"ymanW1 = ${ymanW1.toLong.toBinaryString}")

      // constant to be multiplied is 1-x.
      val c0 = (SafeLong(1) << (manW+1)) - x.manW1
      val c0Shift = manW+1 - binaryWidthSL(c0)
      val c0Shifted = c0 << c0Shift
      assert((SafeLong(1) << manW) <= c0Shifted && c0Shifted < (SafeLong(1) << (manW+1)))

      val zex0 = exBias + yex0 - c0Shift

      (ymanW1, zex0.toInt, c0Shifted)

    } else if(x.ex == exBias) {
      // --------------------------------------------------------------------------
      // 1.0 <= x < 2.0.
      //
      // y = log(x) / (x-1). table approximates log(x) / (x-1) - 1.
      // c = (x-1)
      //

      val dxbp = manW - adrW - 1
      val d    = slice(0,      dxbp+1, x.man) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, adrW,   x.man)

      val res = tSmallPos.interval(adr.toInt).eval(d.toLong, dxbp)
      val ymanW1 = if (res<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        SafeLong(0)
      } else if (res >= (SafeLong(1)<<fracW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        maskSL(fracW+1)
      } else {
        SafeLong(res) << 1 // shift to normalize (make fracW-th bit 1)
      }
      val yex0 = 0
      assert(bit(fracW, ymanW1) == 1, f"ymanW1 = ${ymanW1.toLong.toBinaryString}")

      // constant to be multiplied is x-1.
      val c0      = x.man
      val c0Shift = manW+1 - binaryWidthSL(c0)
      val c0Shifted = c0 << c0Shift
      assert((SafeLong(1) << manW) <= c0Shifted && c0Shifted < (SafeLong(1) << (manW+1)))

      val zex0 = exBias + yex0 - c0Shift - 1 // zres has been shifted

      (ymanW1, zex0.toInt, c0Shifted)

    } else {
      // --------------------------------------------------------------------------
      // polynomial (0.0 < x < 0.5, 2.0 < x < inf).

      val dxbp = manW-adrW-1
      val d    = slice(0,      dxbp+1, x.man) - (SafeLong(1) << dxbp)
      val adr  = slice(dxbp+1, adrW,   x.man)

      val zfrac0Pos0 = t.interval(adr.toInt).eval(d.toLong, dxbp)
      // overflow check
      val zfrac0Pos = if (zfrac0Pos0<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        SafeLong(0)
      } else if (zfrac0Pos0 >= (SafeLong(1)<<fracW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.value.toLong.toBinaryString}(${x.toDouble})")
        maskSL(fracW)
      } else {
        SafeLong(zfrac0Pos0)
      }

      val zfrac0 = if(xexNobias >= 0) {
        zfrac0Pos
      } else if(zfrac0Pos == 0) {
        maskSL(fracW)
      } else {
        (SafeLong(1)<<fracW) - zfrac0Pos
      }

      assert(0L <= zfrac0 && zfrac0 < (SafeLong(1)<<fracW))
      val zfrac  = zfrac0 & maskSL(fracW)
      val zfull0 = (zint0 << fracW) + zfrac

      assert(0L <= zfrac && zfrac < (SafeLong(1)<<fracW)) // avoid overflow in polynomial
      assert(0 <= zfull0, f"zfull0 = ${zfull0} < 0, zint = ${zint0}, zfrac = ${zfrac}")
      val zfullW  = binaryWidthSL(zfull0)
      val zShiftW = exW + fracW - zfullW
//       println(f"sim: zfullW = ${zfullW}")
//       println(f"sim: exW+fracW = ${exW+fracW}")
      assert(0 <= zShiftW)

      val zShifted = zfull0 << zShiftW
      assert(bit(exW + fracW - 1, zShifted) == 1)

      val yman0 = slice(exW-1, fracW, zShifted) // -1 for the hidden bit
      val ymanW1 = (SafeLong(1) << fracW) + yman0

      val zex0  = exBias + (exW-1) - zShiftW - 1 // ln2.ex == -1

      val cmanW1 = Real.log(Real.two)(manW+1)
      assert((SafeLong(1) << manW) < cmanW1 && cmanW1 < (SafeLong(1) << (manW+1)))

      (ymanW1, zex0.toInt, cmanW1)
    }

    // multiply

    val zmanProd = ymanW1 * cmanW1
    val zmanProdMoreThan2 = bit((fracW+1)+(manW+1)-1, zmanProd).toInt
    val zmanRound = slice(fracW + zmanProdMoreThan2, manW, zmanProd) +
                    bit(fracW + zmanProdMoreThan2 - 1, zmanProd)

    val zmanRoundMoreThan2 = bit(manW, zmanRound).toInt
    val lnman = slice(0, manW, zmanRound)
    val lnex = zex0 + zmanProdMoreThan2 + zmanRoundMoreThan2

    return new RealGeneric(x.spec, zsgn, lnex, lnman)
  }
}
