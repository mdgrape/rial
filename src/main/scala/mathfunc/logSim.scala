//% @file logSim.scala
//
// Simulators for log(x)
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.mathfunc

import scala.language.reflectiveCalls
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
import rial.arith.Rounding._
import rial.arith._

object MathFuncLogSim {

  //
  // log_e(x) = log_2(x) / log_2(e)
  //
  def logSimGeneric( t : FuncTableInt, x : RealGeneric ): RealGeneric = {

//     println("==================================================")
    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val adrW      = t.adrW
    val fracW     = t.bp
    val extraBits = fracW - manW

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

    // --------------------------------------------------------------------------

    val xexNobias = x.ex.toLong - exBias.toLong
    val zint0     = if(xexNobias >= 0) {xexNobias} else {-xexNobias - 1}

//     println(f"xexNobias = ${xexNobias}")
//     println(f"zint0     = ${zint0}")

    // --------------------------------------------------------------------------
    // polynomial

    val dxbp = manW-adrW-1
    val d    = slice(0,      dxbp+1, x.man) - (SafeLong(1) << dxbp)
    val adr  = slice(dxbp+1, adrW,   x.man)

    // --------------------------------------------------------------------------
    // FP reconstruction

    val log2 = (a:Double) => {log(a) / log(2.0)}
    val oneOverLog2e = math.round(1.0 / log2(math.E) * (1 << manW)).toLong

    val zfrac0Pos = t.interval(adr.toInt).eval(d.toLong, dxbp)
    val zfrac0 = if(xexNobias >= 0) {zfrac0Pos} else {(1<<fracW) - zfrac0Pos}
//     println(f"s: zfrac0Pos = ${zfrac0Pos.toBinaryString}")
//     println(f"s: zfrac0    = ${zfrac0.toBinaryString   }")

    assert(0L <= zfrac0 && zfrac0 < (1L<<fracW))

    val zfull0 = (zint0 << fracW) + zfrac0.toLong
//     println(f"s: zfull0  = ${zfull0.toBinaryString }")
    assert(0 <= zfull0)

    val zfullConv = zfull0 * oneOverLog2e
    val zfull = zfullConv >> manW

    val zfullW  = zfull.toBinaryString.length
    val zShiftW = exW + fracW - zfullW
//     println(f"s: zfullW  = ${zfullW }")
//     println(f"s: zShiftW = ${zShiftW}")
    assert(0 <= zShiftW)

    val zShifted = zfull << zShiftW
//     println(f"s: zShifted = ${zShifted.toBinaryString }")
    assert(bit(exW + fracW - 1, zShifted) == 1)

    val zman0       = slice(exW-1, fracW, zShifted) // -1 for the hidden bit
    val zmanRounded = slice(extraBits, manW, zman0) + bit(extraBits-1, zman0)
//     println(f"s: zman0 = ${slice(extraBits, manW, zman0).toBinaryString } + ${bit(extraBits-1, zman0)}")
//     println(f"s: zman0 = ${zman0.toBinaryString }")
//     println(f"s: zmanR = ${zmanRounded.toBinaryString }")

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

    val zsgn = if(x.ex < exBias) {1} else {0}

    new RealGeneric(x.spec, zsgn, zex, SafeLong(zman))
  }
}
