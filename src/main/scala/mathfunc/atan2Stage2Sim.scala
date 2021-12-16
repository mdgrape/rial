//% @file atan2Stage2Sim.scala
//
// Simulators for atan2(y, x) stage 2, calculating atan2 from min(x,y)/max(x,y)
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

import rial.math.ATan2Sim

object ATan2Stage2Sim {

  def atan2Stage2SimGeneric(
    ts : Seq[FuncTableInt], x : RealGeneric, status: Int, special: Int, ysgn: Int
    ): RealGeneric = {

//     println("==================================================")

    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val order     = ts(0).nOrder
    val adrW      = ts(0).adrW
    val fracW     = ts(0).bp
    val extraBits = fracW - manW

//     println(f"x = ${x.sgn}|${x.ex}(${xex})|${xman.toLong.toBinaryString}")

    // ------------------------------------------------------------------------
    // check special values

    val ysgnUnit = if(ysgn == 0) {1} else {-1}
    if(special == 1) {
      return RealGeneric.nan(x.spec)
    } else if (special == 2) { // zero
      return new RealGeneric(x.spec, ysgnUnit * 0)
    } else if (special == 3) { // pi
      return new RealGeneric(x.spec, ysgnUnit * Pi)
    } else if (special == 4) { // pi/2
      return new RealGeneric(x.spec, ysgnUnit * Pi * 0.5)
    } else if (special == 5) { // pi/4
      return new RealGeneric(x.spec, ysgnUnit * Pi * 0.25)
    } else if (special == 6) { // 3pi/4
      return new RealGeneric(x.spec, ysgnUnit * Pi * 0.75)
    }

    // ------------------------------------------------------------------------
    // atan table
    val linearThreshold = ATan2Sim.calcLinearThreshold(manW)

    val (atanEx, atanManW1) = if(x.ex - exBias < linearThreshold) {
      // linear approx
      (x.ex, x.manW1.toLong)
    } else {
      // table interpolation
      val exadr = (-(x.ex - exBias) - 1).toInt
      val t     = ts(exadr)

      val (zEx0, zMan0) = if (t.nOrder == 0) {

        val adr   = x.man.toInt
        val res0  = t.interval(adr).eval(0L, 0)

        val scaling = (-(x.ex - exBias) - 1).toInt
        val shift   = fracW+1 - res0.toLong.toBinaryString.length
        val res     = res0 << shift

        (-shift - scaling, res.toLong)

      } else {
        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, x.man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, x.man).toInt
        val res0 = t.interval(adr).eval(d.toLong, dxbp)

        // the result is multiplied by 2^(-ex-1).
        // if x < 1, atan(x) = x - x^3/3 + x^5/5 + O(x^7) is less than x. So
        // in case of x is in [2^-8, 2^-7), atan(x) < 2^-7. So we can multiply
        // atan(x) by 2^7 = 2^(-(-8)-1) to contain more bits. By doing this, we
        // can keep the same precision in all 2^-8 ~ 2^0 regions.
        //   We need to correct the exponent term by this scaling term.

        val scaling = (-(x.ex-exBias) - 1).toInt
        val shift   = fracW+1 - res0.toLong.toBinaryString.length
        val res     = res0 << shift // normalize

        (-shift - scaling, res.toLong)
      }

      // round the result from table
      val zMan = if (extraBits > 0) {
        (zMan0 >> extraBits) + bit(extraBits-1, zMan0)
      } else {
        zMan0
      }

      (zEx0, zMan.toLong)
    }

    // z correction by:
    //   ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|
    //   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y|
    //   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y|
    //   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y|

    val zsgn = ysgn

    if (status == 0) {//   ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|

      return new RealGeneric(x.spec, zsgn, atanEx, atanManW1 - (1<<manW))

    } else if (status == 1) {//   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y|

      val halfpi = new RealGeneric(x.spec, Pi * 0.5)
      val halfpi_manW1sft3 = (halfpi.man + (1 << manW)) << 3
      val atan_shift       = -atanEx
      val atan_aligned     = (atanManW1 << 3) >> atan_shift

      //  ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y|
      val zman0 = halfpi_manW1sft3 - atan_aligned
      // z is in [0.78 ~ 1.57)

      val zman0LessThan1 = if(bit(manW+3, zman0) == 0) { 1 } else { 0 }
      val zmanRound = roundBySpec(RoundSpec.round, 3-zman0LessThan1, zman0)
      val zman = zmanRound - (1<<manW)
      val zex  = exBias - zman0LessThan1

      return new RealGeneric(x.spec, zsgn, zex, zman)

    } else if (status == 2) {//   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y|

      val pi = new RealGeneric(x.spec, Pi)

      // pi.ex = 1
      val pi_manW1   = (pi.man + (1 << manW)) << 3 // pi = 2^1 * 1.57...
      val atan_shift = -atanEx
      val atan_man   = (atanManW1 << 2) >> atan_shift
      val zman0      = pi_manW1 - atan_man // pi ~ pi/2.
      val zman0LessThan2 = if(bit(manW+3, zman0) == 0) { 1 } else { 0 }
      val zman = roundBySpec(RoundSpec.round, 3-zman0LessThan2, zman0)

      return new RealGeneric(x.spec, zsgn, 1-zman0LessThan2 + exBias, zman - (1<<manW))

    } else                  {//   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y|
      assert(status == 3)

      val halfpi = new RealGeneric(x.spec, Pi * 0.5)
      val halfpi_manW1sft3 = (halfpi.man + (1 << manW)) << 3
      val atan_shift       = -atanEx
      val atan_aligned     = (atanManW1 << 3) >> atan_shift

      val zman0 = halfpi_manW1sft3 + atan_aligned
      // z is in [1.57 ~ 2.35), ex is 0 or 1
      val zman0MoreThan2 = bit(manW+1+3, zman0)
      val zmanRound = roundBySpec(RoundSpec.round, 3+zman0MoreThan2, zman0)
      val zman = zmanRound - (1<<manW)
      val zex  = zman0MoreThan2 + exBias

      return new RealGeneric(x.spec, zsgn, zex, zman)
    }
  }
}
