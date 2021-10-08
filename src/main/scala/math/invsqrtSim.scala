//% @file invsqrtSim.scala
//
// Simulators for 1/sqrt function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

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

object InvSqrtSim {

  def invsqrtSimGeneric( t_even : FuncTableInt, t_odd : FuncTableInt, x: RealGeneric ) : RealGeneric = {
    assert(t_even.nOrder == t_odd.nOrder)
    assert(t_even.adrW   == t_odd.adrW)
    assert(t_even.bp     == t_odd.bp)

    val adrW   = t_even.adrW
    val nOrder = t_even.nOrder
    val bp     = t_even.bp

    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias
    val extraBits = if (nOrder==0) {0} else {(bp - manW)}

    val man    = x.man
    val exRaw  = x.ex
    val sgn    = x.sgn
    val ex     = exRaw-exBias

    if (x.isNaN)      return RealGeneric.nan(x.spec)
    if (x.isZero)     return RealGeneric.inf(x.spec, sgn)
    if (x.isInfinite) return RealGeneric.zero(x.spec)
    if (sgn == 1)     return RealGeneric.nan(x.spec)

    // 1.m \in [1.0, 2.0)
    // =>   sqrt(1.m) \in [1.0, 1.414),     sqrt(2*1.m) \in [1.414, 2.0)
    // => 1/sqrt(1.m) \in [1/1.414, 1.0), 1/sqrt(2*1.m) \in [1/2.0, 1/1.414)
    // => 2/sqrt(1.m) \in [2/1.414, 2.0), 2/sqrt(2*1.m) \in [1.0, 2/1.414)
    //
    // 1/sqrt(2^ 5 * 1.m) = 2^-2 * 1/sqrt(2*1.m) = 2^-3 * 2/sqrt(2*1.m)
    // 1/sqrt(2^ 4 * 1.m) = 2^-2 * 1/sqrt(1.m)   = 2^-3 * 2/sqrt(1.m)
    // 1/sqrt(2^-4 * 1.m) = 2^ 2 * 1/sqrt(1.m)   = 2^ 1 * 2/sqrt(1.m)
    // 1/sqrt(2^-5 * 1.m) = 2^ 3 * 1/sqrt(2*1.m) = 2^ 2 * 2/sqrt(2*1.m)
    val zSgn = 0
    val zEx  = exBias - (ex >> 1) - 1

    val order =
      if (adrW>=manW) {
        if (nOrder != 0)
          println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
        0
      } else {
        nOrder
      }
    val calcW = manW + extraBits

    val zman = if (order==0) {
      if (adrW<manW) {
        println("WARNING: table address width < mantissa width, for polynomial order is zero. address width set to mantissa width.")
      }
      val adr = man.toInt
      if (ex % 2 == 0) {
        t_even.interval(adr).eval(0L, 0)
      } else {
        t_odd .interval(adr).eval(0L, 0)
      }
    } else {
      //      +- address part
      //      |        +- fraction part (from 00..0 to 11..1): dxbp
      //   .--+--. .---+---.
      // 1.xxxxxxx|xxxxxxxxx00
      //   <--------------->  : manW
      //   <----->            : adrW
      //           <------->  : manW - adrW
      //                    <>: extraBits
      //           <--------->: calcW

      // To make all the coefficients positive (to reduce complexity of multiplier),
      // we convert invsqrt(x) to sqrt(-x).

      val dxbp = manW-adrW-1
      // we subtract 1 to convert 0 -> 1111(-1), not 0 -> 0
      val d   = (SafeLong(1)<<dxbp) - slice(0, manW-adrW, man) - 1
      val adr = maskI(adrW)-slice(manW-adrW, adrW, man).toInt

      if(ex % 2 == 0) {
        t_even.interval(adr).eval(d.toLong, dxbp)
      } else {
        t_odd .interval(adr).eval(d.toLong, dxbp)
      }
    }
    // Simple rounding
    val zmanRound = if (extraBits>0) { (zman>>extraBits) + bit(extraBits-1, zman)} else zman
    //println(f"zman=${zman}%h")

    val z = if (zman<0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
      0L
    } else if (zmanRound >= (1L<<manW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
      maskL(manW)
    } else {
      zmanRound
    }
    new RealGeneric(x.spec, zSgn, zEx, SafeLong(z))
  }

  def invsqrtTableGenerationEven( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    val tableD = if (order == 0 || order == 1) {
      new FuncTableDouble( x => 2.0 / sqrt(1.0+x) - 1.0, order )
    } else {
      val eps=pow(2.0,-manW)
      new FuncTableDouble( x => 2.0 / sqrt(2.0-(x+eps)) - 1.0, order )
    }
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }
  def invsqrtTableGenerationOdd( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    val tableD = if (order == 0 || order == 1) {
      new FuncTableDouble( x => 2.0 / sqrt(2.0 * (1.0+x)) - 1.0, order )
    } else {
      val eps=pow(2.0,-manW)
      new FuncTableDouble( x => 2.0 / sqrt(2.0 * (2.0-(x+eps))) - 1.0, order )
    }
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  val invsqrtF32TableIEven = InvSqrtSim.invsqrtTableGenerationEven( 2, 8, 23, 23+2 )
  val invsqrtF32TableIOdd  = InvSqrtSim.invsqrtTableGenerationOdd ( 2, 8, 23, 23+2 )
  val invsqrtF32Sim = invsqrtSimGeneric(invsqrtF32TableIEven, invsqrtF32TableIOdd, _ )

  val invsqrtBF16TableIEven = InvSqrtSim.invsqrtTableGenerationEven( 0, 7, 7, 7 )
  val invsqrtBF16TableIOdd  = InvSqrtSim.invsqrtTableGenerationOdd ( 0, 7, 7, 7 )
  val invsqrtBF16Sim = invsqrtSimGeneric(invsqrtBF16TableIEven, invsqrtBF16TableIOdd, _ )
}
