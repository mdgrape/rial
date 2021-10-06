//% @file sqrtSim.scala
//
// Simulators for sqrt function
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

object SqrtSim {

  def sqrtSimGeneric( t_even : FuncTableInt, t_odd : FuncTableInt, x: RealGeneric ) : RealGeneric = {
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

    if (x.isNaN)      return RealGeneric.nan (x.spec)
    if (x.isZero)     return RealGeneric.zero(x.spec)
    if (x.isInfinite) return RealGeneric.inf (x.spec, sgn)
    if (sgn == 1)     return RealGeneric.nan (x.spec)

    val zSgn = 0
    val zEx  = (ex >> 1) + exBias

    val order =
      if (adrW>=manW) {
        if (nOrder != 0)
          println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
        0
      } else {
        nOrder
      }
    val calcW = manW + extraBits

    val zman0 = if (order==0) {
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
      // we convert sqrt(x) to -sqrt(-x). To do that, we invert x and also invert
      // the result. Also, [0, 1) -> (-1, 0], so we need to add 1 to make the range (0, 1]
      //

      // convert d into [-1, 1] in the same way as FuncIntervalDouble.eval -> evalNorm
//       val dxbp  = manW - adrW
//       val d     = (slice(0, dxbp, man) << 1) - SafeLong(1 << dxbp)
//          -invert-> SafeLong(1 << dxbp) - (slice(0, dxbp, man) << 1)
//          -  +1  -> SafeLong(1 << dxbp) - (slice(0, dxbp, man) << 1) + 1
//       val adr   = slice(manW-adrW, adrW, man).toInt

      val dxbp = manW-adrW-1
      // we subtract 1 to convert 0 -> 1111(-1), not 0 -> 0
      val d   = (SafeLong(1)<<dxbp) - slice(0, manW-adrW, man) - 1
      val adr = maskI(adrW)-slice(manW-adrW, adrW, man).toInt

//       println(f"d = ${d.toLong.toBinaryString}(${d.toLong})")

      val res0 = if(ex % 2 == 0) {
//         println(f"cw = ${t_even.interval(adr).cw(0)._1},  ${t_even.interval(adr).cw(1)._1},  ${t_even.interval(adr).cw(2)._1}")
        t_even.interval(adr).eval(d.toLong, dxbp)
      } else {
//         println(f"cw = ${t_odd .interval(adr).cw(0)._1},  ${t_odd .interval(adr).cw(1)._1},  ${t_odd .interval(adr).cw(2)._1}")
        t_odd .interval(adr).eval(d.toLong, dxbp)
      }
      // here we get y = 2 - sqrt(...)
      val rres = res0 - 2 * (SafeLong(1) << calcW) //   y-2  (-sqrt(...))
      val res = -rres                              // -(y-2) ( sqrt(...))
      res.toLong
    }
    // remove leading 1
    val zman = slice(0, calcW, zman0)

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

  //
  // 2^{2k} * 1.m -> 2^k * sqrt(1.m = 1.0 + x)
  //
  // To make all the 0th, 1st, 2nd coefficients positive...
  // sqrt(1+x) -(x:=-x)->  sqrt(1-x) -(x:=x+1)-> sqrt(2-x)
  //           -(y:=-y)-> -sqrt(2-x) -(y:=y+2)-> 2-sqrt(2-x)
  //
  def sqrtTableGenerationEven( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    val tableD = if (order == 0 || order == 1) {
      new FuncTableDouble( x => sqrt(1.0+x), order )
    } else {
      val eps=pow(2.0,-manW)
      new FuncTableDouble( x => (2.0 - sqrt(2.0-(x+eps))), order )
    }
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  //
  // 2^{2k+1} * 1.m -> 2^k * sqrt(2 * 1.m)
  //
  // To make all the 0th, 1st, 2nd coefficients positive...
  // sqrt(2*(1+x)) -(x:=-x)->  sqrt(2*(1-x)) -(x:=x+1)-> sqrt(2*(2-x))
  //               -(y:=-y)-> -sqrt(2*(2-x)) -(y:=y+2)-> 2-sqrt(2*(2-x))
  //
  def sqrtTableGenerationOdd( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    val tableD = if (order == 0 || order == 1) {
      new FuncTableDouble( x => sqrt(2.0 * (1.0 + x)), order )
    } else {
      val eps=pow(2.0,-manW)
      new FuncTableDouble( x => 2.0 - sqrt(2.0 * (2.0 - (x+eps))), order )
    }
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  val sqrtF32TableIEven = SqrtSim.sqrtTableGenerationEven( 2, 8, 23, 23+2 )
  val sqrtF32TableIOdd  = SqrtSim.sqrtTableGenerationOdd ( 2, 8, 23, 23+2 )
  val sqrtF32Sim = sqrtSimGeneric(sqrtF32TableIEven, sqrtF32TableIOdd, _ )

  val sqrtBF16TableIEven = SqrtSim.sqrtTableGenerationEven( 0, 7, 7, 7 )
  val sqrtBF16TableIOdd  = SqrtSim.sqrtTableGenerationOdd ( 0, 7, 7, 7 )
  val sqrtBF16Sim = sqrtSimGeneric(sqrtBF16TableIEven, sqrtBF16TableIOdd, _ )
}
