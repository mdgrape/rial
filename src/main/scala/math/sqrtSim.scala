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

  // [1, 4) -> [1, 2)
  def sqrtSimGeneric( t : FuncTableInt, x: RealGeneric ) : RealGeneric = {
//     println("-----------------------------")
    val adrW   = t.adrW
    val nOrder = t.nOrder
    val bp     = t.bp

    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias
    val extraBits = if (nOrder==0) {0} else {(bp - manW)}

    val exRaw  = x.ex
    val sgn    = x.sgn
    val ex     = exRaw-exBias
    val man    = if (bit(0, ex) == 1) {
      (x.man << 1) + (1 << x.spec.manW)
    } else {
      x.man
    }

//     println(f"x    = ${sgn}|${x.ex}(${ex})|${x.man.toLong.toBinaryString}")
//     println(f"xman = ${man.toLong.toBinaryString}")

    // [1, 2)        [2, 4)
    // 1.xxx -(*2)-> 1x.xx
    //  |             |
    //  v -1          v -10
    // 0.xxx          x.xx
    // [0, 1)        [0, 2) needs +1 -> [1, 3)

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

    val zman = if (order==0) {
      if (adrW<manW) {
        println("WARNING: table address width < mantissa width, for polynomial order is zero. address width set to mantissa width.")
      }
      val adr = man.toInt
      t.interval(adr).eval(0L, 0)
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

      val dxbp = (manW+2)-adrW-1
      val d   = (SafeLong(1)<<dxbp) - slice(0, (manW+2)-adrW, man) - 1
      val adr = maskI(adrW)-slice((manW+2)-adrW, adrW, man).toInt

//       println(f"dxbp = ${dxbp}, manW = ${manW}, adrW = ${adrW}")
//       println(f"d    = ${d.toLong.toBinaryString}, d/dmax = ${d.toDouble/(1<<dxbp)}")
//       println(f"adr  = ${adr.toBinaryString}, d/adrmax = ${adr.toDouble/(1<<adrW)}")

      val res0 = t.interval(adr).eval(d.toLong, dxbp) //   y    (2 - sqrt(...))
      val rres = res0 - (SafeLong(1) << (calcW+1))    //   y-2  (   -sqrt(...))
      val res = -rres                                 // -(y-2) (    sqrt(...))
//       println(f"res0  = ${res0.toLong.toBinaryString}")
//       println(f"rres  = ${rres.toLong.toBinaryString}")
//       println(f"res   = ${res .toLong.toBinaryString}")
      res.toLong
    }
//     println(f"zman      = ${zman.toBinaryString}(${zman}), extrabits = ${extraBits}")

    // Simple rounding
    val zmanRound = if (extraBits>0) { (zman>>extraBits) + bit(extraBits-1, zman)} else zman
//     println(f"zmanRound = ${zmanRound.toLong.toBinaryString}")

    val z = if (zman<0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
      0L
    } else if (zmanRound >= (1L<<manW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
      maskL(manW)
    } else {
      zmanRound
    }
//     println(f"z     = ${z.toLong.toBinaryString}")
    new RealGeneric(x.spec, zSgn, zEx, SafeLong(z))
  }

  def sqrtTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    if (order == 0) {
      // even: 2^2k * 1.m, odd: 2^2k+1 * 1.m = 2^2k * (2 * 1.m)
      // 1.m [= 1.0 + m] is in [1,2), 2*1.m [= 2.0 + 2 * m] is in [2,4)
      // m_even is in [0, 1), m_odd is in [0, 2). to align them, we need to add
      // 1.0 to m_odd to make the range [1, 3). That means that the width of
      // mantissa will increase twice. the first time is 2* in the odd case,
      // the seocnd is +1.0 in the odd case.
      val tableD = new FuncTableDouble( x => sqrt(1.0+x*4.0) - 1.0, order )
      tableD.addRange(0.0, 1.0, 1<<(manW+2))
      new FuncTableInt( tableD, fracW )
    } else {
      val eps=pow(2.0,-(manW+2)) // shift + 1, the range is [0, 4), not [0, 1)
      // To make all the 0th, 1st, 2nd coefficients positive...
      val tableD = new FuncTableDouble( x => 3.0 - sqrt(5.0-4.0*(x+eps)), order )
      tableD.addRange(0.0, 1.0, 1<<adrW)
      new FuncTableInt( tableD, fracW )
    }
  }

  val sqrtF32TableI = SqrtSim.sqrtTableGeneration( 2, 8, 23, 23+2 )
  val sqrtF32Sim = sqrtSimGeneric(sqrtF32TableI, _ )

  val sqrtBF16TableI = SqrtSim.sqrtTableGeneration( 0, 7, 7, 7 )
  val sqrtBF16Sim = sqrtSimGeneric(sqrtBF16TableI, _ )
}
