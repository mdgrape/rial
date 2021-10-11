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

  // [1, 4) -> (1/2, 1]
  def invsqrtSimGeneric( t : FuncTableInt, x: RealGeneric ) : RealGeneric = {
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
    val emanW = manW + 2

    if (x.isNaN)      return RealGeneric.nan(x.spec)
    if (x.isZero)     return RealGeneric.inf(x.spec, sgn)
    if (x.isInfinite) return RealGeneric.zero(x.spec)
    if (sgn == 1)     return RealGeneric.nan(x.spec)

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
      t.interval(adr).eval(0L, 0) - (1 << calcW)
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

      val dxbp = emanW-adrW-1
      // we subtract 1 to convert 0 -> 1111(-1), not 0 -> 0
      val d   = (SafeLong(1)<<dxbp) - slice(0, emanW-adrW, man) - 1
      val adr = maskI(adrW)-slice(emanW-adrW, adrW, man).toInt

      t.interval(adr).eval(d.toLong, dxbp) - (1 << calcW)
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

  def invsqrtTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    if (order == 0) {
      val tableD = new FuncTableDouble( x => 2.0 / sqrt(1.0+x*4), order )
      tableD.addRange(0.0, 1.0, 1<<(manW+2))
      new FuncTableInt( tableD, fracW )
    } else {
      val eps=pow(2.0,-(manW+2))
      val tableD = new FuncTableDouble( x => 2.0 / sqrt(4.0-(x+eps)*4.0 + 1.0), order )
      tableD.addRange(0.0, 1.0, 1<<adrW)
      new FuncTableInt( tableD, fracW )
    }
  }

  val invsqrtF32TableI = InvSqrtSim.invsqrtTableGeneration( 2, 8, 23, 23+2 )
  val invsqrtF32Sim = invsqrtSimGeneric(invsqrtF32TableI, _ )

  val invsqrtBF16TableI = InvSqrtSim.invsqrtTableGeneration( 0, 7, 7, 7 )
  val invsqrtBF16Sim = invsqrtSimGeneric(invsqrtBF16TableI, _ )
}
