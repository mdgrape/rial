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
    val manW1  = x.man + (1<<x.spec.manW)
    val man    = if (bit(0, ex) == 1) { manW1 >> 1 } else { manW1 >> 2 }

    if (x.isNaN)      {return RealGeneric.nan(x.spec)}
    if (x.isZero)     {return RealGeneric.inf(x.spec, sgn)}
    if (x.isInfinite) {return RealGeneric.zero(x.spec)}
    if (sgn == 1)     {return RealGeneric.inf(x.spec, sgn)}

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
      t.interval(adr).eval(0L, 0)
    } else {
      val dxbp = manW-adrW-1
      val d    = slice(0, manW-adrW, man) - (SafeLong(1)<<dxbp)
      val adr  = slice(manW-adrW, adrW, man).toInt
      val res  = t.interval(adr).eval(d.toLong, dxbp)
      res.toLong
    }
    // Simple rounding
    val zmanRound = if (extraBits>0) { (zman>>extraBits) + bit(extraBits-1, zman)} else zman
    //println(f"zman=${zman}%h")

    val z = if (zman<0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=${x.toDouble}")
      0L
    } else if (zmanRound >= (1L<<manW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=${x.toDouble}")
      maskL(manW)
    } else {
      zmanRound
    }
    new RealGeneric(x.spec, zSgn, zEx, SafeLong(z))
  }

  def invsqrtTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    val tableD = new FuncTableDouble( x => if(x<0.25) {0} else {2.0/sqrt(x*4) - 1.0}, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  val invsqrtF32TableI = InvSqrtSim.invsqrtTableGeneration( 2, 8, 23, 23+2 )
  val invsqrtF32Sim = invsqrtSimGeneric(invsqrtF32TableI, _ )

  val invsqrtBF16TableI = InvSqrtSim.invsqrtTableGeneration( 0, 7, 7, 7 )
  val invsqrtBF16Sim = invsqrtSimGeneric(invsqrtBF16TableI, _ )
}
