//% @file reciprocalSim.scala
//
// Simulators for reciprocal function
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.table._
import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric

private[rial] object ReciprocalSim {

  def reciprocalSimGeneric( t : FuncTableInt, x: RealGeneric ) : RealGeneric = {

    val adrW = t.adrW
    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias
    val extraBits = t.bp - manW

    val man    = x.man
    val exRaw  = x.ex
    val sgn    = x.sgn
    val ex     = exRaw-exBias

    // Check NaN, Inf, Zero
    if (x.isNaN) return RealGeneric.nan(x.spec)
    if (x.isZero) return RealGeneric.inf(x.spec,0)
    if (x.isInfinite) return RealGeneric.zero(x.spec)

    val zEx = exBias - ex + ( if (man==0) 0 else -1 )
    if ( zEx >= maskI(expW) ) return RealGeneric.inf(x.spec,0)
    else if ( zEx <= 0 ) return RealGeneric.zero(x.spec)

    val order =
      if (adrW>=manW) {
        if (t.nOrder != 0)
          println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
        0
      } else {
        t.nOrder
      }
    val calcW = manW + extraBits

    val zman0 = if (order==0) {
      if (adrW<manW)
        println("WARNING: table address width < mantissa width, for polynomial order is zero. address width set to mantissa width.")
      val adr = man.toInt
      if (man==0) {0} else {t.interval(adr).eval(0L, 0)}
    } else {
      val dxbp = manW-adrW-1
      val d   = (SafeLong(1)<<dxbp) - slice(0, manW-adrW, man)-1
      val adr = maskI(adrW)-slice(manW-adrW, adrW, man).toInt
      //println(f"d=$d%x")
      // From here Long is used instead of SafeLong; should be fixed
      if (man==0) {0} else {t.interval(adr).eval(d.toLong, dxbp)}
    }

    val zman = if (zman0 < 0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
      0L
    } else if (zman0 >= (SafeLong(1)<<calcW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
      maskL(calcW)
    } else {
      zman0
    }

    // Simple rounding
    val z = if (extraBits>0) {(zman>>extraBits) + bit(extraBits-1, zman)} else {zman}

    new RealGeneric(x.spec, sgn, zEx, SafeLong(z))
  }

  // fracW include extra bits added during calc.
  def reciprocalTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    val tableD = if (order==0) {
      new FuncTableDouble( x => 2.0/(1.0+x)-1.0, order )
    } else {
      val eps=pow(2.0,-manW)
      new FuncTableDouble( x => 2.0/(2.0-(x+eps))-1.0, order )
    }
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting  )
  }

  val reciprocalF32TableI = ReciprocalSim.reciprocalTableGeneration( 2, 8, 23, 23+2 )
  val reciprocalF32Sim = reciprocalSimGeneric(reciprocalF32TableI, _ )

  val reciprocalBF16TableI = ReciprocalSim.reciprocalTableGeneration( 0, 7, 7, 7 )
  val reciprocalBF16Sim = reciprocalSimGeneric( reciprocalBF16TableI, _ )
}
