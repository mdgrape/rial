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

object ReciprocalSim {

  def reciprocalSimGeneric( t : FuncTableInt, x: RealGeneric ) : RealGeneric = {

    val adrW = t.adrW
    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias
    val extraBits = if (t.nOrder==0) 0 else (t.bp - manW)

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

    val zman = if (order==0) {
      if (adrW<manW) 
        println("WARNING: table address width < mantissa width, for polynomial order is zero. address width set to mantissa width.")
      val adr = man.toInt
      t.interval(adr).eval(0L, 0)
    } else {
      val dxbp = manW-adrW-1
      val d   = (SafeLong(1)<<dxbp) - slice(0, manW-adrW, man)
      val adr = slice(manW-adrW, adrW, man).toInt
      //println(f"d=$d%x")
      // From here Long is used instead of SafeLong; should be fixed
      t.interval(adr).eval(d.toLong, dxbp)
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
    new RealGeneric(x.spec, sgn, zEx, SafeLong(z))
  }

  val reciprocalF32ExtraBits = 2
  val eps = pow(2.0, -23)
  val reciprocalF32TableD = new FuncTableDouble( x => 2.0/(2.0-(x+eps))-1.0, 2 )
  reciprocalF32TableD.addRange(0.0, 1.0, 1<<8)
  val reciprocalF32TableI = new FuncTableInt( reciprocalF32TableD, 23+reciprocalF32ExtraBits )

  val reciprocalF32Sim = reciprocalSimGeneric(reciprocalF32TableI, _ )

  val reciprocalBF16TableD = new FuncTableDouble( x => 2.0/(1.0+x)-1.0, 0 )
  reciprocalBF16TableD.addRange(0.0, 1.0, 1<<7)
  val reciprocalBF16TableI = new FuncTableInt( reciprocalBF16TableD, 7 )
  val reciprocalBF16Sim = reciprocalSimGeneric( reciprocalBF16TableI, _ )
}
