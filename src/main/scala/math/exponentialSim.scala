//% @file exponential.scala
//
// Simulators for Exponential functions
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

//
// First we write fixed precision version, then generalize it
//
// Pow2_f32 : Power of 2 function, 32-bit single input / output
//            x = (-1)^s * 2^e * 1.m
//            -> xi = fixed point representation of x
//               overflow     if xi >=  0x80
//               underflow(0) if xi <= -0x7f
//               e >= 7, s = 0 -> overflow
//               e >= 7, s = 1 -> 0
//            z = 2^xi = 2^(xi.int) * 2^(xi.frac)
//               1 <= 2^(xi.frac) < 2

// Stage : Insert pipeline stages
//         Pipeline stages are inserted after the logics
//         as shift registers. If multiple registers are
//         inserted, they should be balanced in synthesis.

object ExponentialSim {

  def pow2simGeneric( t : FuncTableInt, extraBits: Int, x: RealGeneric ) : RealGeneric = {

    val adrW = t.adrW
    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias

    val man    = x.man
    val exRaw  = x.ex
    val sgn    = x.sgn
    val ex     = exRaw-exBias
    val overflow_value = RealGeneric.inf(x.spec,0)
    // Check NaN
    if (x.isNaN) return RealGeneric.nan(x.spec)

    if ( ex >= expW-1 ) { // Overflow or Underflow
      return ( if (sgn==0) overflow_value else RealGeneric.zero(x.spec) )
    } 
    val calcW0 = manW+extraBits
    val extraManW = if (calcW0<adrW) {
      if (t.nOrder != 0) {
        sys.error(f"ERROR: Address width $adrW > calculation width $calcW0 for polynomial order ${t.nOrder}")
        0
      } else {
        adrW-manW
      }
    } else {
      extraBits
    }
    val calcW = manW + extraManW

    val manWith1 = x.manW1<<(expW-2+extraManW)
    // right shift amount :
    //   0  if ex==expW-2
    val shift = expW-2-ex
    val xu = manWith1 >> shift
    //println(f"$xu%08x, $manWith1%08x, ", expW-2-ex)
    // binary point at manW+extraBits
    val xi = if (sgn==0) xu else -xu
    //println(f"xisim=$xi%x")
    // Result exponent
    val zEx = exBias + (xi>>calcW).toInt
    //println(f"ex=$ex%d, xu=$xu%08x, zEx=$zEx%d, $tmp%d")
    if   (zEx>=maskI(expW)) { return overflow_value }
    else if (zEx<=0) { return RealGeneric.zero(x.spec) }

    val dxbp = calcW-adrW-1
    val d    = slice(0, dxbp+1, xi) - (SafeLong(1)<<dxbp)
    //println(f"d=$d%x")
    val adr  = slice(dxbp+1, adrW, xi).toInt
    // From here Long is used instead of SafeLong; should be fixed
    val zman = t.interval(adr).eval(d.toLong, dxbp)
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
    new RealGeneric(x.spec, 0, zEx, SafeLong(z))
  }

  val pow2F32Order = 2
  val pow2F32AdrW  = 8
  val pow2F32ExtraBits = 2

  val pow2F32TableD = new FuncTableDouble( x => pow(2.0,x)-1.0, pow2F32Order )
  pow2F32TableD.addRange(0.0, 1.0, 1<<pow2F32AdrW)
  val pow2F32TableI = new FuncTableInt( pow2F32TableD, 23+pow2F32ExtraBits )

  val pow2F32Sim = pow2simGeneric(pow2F32TableI, pow2F32ExtraBits, _)

  val pow2BF16Order = 0
  val pow2BF16AdrW  = 7
  val pow2BF16ExtraBits = 1

  val pow2BF16TableD = new FuncTableDouble( x => pow(2.0,x)-1.0, pow2BF16Order )
  pow2BF16TableD.addRange(0.0, 1.0, 1<<pow2BF16AdrW)
  val pow2BF16TableI = new FuncTableInt( pow2BF16TableD, 7+pow2BF16ExtraBits )

  val pow2BF16Sim = pow2simGeneric(pow2BF16TableI, pow2BF16ExtraBits, _)

  def expSimGeneric( t : FuncTableInt, extraBits: Int, x: RealGeneric ) : RealGeneric = {
    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias

    val man    = x.manW1
    val xe     = x.ex
    val sgn    = x.sgn

    if (x.isNaN || x.isInfinite || x.isZero) {
      pow2simGeneric(t, extraBits, x )
    } else {
      val log2e = 1.0d/java.lang.Math.log(2.0d)
      val log2eMan = BigInt(round(java.lang.Math.scalb(log2e, manW)))
      val prod = man * log2eMan
      val moreThan2 = bit(manW*2+1, prod)
      val ym = slice(manW+moreThan2, manW, prod)
      val r  = bit(manW-1+moreThan2, prod)
      val stickey = if (slice(0, manW-2+moreThan2, prod) != 0) 1 else 0 
      val inc = r & (stickey | bit(0,ym))
      val ymRound = ym + inc
      val ye = if ((bit(manW,ymRound)|moreThan2) !=0 ) xe+1 else xe
      val yman = ymRound & maskSL(manW)
      val y = new RealGeneric(x.spec, sgn, ye, yman)
      //println(f"${x.value}%h ${y.value}%h")
      val z = pow2simGeneric(t, extraBits, y)
      //println(f"${z.value}%h ${x.toDouble}%f ${y.toDouble}%f ${z.toDouble}%f")
      z
    }
  }

}
