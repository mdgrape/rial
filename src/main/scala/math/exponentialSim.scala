//% @file exponential.scala
//
// Simulators for Exponential functions
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import rial.table._
import rial.util._
import rial.util.ScalaUtil._

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

  def pow2simGeneric( t : FuncTableInt, expW : Int, manW : Int, adrW: Int, extraBits: Int, x: Long ) : Long = {
    val man = slice(x, 0, manW)
    val exRaw  = slice(x, manW, expW)
    val sgn = bit(x, manW+expW)
    val exBias = mask(expW-1)
    val ex = exRaw-exBias
    val overflow_value = mask(expW)<<manW
    // Check NaN
    if ( (exRaw == mask(expW)) && (man != 0L) )  {
      return mask(expW+1) << (manW-1)
    }
    if ( ex >= expW-1 ) { // Overflow or Underflow
      return ( if (sgn==0) overflow_value else 0L )
    } 
    val manWith1 = (man+(1L<<manW))<<(expW-2+extraBits)
    // right shift amount :
    //   0  if ex==expW-2
    val shift = expW-2-ex
    val xu = if (shift>=64) 0L else (manWith1 >> (expW-2-ex))
    //println(f"$xu%08x, $manWith1%08x, ", expW-2-ex)
    // binary point at manW+extraBits
    val calcW = manW+extraBits
    val xi = if (sgn==0) xu else -xu

    // Result exponent
    val tmp = (xi>>calcW)
    val zEx = exBias + (xi>>calcW)
    //println(f"ex=$ex%d, xu=$xu%08x, zEx=$zEx%d, $tmp%d")
    if   (zEx>=mask(expW)) { return overflow_value }
    else if (zEx<=0) { return 0L }

    val dxbp = calcW-adrW-1
    val d    = slice(xi, 0, dxbp+1) - (1L<<dxbp)
    val adr  = slice(xi, dxbp+1, adrW).toInt
    val zman = t.interval(adr).eval(d, dxbp)
    // Simple rounding
    val zmanRound = (zman>>extraBits) + bit(zman, extraBits-1)
    val z = if (zman<0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
      0L
    } else if (zmanRound >= (1L<<manW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
      mask(manW)
    } else {
      zmanRound
    }
    (zEx << manW) + z
  }

  val pow2F32Order = 2
  val pow2F32AdrW  = 8
  val pow2F32ExtraBits = 2

  val pow2F32TableD = new FuncTableDouble( x => pow(2.0,x)-1.0, pow2F32Order )
  pow2F32TableD.addRange(0.0, 1.0, 1<<pow2F32AdrW)
  val pow2F32TableI = new FuncTableInt( pow2F32TableD, 23+pow2F32ExtraBits )

  val pow2F32Sim = pow2simGeneric(pow2F32TableI, 8, 23, pow2F32AdrW, pow2F32ExtraBits, _)
}
