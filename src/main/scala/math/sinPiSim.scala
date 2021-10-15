//% @file sinPiSim.scala
//
// Simulator for sinPi function
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
import rial.arith.Rounding._
import rial.arith._

// This only considers x is in [-2, 2), scaled by pi
// sin(pix): floating -> floating
//
// sin(pix) = pix - pi^3/3! x^3 + pi^5/5! x^5 - ...
//          < pix - 5.168 x^3   + 2.551 x^5
//          < pix - 8 x^3       + 4 x^5
//
// pi^3/3! x^3 < pix * 2^-23
// pi^2/3! x^2 < 2 x^2 < 2^-23
//         x^2 < 2^-24
//         x   < 2^-12
//
// pi^5x^5/5! < pix * 2^-23
// pi^4/5! x^4 < 2^-23
// pi^4/5! x^4 < x^4 < 2^-23
//               x   < 2^-6
//
// if x < 2^-12 : sin(pix) = pix
// if x < 2^-5 : 4x^3 < 2^-23, sin(pix) = pix - 8x^3
// otherwise   : interpolation
//
// generic:
// ex*3 + 3 < -manW <=> ex < -ceil((manW + 3) / 3)
// ex*5 + 2 < -manW <=> ex < -ceil((manW + 2) / 5)
//
object SinPiSim {

  def sinPiSimGeneric( ts : Seq[FuncTableInt], x: RealGeneric ) : RealGeneric = {

    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias

    if (x.isNaN)      return RealGeneric.nan(x.spec)
    if (x.isInfinite) return RealGeneric.nan(x.spec)

    val xexNobias = x.ex - x.spec.exBias

    // skip large x (2 < x)
    if (1 <= xexNobias) {
      if (xexNobias == 0 && x.man == 0) {
        return RealGeneric.zero(x.spec)
      } else {
        return RealGeneric.nan(x.spec)
      }
    }

    // here, x.ex <= 0

    //         y
    //         ^
    //         | .--.
    // _.______|'____'.______.__ x
    // -1'.__.'|      1'.__.' 2
    //         |
    //           '-'--'------'
    //            | -1  ex=0
    //            -2

    val zSgn = if ((x.sgn == 1) || (xexNobias == 0)) { 1 } else { 0 }

    // round everything into ex = -inf to -2

    val (xex, xman) = if (xexNobias == 0) { // 1 ~ 2
      val halfrange = if (bit(manW-1, x.man) == 1) {
        (1 << manW) - x.man
      } else {
        x.man
      }
      assert(halfrange >= 0)

      val shiftW = manW + 1 - halfrange.toLong.toBinaryString.length()
      val newex  = if(halfrange == 0) { -exBias } else { -shiftW }
      val newman = (halfrange << shiftW) - (1 << manW)

      (newex.toLong, newman.toLong)

    } else if (xexNobias == -1) { // 0.5 ~ 1
      val halfrange = (1 << manW) - x.man

      val shiftW = manW+1 - halfrange.toLong.toBinaryString.length()
      val newex  = -1-shiftW
      val newman = (halfrange << shiftW) - (1 << manW)

      (newex.toLong, newman.toLong)

    } else { // already small enough
      (xexNobias.toLong, x.man.toLong)
    }
    assert(xex  <= -2)
    assert(xman <  (1<<manW))
    // {xex, xman} is in [0, 1/2).

    val linearThreshold = (-math.ceil((manW + 1) / 2)).toLong // -12, if FP32
    val coef1 = new RealGeneric(x.spec, Pi)

    if (xex == -exBias && xman == 0) { // sin(0) = 0

      return RealGeneric.zero(x.spec)

    } else if (xex < linearThreshold) { // sin(pix) = pix

      val prodEx        = coef1.ex-exBias + xex
      val prodMan       = (coef1.man + (1<<manW)).toLong * (xman + (1<<manW)).toLong
      val prodbp        = manW + manW
      val prodMoreThan2 = bit(prodbp+1, prodMan)
      val prodRoundBits = prodbp - manW + prodMoreThan2
      val prodRound     = roundBySpec(RoundSpec.roundToEven, prodRoundBits.toInt, SafeLong(prodMan))
      val prodMoreThan2AfterRound = bit(manW + 1, prodRound)
      val prodExInc     = if(prodMoreThan2 == 1 || prodMoreThan2AfterRound == 1) {1} else {0}

      val zMan = if (prodMoreThan2AfterRound == 1) {prodRound >> 1} else {prodRound}
      val zEx  = prodEx + exBias + prodExInc
      return new RealGeneric(x.spec, zSgn, zEx.toInt, zMan)

    } else { // table interpolation

      // table interpolate for x in [0, 1/2) (x.ex = -2, -3, -4, -5 if FP32)
      val exadr = (-xex - 2).toInt
      val t = ts(exadr)

      val adrW   = t.adrW
      val nOrder = t.nOrder
      val bp     = t.bp
      val extraBits = if (nOrder==0) {0} else {(bp - manW)}
      val calcW  = manW + extraBits

      val dxbp = manW-adrW-1
      val d   = (SafeLong(1)<<dxbp) - slice(0, manW-adrW, xman) - 1
      val adr = maskI(adrW)-slice(manW-adrW, adrW, xman).toInt

      val res0 = t.interval(adr).eval(d.toLong, dxbp)
      val res  = (1 << calcW) - res0

      val shift = calcW+1 - res.toLong.toBinaryString.length
      val zEx   = -shift;

      val zman = (res << shift).toLong - (1 << calcW)
      val zmanRound = if (extraBits>0) { (zman>>extraBits) + bit(extraBits-1, zman)} else {zman}

      val z = if (zman<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
        0L
      } else if (zmanRound >= (1L<<manW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
        maskL(manW)
      } else {
        zmanRound
      }
      new RealGeneric(x.spec, zSgn, zEx.toInt + exBias, SafeLong(z))
    }
  }

  def sinPiTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    val linearThreshold = (-math.ceil((manW + 1) / 2)).toInt // -12, if FP32
//     val cubicThreshold  = (-math.ceil( manW / 4)).toInt

    // ex=-2 means the value is in [1/4, 1/2).
    (-2 to linearThreshold by -1).map( i => {
      val eps = pow(2.0, -manW)
      val tableD = new FuncTableDouble( x => 1.0 - sin(Pi * scalb(2.0 - (x+eps), i)), order )
      tableD.addRange(0.0, 1.0, 1<<adrW)
      new FuncTableInt( tableD, fracW )
    })
  }

  val sinPiF32TableI = SinPiSim.sinPiTableGeneration( 2, 8, 23, 23+2 )
  val sinPiF32Sim = sinPiSimGeneric(sinPiF32TableI, _ )

//   val sinPiBF16TableI = sinPiSim.sinPiTableGeneration( 0, 7, 7, 7 )
//   val sinPiBF16Sim = sinPiSimGeneric(sinPiBF16TableI, _ )
}
