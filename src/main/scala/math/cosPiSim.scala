//% @file cosPiSim.scala
//
// Simulator for cosPi function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import chisel3.util.log2Up

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

// This only considers x is in [-1, 2), scaled by pi
// cos(pix): floating -> floating
//
// cos(pix) = 1 - pi^2/2 x^2 + pi^4/24 x^4
//          < 1 - 4.935  x^2   + 4.059 x^4
//          < 1 - 2^3    x^2   + 2^3   x^4
//
// if x < 2^-13 : cos(pix) = 1
// if x < 2^-7  : cos(pix) = 1 - (pi^2/2)x^2
// otherwise    : interpolation
//
// generic:
// ex*2 + 3 < -manW <=> ex < -ceil((manW + 3) / 2)
// ex*4 + 3 < -manW <=> ex < -ceil((manW + 3) / 4)
//
object CosPiSim {

  def cosPiSimGeneric( ts : Seq[FuncTableInt], x: RealGeneric ) : RealGeneric = {

    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias

    if (x.isNaN)      return RealGeneric.nan(x.spec)
    if (x.isInfinite) return RealGeneric.nan(x.spec)

    val xexNobias = x.ex - x.spec.exBias

    // skip large x (2 < x)
    if (1 <= xexNobias) {
      if (xexNobias == 1 && x.man == 0) {
        return new RealGeneric(x.spec, 1.0)
      } else {
        return RealGeneric.nan(x.spec)
      }
    }

    // here, x.ex <= 0

    //         y
    //         ^
    //         |-.--.      .--.
    // _.______|'_'._'.__.'__.'_ x
    // -1'.__.'|    '._''__.'
    //         |       1     2
    //           '-'---'-----'
    //            |  -1  ex=0
    //            -2

    // convert cos to sin
    val (zSgn, xex, xman) = if (xexNobias == 0) { // 1 ~ 2

      val (sgn, from0) = if (bit(manW-1, x.man) == 1) {
        (0, x.man - (1 << (manW-1)))
      } else {
        (1, (1 << (manW-1)) - x.man)
      }

      val shift = (manW+1) - from0.toLong.toBinaryString.length
      val norm  = from0 << shift

      val newex  = if(from0 == 0) {-exBias.toLong} else {-shift.toLong}
      val newman = if(from0 == 0) { 0     .toLong} else {(norm - (1 << manW)).toLong}

      (sgn, newex, newman)

    } else if (xexNobias == -1) { // 0.5 ~ 1

      val from0 = x.man
      val shift = (manW+1) - from0.toLong.toBinaryString.length
      val norm  = from0 << shift
      val newex  = if(from0 == 0) {-exBias.toLong} else {(-shift-1).toLong}
      val newman = if(from0 == 0) { 0     .toLong} else {(norm - (1 << manW)).toLong}

      (1, newex, newman)

    } else if (xexNobias >= -manW) { // 0 ~ 0.5

      val from0 = (1 << manW) - (((1<<manW) + x.man) >> (-1 - xexNobias))
      val shift = (manW+1) - from0.toLong.toBinaryString.length
      val norm  = from0 << shift

      val newex  = if(from0 == 0) {-exBias.toLong} else {(-shift-1).toLong}
      val newman = if(from0 == 0) { 0     .toLong} else {(norm - (1 << manW)).toLong}

      (0, newex, newman)
    } else {
      (0, -2.toLong, ((1<<manW)-1).toLong)
    }
//     println(f"x      = ${x.toDouble}")
//     println(f"xexNobias = ${xexNobias}")
//     println(f"xman      = ${x.man.toInt.toBinaryString}")
//     println(f"newx   = ${new RealGeneric(x.spec, 0, xex.toInt + exBias, xman).toDouble}")
//     println(f"newex  = ${xex}")
//     println(f"newman = ${xman.toInt.toBinaryString}")
    assert(xex  <= -1)
    if (xex == -1) {assert(xman == 0)}
    assert(xman <  (1<<manW))

    // {xex, xman} is in [0, 1/2].

    val linearThreshold = calcLinearThreshold(manW)
    val coef1 = new RealGeneric(x.spec, Pi)

    if (xex == -exBias && xman == 0) { // sin(0) = 0

      return RealGeneric.zero(x.spec)

    } else if (xex == -1 && xman == 0) { // sin(pi/2) = 1

      return new RealGeneric(x.spec, zSgn, exBias, 0)

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

      val adrW      = t.adrW
      val nOrder    = t.nOrder
      val bp        = t.bp
      val extraBits = bp - manW
      val calcW     = manW + extraBits

      val order = if(adrW >= manW) {
          if (nOrder != 0)
            println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
          0
        } else {
          nOrder
        }

      val (zEx, zman) = if (order == 0) {
        val adr   = xman.toInt
        val res0  = t.interval(adr).eval(0L, 0)
        val lessThanHalf = if(bit(calcW-1, res0) == 0) { 1 } else { 0 }
        val ex    = xex+2-lessThanHalf
        val res   = (res0 << (1+lessThanHalf)).toLong - (1 << calcW)

        (ex.toInt, res)

      } else {
        val dxbp = manW-adrW-1
        val d    = slice(0, manW-adrW, xman) - (SafeLong(1)<<dxbp)
        val adr  = slice(manW-adrW, adrW, xman).toInt

        val res = t.interval(adr).eval(d.toLong, dxbp)
        val lessThanHalf = if(bit(calcW-1, res) == 0) { 1 } else { 0 }
        ((xex+2-lessThanHalf).toInt, (res << (1+lessThanHalf)).toLong - (1L<<calcW))
      }
      val zmanRound = if (extraBits>0) {(zman>>extraBits) + bit(extraBits-1, zman)} else {zman}

      val z = if (zman<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
        0L
      } else if (zmanRound >= (1L<<manW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
        println(f"WARNING (${this.getClass.getName}) : x = ${x.toDouble}, sin(x) = ${sin(Pi * x.toDouble)}, z = ${new RealGeneric(x.spec, zSgn, zEx.toInt + exBias, SafeLong(maskL(manW))).toDouble}, zman = ${zmanRound.toBinaryString}")
        maskL(manW)
      } else {
        zmanRound
      }
      new RealGeneric(x.spec, zSgn, zEx.toInt + exBias, SafeLong(z))
    }
  }

  def calcLinearThreshold(manW: Int): Int = {
    // -12 for FP32
    math.floor(log2D(  6.0 * math.pow(2, -manW) / (Pi * Pi)) * 0.5).toInt
  }
  def calcCubicThreshold(manW: Int): Int = {
    // -6 for FP32
    math.floor(log2D(120.0 * math.pow(2, -manW) / math.pow(Pi, 4)) * 0.25).toInt
  }

  // number of tables depending on the exponent and linearThreshold
  def calcExAdrW(spec: RealSpec, allowCubicInterpolation: Boolean = false): Int = {
    // .--- table interp ---. .--------- cubic approx ---------.  .-- linear approx --.
    // -2 ~ cubicThreshold+1, cubicThreshold ~ linearThreshold+1, linearThreshold, -inf
    if (allowCubicInterpolation) {
      val cubicThreshold = calcLinearThreshold(spec.manW)
      val nTables = -cubicThreshold - 2 // = -2 - (cubicThreshold+1) + 1
      log2Up(nTables)
    } else {
      val linearThreshold = calcLinearThreshold(spec.manW)
      val nTables = -linearThreshold - 2 // = -2 - (linearThreshold+1) + 1
      log2Up(abs(linearThreshold)+1)
    }
  }

  def cosPiTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    val linearThreshold = calcLinearThreshold(manW)

    if(adrW >= manW) {assert(order == 0)}

    val maxCalcWidth = (-2 to linearThreshold by -1).map(exponent => {
      val tableD = new FuncTableDouble( x => scalb(sin(Pi * scalb(1.0+x, exponent)), -exponent-3), order )
      tableD.addRange(0.0, 1.0, 1<<adrW)
      val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
      tableI.calcWidth
    }).reduce( (lhs, rhs) => {
      lhs.zip(rhs).map( x => max(x._1, x._2))
    })

    (-2 to linearThreshold by -1).map( i => {
      val tableD = new FuncTableDouble( x => scalb(sin(Pi * scalb(1.0+x, i)), -i-3), order )
      tableD.addRange(0.0, 1.0, 1<<adrW)
      new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
    })
  }

  val cosPiF32TableI = CosPiSim.cosPiTableGeneration( 2, 8, 23, 23+2 )
  val cosPiF32Sim = cosPiSimGeneric(cosPiF32TableI, _ )

  val cosPiBF16TableI = CosPiSim.cosPiTableGeneration( 0, 7, 7, 7 )
  val cosPiBF16Sim = cosPiSimGeneric(cosPiBF16TableI, _ )
}
