//% @file sinPiSim.scala
//
// Simulator for sinPi function Copyright (C) Toru Niina RIKEN BDR 2021
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
// sin(pi*x): floating -> floating

object SinPiSim {

  def sinPiSimGeneric(
    ts: Seq[FuncTableInt],
    x:  RealGeneric,
    useCubicTerm: Boolean = false
  ) : RealGeneric = {

    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias

    if (x.isNaN)      return RealGeneric.nan(x.spec)
    if (x.isInfinite) return RealGeneric.nan(x.spec)

    val xexNobias = x.ex - x.spec.exBias

    // skip large x (2 <= x)
    if (1 <= xexNobias) {
      if (xexNobias == 1 && x.man == 0) { // x = 2.0
        return RealGeneric.zero(x.spec)
      } else {
        return RealGeneric.nan(x.spec)
      }
    }

    if(0 <= xexNobias && x.man == 0) {
      return RealGeneric.zero(x.spec)
    }

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

    assert(xex  <= -1)
    assert(xman <  (1<<manW))
    if (xex == -1) {assert(xman == 0)}

    // {xex, xman} is in [0, 1/2].

    val linearThreshold = calcLinearThreshold(manW).toLong // -12
    val cubicThreshold  = calcCubicThreshold(manW).toLong  // -6
    val pi = new RealGeneric(x.spec, Pi)

    if (xex == -exBias && xman == 0) { // sin(0) = 0

      return RealGeneric.zero(x.spec)

    } else if (xex == -1 && xman == 0) { // sin(pi/2) = 1

      return new RealGeneric(x.spec, zSgn, exBias, 0)

    } else if (useCubicTerm && linearThreshold < xex && xex < cubicThreshold) {
      //
      // sin(pi*x) = pi*x - pi^3 * x^3 / 3!
      //           = pi*x * (1 - pi^2 * x^2 / 6)
      //
      val coef1 = pi
      val coef3 = new RealGeneric(x.spec, Pi * Pi / 6.0)
      assert(xex < cubicThreshold)

      // TODO
      // pi^2 / 6 ~ 1.5, and x < cubicThreshold, we can reduce the width of term2

      val xmanW1  = xman + (1<<manW)
      val xsqMan0 = (xmanW1 * xmanW1)
      val (xsqMan, xsqEx)  = if(bit(manW + manW + 1, xsqMan0) == 1L) {
        ((xsqMan0 >> (manW+1)) + bit(manW, xsqMan0), xex + xex + 1)
      } else {
        ((xsqMan0 >>  manW) + bit(manW-1, xsqMan0), xex + xex)
      }
      assert(xsqEx < cubicThreshold * 2)

      val term2Prod = coef3.manW1 * xsqMan
      val (term2Man, term2Ex) = if(bit(manW+manW+1, term2Prod) == 1L) {
        ((term2Prod >> (manW+1)) + bit(manW, term2Prod), xsqEx + coef3.exNorm + 1)
      } else {
        ((term2Prod >>  manW) + bit(manW-1, term2Prod),  xsqEx + coef3.exNorm)
      }
      assert(term2Ex <= cubicThreshold*2)

      val term2Shifted = term2Man >> (abs(term2Ex)-1).toInt
      val term2SubMan0 = (1 << (manW+1)) - term2Shifted
      val (term2SubMan, term2SubEx) = if(term2Shifted == 0) {
        ((1L<<manW), 0)
      } else {
        (term2SubMan0.toLong, -1)
      }

      val pixMan0 = pi.manW1.toLong * (xman + (1<<manW)).toLong
      val (pixMan, pixEx) = if(bit(manW+manW+1, pixMan0) == 1L) {
        ((pixMan0 >> (manW+1)) + bit(manW, pixMan0), pi.ex - exBias + xex + 1)
      } else {
        ((pixMan0 >> manW) + bit(manW-1, pixMan0), pi.ex - exBias + xex)
      }

      val cubicApproxMan0 = term2SubMan * pixMan
      val (cubicApproxMan, cubicApproxEx) = if(bit(manW+manW+1, cubicApproxMan0) == 1L) {
        ((cubicApproxMan0 >> (manW+1)) + bit(manW, cubicApproxMan0), pixEx + term2SubEx + 1)
      } else {
        ((cubicApproxMan0 >> manW) + bit(manW-1, cubicApproxMan0),   pixEx + term2SubEx)
      }

      return new RealGeneric(x.spec, zSgn, cubicApproxEx.toInt + x.spec.exBias, cubicApproxMan)

    } else if (xex < linearThreshold) { // sin(pix) = pix

      val prodEx        = pi.ex-exBias + xex
      val prodMan       = (pi.man + (1<<manW)).toLong * (xman + (1<<manW)).toLong
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

  // sin(pi x) = pi x - pi^3 x^3 / 3! + pi^5 x^5 / 5! + O(x^7)
  def calcLinearThreshold(manW: Int): Int = {
    // -12 for FP32
    // pi * x * 2^-manW > (pi * x)^3 / 3!
    math.floor(log2D(  6.0 * math.pow(2, -manW) / (Pi * Pi)) * 0.5).toInt
  }
  def calcCubicThreshold(manW: Int): Int = {
    // -6 for FP32
    //     pi * x * 2^-manW   > (pi * x)^5 / 5!
    // <=> 5! * 2^-manW       > (pi * x)^4
    // <=> 5! * 2^-manW /pi^4 > x^4
    // <=> log2(5! * 2^-manW / pi^4) / 4 > log2(x)
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

  def sinPiTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    val linearThreshold = calcLinearThreshold(manW)

    if(adrW >= manW) {assert(order == 0)}

    val maxCalcWidth = (-2 to linearThreshold by -1).map(exponent => {
        val tableD = new FuncTableDouble( x => scalb(sin(Pi * scalb(1.0 + x, exponent)), -exponent-3), order )
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

  val sinPiF32TableI = SinPiSim.sinPiTableGeneration( 2, 8, 23, 23+2 )
  val sinPiF32Sim = sinPiSimGeneric(sinPiF32TableI, _, false )

  val sinPiBF16TableI = SinPiSim.sinPiTableGeneration( 0, 7, 7, 7 )
  val sinPiBF16Sim = sinPiSimGeneric(sinPiBF16TableI, _, false )
}
