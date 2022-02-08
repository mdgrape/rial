//% @file sinSim.scala
//
// Simulator for sin function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

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

object MathFuncSinSim {

  def sinSimGeneric(
    ts: Seq[FuncTableInt],
    x:  RealGeneric,
    useCubicTerm: Boolean = false
  ) : RealGeneric = {

    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias

    if (x.isNaN)      return RealGeneric.nan(x.spec)
    if (x.isInfinite) return RealGeneric.nan(x.spec)

    // ------------------------------------------------------------------------
    // round everything into ex = -inf to -2

    // 1/2 > 1/pi > 1/4, (1/pi).exNobias == -2
    val oneOverPi = (1.0 / math.Pi * (1L << (manW+2))).toLong

    val xOverPiProd = x.manW1 * oneOverPi
    val xOverPiProdMoreThan2 = bit((manW+1)*2-1, xOverPiProd)
    val xOverPiRounded = (xOverPiProd >> (manW + xOverPiProdMoreThan2)) +
                         bit(manW + xOverPiProdMoreThan2 - 1, xOverPiProd)
    val xOverPiProdMoreThan2AfterRound = bit(manW+1, xOverPiRounded)

    assert((xOverPiProdMoreThan2AfterRound == 1) || (bit(manW, xOverPiRounded) == 1))

    val xOverPiMan      = slice(0, manW, xOverPiRounded)
    val xOverPiExNobias = x.ex - x.spec.exBias - 2 +
      xOverPiProdMoreThan2 + xOverPiProdMoreThan2AfterRound

    // ------------------------------------------------------------------------
    // check its value

    // skip large x (2 <= x)
    if (1 <= xOverPiExNobias) {
      if (xOverPiExNobias == 1 && xOverPiMan == 0) { // x = 2.0
        return RealGeneric.zero(x.spec)
      } else {
        return RealGeneric.nan(x.spec)
      }
    }

    if(0 <= xOverPiExNobias && xOverPiMan == 0) {
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

    val zSgn = if ((x.sgn == 1) || (xOverPiExNobias == 0)) { 1 } else { 0 }

    // ------------------------------------------------------------------------
    // round everything into ex = -inf to -2

    val (xex, xman) = if (xOverPiExNobias == 0) { // 1 ~ 2

      assert(xOverPiMan != 0)

      val halfrange = if (bit(manW-1, xOverPiMan) == 1) {
        (1 << manW) - xOverPiMan
      } else {
        xOverPiMan
      }
      assert(0 <= halfrange && halfrange <= (1<<manW))

      // normalize
      val halfrangeMSB = halfrange.toLong.toBinaryString.length()
      val shiftW = manW + 1 - halfrangeMSB
      // if xOverPiMan == 1, we need 1<<manW for normalization
      assert((1 <= shiftW && shiftW <= manW))

      val newex  = if(halfrange == 0) { -exBias } else { -shiftW }
      val newmanW1 = (halfrange << shiftW)
      val newman   = newmanW1 - (1 << manW)

      assert((1<<manW) <= newmanW1)

      (newex.toLong, newman.toLong)

    } else if (xOverPiExNobias == -1) { // 0.5 ~ 1

      val halfrange = (1 << manW) - xOverPiMan
      assert(0 <= halfrange && halfrange <= (1<<manW))

      val halfrangeMSB = halfrange.toLong.toBinaryString.length()
      val shiftW = manW+1 - halfrangeMSB

      // here, man can be zero. if man==0, halfrange == 1<<manW and
      // no shift will be required.
      assert(0 <= shiftW && shiftW <= manW)

      val newex  = -1-shiftW
      val newman = (halfrange << shiftW) - (1 << manW)

      (newex.toLong, newman.toLong)

    } else { // already small enough
      (xOverPiExNobias.toLong, xOverPiMan.toLong)
    }

    assert(xex  <= -1)
    assert(xman <  (1<<manW))
    if (xex == -1) {assert(xman == 0)}

    // ------------------------------------------------------------------------
    // now, {xex, xman} is in [0, 1/2].
    // Do polynomial or taylor approximation.

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

      // 2-bit error here!

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
    // sin(pi*x) = pi*x - pi^3 * x^3 / 6
    // so the condition is:
    //               pi*x * 2^-manW     > pi^3 * x^3 / 6
    //                      2^-manW     > pi^2 * x^2 / 6
    //           6 / pi^2 * 2^-manW     > x^2
    //   log2(6) - 2log2(pi) - manW     > 2*log2(x)
    //   log2(6) - 2log2(pi) - manW     > 2*(log2(2^x.ex) + log2(1.0 + x.man))
    //   log2(6) - 2log2(pi) - manW / 2 > x.ex + log2(1.0 + x.man)
    //
    // Here, the maximum value of log2(1.0+x.man) is log2(2.0 - 2^-manW).
    //   Thus the result is
    //
    //   log2(6) - 2log2(pi) - manW / 2 - log2(2.0 - 2^-manW) > x.ex
    //
    // The left hand side value is the linear-threshold. In case of FP32, the
    // value is approximately -12.859.
    //
    math.floor(
      (log2D(6.0) - 2 * log2D(Pi) - manW) / 2 - log2D(2.0 - math.pow(2.0, -manW))
    ).toInt
  }
  def calcCubicThreshold(manW: Int): Int = {
    // sin(pi*x) = pi*x - pi^3 * x^3 / 6 + pi^5 * x^5 / 120
    // so the condition is:
    //                 pi*x * 2^-manW     > pi^5 * x^5 / 120
    //                        2^-manW     > pi^4 * x^4 / 120
    //           120 / pi^4 * 2^-manW     > x^4
    //   log2(120) - 4log2(pi) - manW     > 4*log2(x)
    //   log2(120) - 4log2(pi) - manW     > 4*(log2(2^x.ex) + log2(1.0 + x.man))
    //   log2(120) - 4log2(pi) - manW / 4 > x.ex + log2(1.0 + x.man)
    //
    // Here, the maximum value of log2(1.0+x.man) is log2(2.0 - 2^-manW).
    //   Thus the result is
    //
    //   log2(120) - 4log2(pi) - manW / 4 - log2(2.0 - 2^-manW) > x.ex
    //
    // The left hand side value is the linear-threshold. In case of FP32, the
    // value is approximately -6.674.

    math.floor(
      (log2D(120.0) - 4 * log2D(Pi) - manW) / 4 - log2D(2.0 - math.pow(2.0, -manW))
    ).toInt
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

  def sinTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int,
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
}
