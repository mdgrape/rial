//% @file acosSim.scala
//
// Simulators for acos(x) function
// Copyright (C) Toru Niina RIKEN BDR 2022
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

object MathFuncACosSim {

  // assuming x is in range [-1, 1]. otherwise, NaN.
  //
  // acos(x) = pi/2 - x - x^3/6 - 3x^5/40
  //         < pi/2 - x - x^3/4 -  x^5/8
  //         = pi/2 - x - 2^-2 x^3 - 2^-3 x^5
  //
  // pi/2 = 2^0 * 1.57..
  //
  // if x < 2^-23 ... acos(x) = pi/2
  // if x < 2^-8  ... acos(x) = pi/2 - x
  // if x < 2^-4  ... acos(x) = pi/2 - x - x^3/6
  // otherwise    ... use tables
  //
  // Note: acos(-x) = pi - acos(x)
  //
  def acosSimGeneric( ts : Seq[FuncTableInt],
    tEdge1: FuncTableInt, tEdge2: FuncTableInt, tEdge3: FuncTableInt,
    x: RealGeneric ) : RealGeneric = {
    val adrW   = ts(0).adrW
    val nOrder = ts(0).nOrder
    val bp     = ts(0).bp

    val expW   = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val extraBits = if (nOrder==0) {0} else {(bp - manW)}

    val sgn    = x.sgn
    val exRaw  = x.ex
    val ex     = exRaw-exBias
    val man    = x.man

//     println(f"x = ${x.sgn}|${x.ex}(${ex})|${man.toLong.toBinaryString}")
//     println(f"x = ${x.toDouble}, acos(x) = ${acos(x.toDouble)}")

    if(x.isNaN) {
      return RealGeneric.nan(x.spec)
    }
    if(x.isInfinite) {
      return RealGeneric.nan(x.spec)
    }
    if(0 <= ex) { // 1 <= |x|
      if (sgn == 0) {
        return new RealGeneric(x.spec, 0.0)
      } else {
        return new RealGeneric(x.spec, Pi)
      }
    }

    val zSgn = 0 // returns [0, 2pi)

    // here ex < 0.

    val constThreshold  = -manW                   // -23, if FP32
    val linearThreshold = calcLinearThreshold(manW) // 8
    val halfPi = new RealGeneric(x.spec, Pi * 0.5)

//     println(f"constant threshold = ${constThreshold}")
//     println(f"linear   threshold = ${linearThreshold}")

    // here the x.sgn has no effect because |x| << 1.
    if (ex < constThreshold) {
//       println(f"x < constantThreshold(${pow(2.0, -constThreshold)})")
      return halfPi
    } else if (ex < linearThreshold) { // Pi/2 - x
//       println(f"x < linearThreshold(${pow(2.0, linearThreshold)})")
      // Pi / 2 ~ 1.57... and x < 2^(linearThreshold). we don't need to check MoreThan2.

      val shiftOffset = 3
      val calcW       = (manW + 1) + shiftOffset
      val halfPiMan   = ((1 << manW) + halfPi.man) << shiftOffset
      val xManAligned = ((1 << manW) + x.man) >> (- ex - shiftOffset)
      val xManSigned  = if (sgn == 0) { ~(xManAligned.toLong) + 1 } else { xManAligned.toLong }
      val zMan0       = slice(0, calcW, halfPiMan + xManSigned) // rm the leading 1

      if(bit(calcW, zMan0) == 1) {
        val zEx  = exBias + 1
        val zMan = slice(0, manW, roundBySpec(RoundSpec.roundToEven, shiftOffset+1, zMan0))
        return new RealGeneric(x.spec, 0, zEx, zMan)
      }
      if(bit(calcW-1, zMan0) == 1) {
        val zEx  = exBias
        val zMan = slice(0, manW, roundBySpec(RoundSpec.roundToEven, shiftOffset, zMan0))
        return new RealGeneric(x.spec, 0, zEx, zMan)
      } else {
        val zEx  = exBias - 1
        val zMan = slice(0, manW, roundBySpec(RoundSpec.roundToEven, shiftOffset, zMan0))
        return new RealGeneric(x.spec, 0, zEx, zMan)
      }
    } else if(ex == -1 && x.sgn == 0) {
      // x is in [0.5, 1).

      val calcW = tEdge1.bp
      assert(tEdge1.bp == tEdge2.bp)
      assert(tEdge1.bp == tEdge3.bp)

      val xNegMan   = (1L<<manW) - x.man
      val xNegManW  = xNegMan.toLong.toBinaryString.length
      val xNegManEx = manW + 1 - xNegManW
      val baseEx    = 1 + floor((xNegManEx+1)/2).toInt

      val res0 = if (xNegManW < 8) {
        assert(adrW == 8)
        val dxbp = 8 - adrW - 1
        val d    = 0
        val adr  = slice(dxbp+1, adrW, man).toInt
        tEdge1.interval(adr).eval(d.toLong, dxbp)
      } else if (xNegManW < 16) {
        val dxbp = 16 - adrW - 1
        val d    = slice(0, dxbp+1, man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, man).toInt
        tEdge2.interval(adr).eval(d.toLong, dxbp)
      } else {
        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, man).toInt
        tEdge3.interval(adr).eval(d.toLong, dxbp)
      }

      val res0MoreThanHalf = bit(calcW-1, res0)
      val res0MoreThanQuat = bit(calcW-2, res0)

      val (zExInc, zManRounded) = if(res0MoreThanHalf == 1) {
        ( 0, (res0 >> extraBits) + bit(extraBits-1, res0))
      } else if (res0MoreThanQuat == 1) {
        val res = res0 << 1
        (-1, (res >> extraBits) + bit(extraBits-1, res))
      } else {
        val res = res0 << 2
        (-2, (res >> extraBits) + bit(extraBits-1, res))
      }

      val zManRoundedMoreThan2 = bit(manW+1, zManRounded)
      val zEx  = baseEx + zExInc + zManRoundedMoreThan2
      val zMan = slice(0, manW, zManRounded)

      assert(bit(manW, zManRounded) == 1 || zManRoundedMoreThan2 == 1)

      return new RealGeneric(x.spec, 0, zEx.toInt + exBias, zMan)

    } else { // use table.
//       println(f"linearThreshold < x. use table.")

      // ex = -1 -> idx = 0,
      // ex = -2 -> idx = 1,
      val exAdr = -ex-1
      val t = ts(exAdr)

      val adrW      = t.adrW
      val nOrder    = t.nOrder
      val bp        = t.bp
      val extraBits = bp - manW
      val calcW     = manW + extraBits

      val order =
        if (adrW>=manW) {
          if (nOrder != 0)
            println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
          0
        } else {
          nOrder
        }

      val (zEx, zman) = if (order == 0) {
        val adr   = man.toInt

        if (sgn == 1) {
          val piFixed = (math.round(Pi * (1<<calcW))).toLong
          val res0    = piFixed - t.interval(adr).eval(0L, 0)
          // res0 range is [1.57, 3.14]. exponent should be 0 or 1

          val shift = bit(calcW+1, res0)
          val resShifted = res0 >> shift
          val res = resShifted - (1<<calcW)

          (shift.toInt, res)
        } else {
          val res0  = t.interval(adr).eval(0L, 0)
          val shift      = calcW+1 - res0.toLong.toBinaryString.length
          val resShifted = if(shift > 0) {res0 << shift} else {res0}
          val res = resShifted - (1<<calcW)

          (-shift.toInt, res)
        }

      } else { // non-zero order

        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, man).toInt

        val halfPiFixed = math.round(Pi * 0.5 * (1<<calcW))

        // pi/2 - acos(x)
        val res0 = t.interval(adr).eval(d.toLong, dxbp) << 1
        val res  = if (sgn == 1) {
          halfPiFixed + res0
        } else {
          halfPiFixed - res0
        }
        val shift = calcW+2 - res.toLong.toBinaryString.length
        val resShifted = ((res << shift).toLong) >> 1

        ((-shift+1).toInt, resShifted - (1 << calcW))
      }
//       println(f"zex    = ${zEx}")
//       println(f"zman   = ${zman.toLong.toBinaryString}")

      val zmanRound = if (extraBits>0) {(zman >> extraBits) + bit(extraBits-1, zman)} else {zman}

//       println(f"zmanR  = ${zmanRound.toLong.toBinaryString}")

      val z = if (zman<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
        0L
      } else if (zmanRound >= (1L<<manW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
        maskL(manW)
      } else {
        zmanRound
      }
//       println(f"Sim: zman = ${z.toBinaryString}")

      return new RealGeneric(x.spec, zSgn, zEx + exBias, SafeLong(z))
    }
  }

  def calcLinearThreshold(manW: Int): Int = {
    math.ceil((-manW - 1) / 3.0).toInt // (-23 - 1) / 3.0 = 8
  }

  // number of tables depending on the exponent and linearThreshold
  def calcExAdrW(spec: RealSpec): Int = {
    val linearThreshold = calcLinearThreshold(spec.manW)
    // from -1 to linearThreshold (-8 in FP32), 126 to 119 if biased.
    // 0b01111110 to 0b01110111, 4 bits required
    log2Up(abs(linearThreshold)+1)
  }

  def acosTableEdge1( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
      // xman in [1.111'1111'1111'1111'0000'0000, 1.111'1111'1111'1111'1111'1111]
      val f = (x: Double) => {
        val f64      = RealSpec.Float64Spec
        val nx       = x * pow(2.0, 8 - manW)
        val xi       = java.lang.Double.doubleToRawLongBits(nx)
        val xex      = slice(f64.manW, f64.exW, xi) - f64.exBias.toLong
        val baseline = pow(2.0, -(1+floor((xex+1)/2).toInt))
        val z        = acos(1.0 - nx)
        z * baseline
      }
      val tableD = new FuncTableDouble( f, order )
      tableD.addRange(0.0, 1.0, 1<<manW)
      new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }
  def acosTableEdge2( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
      // xman in [1.111'1111'1111'1111'0000'0000, 1.111'1111'1111'1111'1111'1111]
      val f = (x: Double) => {
        val f64      = RealSpec.Float64Spec
        val nx       = x * pow(2.0, 16 - manW)
        val xi       = java.lang.Double.doubleToRawLongBits(x)
        val xex      = slice(f64.manW, f64.exW, xi) - f64.exBias.toLong
        val baseline = pow(2.0, -(1+floor((xex+1)/2).toInt))
        val z        = acos(1.0 - x)
        z * baseline
      }
      val tableD = new FuncTableDouble( f, order )
      tableD.addRange(0.0, 1.0, 1<<manW)
      new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }
  def acosTableEdge3( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
      // xman in [1.111'1111'1111'1111'0000'0000, 1.111'1111'1111'1111'1111'1111]
      val f = (x: Double) => {
        val f64      = RealSpec.Float64Spec
        val nx       = x
        val xi       = java.lang.Double.doubleToRawLongBits(x)
        val xex      = slice(f64.manW, f64.exW, xi) - f64.exBias.toLong
        val baseline = pow(2.0, -(1+floor((xex+1)/2).toInt))
        val z        = acos(1.0 - x)
        z * baseline
      }
      val tableD = new FuncTableDouble( f, order )
      tableD.addRange(0.0, 1.0, 1<<manW)
      new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }

  def acosTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    val linearThreshold = calcLinearThreshold(manW)

    if (order == 0 || adrW >= manW) {
      val maxCalcWidth = (-1 to linearThreshold by -1).map(i => {
        val tableD = new FuncTableDouble( x => acos(scalb(1.0 + x, i)), 0)
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

      (-1 to linearThreshold by -1).map( i => {
        val tableD = new FuncTableDouble( x => acos(scalb(1.0 + x, i)), 0)
        tableD.addRange(0.0, 1.0, 1<<manW)
        new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
      })
    } else {
      val maxCalcWidth = (-1 to linearThreshold by -1).map(i => {
        val tableD = new FuncTableDouble( x => ((Pi * 0.5) - acos(scalb(1.0 + x, i))) * 0.5, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

      // ex == -1 corresponds to the range [0.5, 1).
      (-1 to linearThreshold by -1).map( i => {
        val tableD = new FuncTableDouble( x => ((Pi * 0.5) - acos(scalb(1.0 + x, i))) * 0.5, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
      })
    }
  }
}
