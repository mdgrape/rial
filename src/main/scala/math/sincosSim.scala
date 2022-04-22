//% @file sincosSim.scala
//
// Simulator for sin/cos function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import chisel3.util.log2Up

import spire.math.SafeLong
import spire.math.Numeric
import spire.math.Real
import spire.implicits._

import rial.table._
import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.Rounding._
import rial.arith._

object SinCosSim {

  def sincosSimGeneric(
    isSin: Boolean,
    ts: Seq[FuncTableInt],
    t:  FuncTableInt,
    x:  RealGeneric,
    taylorOrder: Int
  ) : RealGeneric = {

    val spec = x.spec
    val exW  = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias

    if (x.isNaN)      {return RealGeneric.nan(spec)}
    if (x.isInfinite) {return RealGeneric.nan(spec)}

    // ------------------------------------------------------------------------
    // calc x/pi

    val oneOverPiPad = 23
    // 1/2 > 1/pi > 1/4, (1/pi).exNobias == -2, and 2 extra bits for rounding
//     val oneOverPi = math.round(1.0 / math.Pi * (1L << (manW+2+oneOverPiPad))).toBigInt
    val oneOverPi = (Real.one / Real.pi)(manW+2+oneOverPiPad)
    assert(bit(manW+oneOverPiPad, oneOverPi) == 1L)

    // no rounding! We will subtract 0.5 or 1.5 from this x/pi.
    // To keep precision, we should have enough bits here.
    val xOverPiProd          = x.manW1 * oneOverPi
    val xOverPiProdMoreThan2 = bit((1+manW)+(1+manW+oneOverPiPad)-1, xOverPiProd)
    val xOverPiEx            = x.ex - 2 + xOverPiProdMoreThan2

    // here we extend the fraction part one more bit.
    val xOverPi              = xOverPiProd << (1 - xOverPiProdMoreThan2)
    val xOverPiFracW         = (1+manW) + (1+manW+oneOverPiPad) - 1
    //                              remove the top, hidden bit  ^^^

    assert((SafeLong(1) << (xOverPiFracW)) < xOverPi && xOverPi < (SafeLong(1) << (xOverPiFracW+1)))
//     println(f"xOverPi   = ${xOverPi.toDouble * pow(2.0, -xOverPiFracW) * pow(2.0, xOverPiEx-exBias)}")
//     println(f"|x|/Pi    = ${x.toDouble.abs / Pi}")
//     println(f"xOverPiEx = ${xOverPiEx}")

    // ------------------------------------------------------------------------
    // convert full range x/pi into (0, 2)

    // now this is in [0, 2)
    val xOverPiAligned = if(xOverPiEx >= exBias) {
      // remove bits that represents larger than 2
      slice(0, 1+xOverPiFracW, xOverPi << (xOverPiEx - exBias))
    } else {
      xOverPi >> (exBias - xOverPiEx)
    }
//     println(f"xOverPi        = ${xOverPi}, W = ${log2Up(xOverPi)}")
//     println(f"xOverPiFracW   = ${xOverPiFracW}")
//     println(f"xOverPiAligned = ${xOverPiAligned}, W = ${log2Up(xOverPiAligned)}")

    val xOverPiAligned2MSBs = slice(1+xOverPiFracW-2, 2, xOverPiAligned)
//     println(f"xOverPiAligned2MSBs = ${xOverPiAligned2MSBs.toLong.toBinaryString}")
    val xOverPiAlignedMoreThan3over2 = xOverPiAligned2MSBs == 3
    val xOverPiAlignedMoreThan1      = xOverPiAligned2MSBs == 2
    val xOverPiAlignedMoreThan1over2 = xOverPiAligned2MSBs == 1

    // we can already calculate the sign of return value from its position in [0, 2pi)

    val zSgn = if (isSin) {
      val zSgnXAbs = if(xOverPiAlignedMoreThan1 || xOverPiAlignedMoreThan3over2) {1} else {0}
      if(x.sgn == 1) {
        1 - zSgnXAbs
      } else {
        zSgnXAbs
      }
    } else {
      if(xOverPiAlignedMoreThan1 || xOverPiAlignedMoreThan1over2) {1} else {0}
    }

    // ------------------------------------------------------------------------
    // convert full range x/pi into (0, 1/2)

//     println(f"sim:xOverPiAligned = ${xOverPiAligned}")

    // than3/2: 11, than1: 10, than1/2: 01, else: 00
    val yman0 = if(isSin) {
      if (xOverPiAlignedMoreThan3over2 || xOverPiAlignedMoreThan1over2) {
        (SafeLong(1) << (1+xOverPiFracW-1)) - slice(0, 1+xOverPiFracW-1, xOverPiAligned)
      } else {
        slice(0, 1+xOverPiFracW-1, xOverPiAligned)
      }
    } else {
      if (xOverPiAlignedMoreThan3over2 || xOverPiAlignedMoreThan1over2) {
        slice(0, 1+xOverPiFracW-2, xOverPiAligned)
      } else {
        (SafeLong(1) << (1+xOverPiFracW-2)) - slice(0, 1+xOverPiFracW-2, xOverPiAligned)
      }
    }

//     println(f"sim:ymanPos0 = ${slice(0, 1+xOverPiFracW-1, xOverPiAligned)}")
//     println(f"sim:ymanNeg0 = ${(1.toBigInt << (1+xOverPiFracW-1)) - slice(0, 1+xOverPiFracW-1, xOverPiAligned)}")
//     println(f"sim:yman0 = ${yman0}")

    val (yex, yman) = if (yman0 == 0) {
      (0, SafeLong(0))
    } else {
      val yman0W        = binaryWidthSL(yman0)
      val yman0Shift    = 1+xOverPiFracW - yman0W
      val yman0Shifted  = yman0 << yman0Shift
      val yman0RoundBit = xOverPiFracW - manW
      val yman0Rounded  = (yman0Shifted >> yman0RoundBit) + bit(yman0RoundBit-1, yman0Shifted)
      val yman0MoreThan2 = bit(manW+1, yman0Rounded)
      assert((yman0MoreThan2 == 1) || (bit(manW, yman0Rounded) == 1))
//       println(f"sim:yman0W       = ${yman0W}")
//       println(f"sim:yman0Shift   = ${yman0Shift}")
//       println(f"sim:yman0Shifted = ${yman0Shifted}")
//       println(f"sim:yman0Rounded = ${yman0Rounded.toLong.toBinaryString}")

      ((exBias-yman0Shift+yman0MoreThan2).toInt, slice(0, manW, yman0Rounded))
    }

    assert(yex  <= exBias-1)
    assert(yex  != exBias-1 || yman == 0)
    assert(yman < (SafeLong(1)<<manW))

//     println(f"sim:yex  = ${yex}")
//     println(f"sim:yman = ${yman.toLong.toBinaryString}")
//     println(f"y        = ${new RealGeneric(x.spec, 0, yex, yman).toDouble}")
//     println(f"|x|/pi   = ${x.toDouble.abs / Pi}")
// 
//     println(f"0.5 - |x|/pi = ${(0.5 - x.toDouble.abs / Pi)}")
//     println(f"|x|/pi - 0.5 = ${(x.toDouble.abs / Pi - 0.5)}")
//     println(f"|x|/pi - 0.5 = ${(x.toFloat.abs / Pi.toFloat - 0.5f)}")
//     println(f"1.5 - |x|/pi = ${(1.5 - x.toDouble.abs / Pi)}")
//     println(f"|x|/pi - 1.5 = ${(x.toDouble.abs / Pi - 1.5)}")
// 
//     println(f"cos(x)   = ${cos(x.toDouble)}")
//     println(f"cos(yPi) = ${(if(zSgn == 1) {-1} else {1}) * sin(new RealGeneric(x.spec, 0, yex, yman).toFloat * Pi.toFloat)}")


    // -----------------------------------------------------------------------
    // calculate sin(yPi) using a table approx of `sin(yPi)/4y`
    //
    // sin(yPi)/4y is in [0.5, 1) if y is in [0, 1/2)

    // align y [0, 1/2) into [0, 1)

    if(yex == exBias-1) {
      assert(yman == 0)
      return new RealGeneric(spec, zSgn, exBias, 0)
    }
    if(yex == 0) {
      assert(yman == 0)
      return new RealGeneric(spec, 0, 0, 0)
    }

    val ymanW1   = yman + (SafeLong(1) << manW)
    val yAligned = ymanW1 >> (exBias-1 - yex) // we already checked y != 1/2
    assert(bit(manW, yAligned) != 1)

    val adrW      = t.adrW
    val nOrder    = t.nOrder
    val bp        = t.bp
    val extraBits = bp - manW
    val fracW     = manW + extraBits
    val order     = if(adrW >= manW) { 0 } else { nOrder }

    if(manW <= adrW && nOrder != 0) {
      println("WARNING: table address width >= mantissa width, but polynomial" +
              " order is not zero. Polynomial order is set to zero.")
    }
    assert(extraBits >= 1, "The implementation requires extraBits >= 1.")

    // w = sin(Piy)/y. later we multiply this with w.

    val wex = exBias-1 // w is in [0.5, 1).
    val wmanW1 = if (order == 0) {
      val adr   = yAligned.toInt
      val res0  = t.interval(adr).eval(0L, 0)
      val res = if (res0<0) {
          println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x = ${x.toDouble}, sin(x) = ${sin(x.toDouble)}")
          SafeLong(0)
        } else if (res0 >= (SafeLong(1)<<fracW)) {
          println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x = ${x.toDouble}, sin(x) = ${sin(x.toDouble)}")
          maskSL(fracW)
        } else {
          SafeLong(res0)
        }

      val wmanW1 = slice(fracW-manW-1, manW+1, res)
      wmanW1

    } else {
      val dxbp = manW-adrW-1
      val d    = slice(0, manW-adrW, yAligned) - (SafeLong(1)<<dxbp)
      val adr  = slice(manW-adrW, adrW, yAligned).toInt

      val res0 = t.interval(adr).eval(d.toLong, dxbp)
      val res = if (res0 < 0) {
          println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x = ${x.toDouble}, sin(x) = ${sin(x.toDouble)}")
          SafeLong(0)
        } else if (res0 >= (SafeLong(1)<<fracW)) {
          println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x = ${x.toDouble}, sin(x) = ${sin(x.toDouble)}")
          maskSL(fracW)
        } else {
          SafeLong(res0)
        }

      val wmanW1 = slice(fracW-manW-1, manW+1, res)
      wmanW1
    }

    // restore z = w * y = sin(Pi*y)/y * y = sin(Pi*y) = sin(x)

    val zProd = wmanW1 * ymanW1

    val zMoreThan2 = bit(((manW+1)*2-1).toInt, zProd)
    val zRounded   = (zProd >> (manW+zMoreThan2)) +
                     bit((manW+zMoreThan2-1).toInt, zProd)
    val zMoreThan2AfterRound = bit(manW+2-1, zRounded)
    val zExInc     = zMoreThan2 + zMoreThan2AfterRound
    val zManW1     = if(zMoreThan2AfterRound == 1) {SafeLong(1) << manW} else {zRounded}

    val zMan = slice(0, manW, zManW1)
    val zEx  = wex - exBias + yex - exBias + zExInc + exBias + 2 // 2 for table

    new RealGeneric(x.spec, zSgn, zEx, zMan)
  }

  def calcTaylorThreshold(manW: Int, taylorOrder: Int): Int = {
    if (taylorOrder <= 2) {
      // sin(pix) = pix - pi^3x^3 / 6 + ...
      //          = pix (1 - pi^2x^2/6 + ...)
      // pi^2x^2/6 < 2^-manW
      // pi^2x^2/6 < 1.65 x^2 < 2x^2 <= 2^-manW
      // x <= 2^(-manW+1)/2
      // x.ex <= (-manW+1)/2-1
      math.floor((-manW+1)/2 - 1).toInt
    } else if (taylorOrder <= 4) {
      // sin(pix) = pix - pi^3x^3 / 6 + pi^5x^5 / 120 - ...
      //          = pix (1 - pi^2x^2/6 + pi^4x^4 / 120)
      // pi^4x^4 / 120 < 2^-manW
      // pi^4x^4 / 120 < 0.812 x^4 <= x^4 <= 2^-manW
      // 2^xex * 1.xman <= 2^-manW/4
      // x.ex <= -manW/4 - 1
      math.floor(-manW/4 - 1).toInt
    } else {
      assert(taylorOrder == 5, "taylorOrder should be <= 5")
      // sin(pix) = pix - pi^3x^3 / 6 + pi^5x^5 / 120 - pi^7x^7/7!
      //          = pix (1 - pi^2x^2 / 6 + pi^4x^4 / 120 - pi^6x^6/7!)
      // pi^6x^6 / 5040 < 2^-manW
      // pi^6x^6 / 5040 < 0.190x^6 < 2^-2 x^6 < 2^-manW
      // x < 2^(-(manW+2)/6)
      math.floor(-(manW+2) / 6).toInt
    }
  }

  // number of tables depending on the exponent and linearThreshold
  def calcExAdrW(spec: RealSpec, taylorOrder: Int): Int = {
    //      .--- table interp --. .-----taylor------.
    // ex = -2 ~ taylorThreshold, taylorThreshold-1 ~ 0

    val taylorThreshold = calcTaylorThreshold(spec.manW, taylorOrder)
    val nTables = -2 - taylorThreshold + 1
    log2Up(nTables)
  }

  def sincosTableGeneration1(
      order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    val tableD = new FuncTableDouble( x0 => {
      val x = x0 / 2.0 // convert [0, 1) into [0, 1/2)
      val z = if(x == 0) {Pi/4.0} else {sin(Pi*x)/(4.0*x)}
      assert(0.5 <= z && z < 1.0, f"x = ${x}, sin(Pi*x) = ${sin(Pi*x)}, z = ${z}")
      z
    }, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  def sincosTableGeneration(
      order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None,
      taylorOrder: Int
    ) = {
    val taylorThreshold = calcTaylorThreshold(manW, taylorOrder)

    if(adrW >= manW) {assert(order == 0)}

    val maxCalcWidth = (-2 to taylorThreshold by -1).map(exponent => {
        val tableD = new FuncTableDouble( x => scalb(sin(Pi * scalb(1.0 + x, exponent)), -exponent-3), order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
      val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

    (-2 to taylorThreshold by -1).map( i => {
      val tableD = new FuncTableDouble( x => scalb(sin(Pi * scalb(1.0+x, i)), -i-3), order )
      tableD.addRange(0.0, 1.0, 1<<adrW)
      new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
    })
  }
}
