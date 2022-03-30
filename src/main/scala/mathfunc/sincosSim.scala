//% @file sincosSim.scala
//
// Simulator for sin/cos function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

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

object MathFuncSinCosSim {

  def sincosSimGeneric(
    isSin: Boolean,
    ts: Seq[FuncTableInt],
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

    // ------------------------------------------------------------------------
    // calculate sin(y*Pi).

    val taylorThreshold = calcTaylorThreshold(manW, taylorOrder)

    if (yex == 0) { // sin(0) = 0
//       println(f"x = ${x.toDouble}, y = ${x.toDouble / Pi} equiv 0")

      return RealGeneric.zero(x.spec)

    } else if (yex == exBias-1 && yman == 0) { // sin(pi/2) = 1
//       println(f"x = ${x.toDouble}, y = ${x.toDouble / Pi} equiv 1/2")

      return new RealGeneric(x.spec, zSgn, exBias, 0)

    } else if (yex - exBias < taylorThreshold) {
//       println("using Taylor")
      //
      // sin(pi*y) = pi*y - pi^3 * y^3 / 6 + pi^5y^5/120
      //           = pi*y * (1 - pi^2 * x^2 / 6 + pi^4y^4/120)
      //
      val coefPad    = 2 // for precision after rounding
      val fracW      = manW+coefPad
      val coef1Ex    = 1 + exBias
      val coef1ManW1 = Real.pi(fracW-(coef1Ex-exBias))
//       val coef1ManW1 = math.round(Pi * (1<<(fracW-(coef1Ex-exBias)))).toLong
      val coef3Ex    = 0 + exBias
//       val coef3ManW1 = math.round(Pi * Pi / 6.0 * (1<<(fracW-(coef3Ex-exBias)))).toLong
      val coef3ManW1 = (Real.pi * Real.pi / Real(6))(fracW-(coef3Ex-exBias))
      val coef5Ex    = -1 + exBias
//       val coef5ManW1 = math.round(pow(Pi, 4) / 120.0 * (1<<(fracW-(coef5Ex-exBias)))).toLong
      val coef5ManW1 = (Real.pi.pow(4) / Real(120))(fracW-(coef5Ex-exBias))
      assert(bit(fracW, coef1ManW1) == 1)
      assert(bit(fracW, coef3ManW1) == 1)
      assert(bit(fracW, coef5ManW1) == 1)

      // returning W = 1+fracW
      val multiply = (xFracW: Int, xmanW1: SafeLong, yFracW: Int, ymanW1: SafeLong) => {
        assert(bit(xFracW, xmanW1) == 1)
        assert(bit(yFracW, ymanW1) == 1)
        val zProd          = xmanW1 * ymanW1
        val zProdMoreThan2 = bit((1+xFracW)+(1+yFracW)-1, zProd)
        val zProdShift     = (xFracW + yFracW - fracW + zProdMoreThan2).toInt
        val zProdRounded   = (zProd >> zProdShift) + bit(zProdShift-1, zProd)
        val zProdMoreThan2AfterRound = bit(2+fracW-1, zProdRounded)
        val zExInc = zProdMoreThan2 + zProdMoreThan2AfterRound
        val zManW1 = if(zProdMoreThan2AfterRound == 1) {SafeLong(1) << fracW} else {zProdRounded}
        assert(bit(fracW, zManW1) == 1)
        assert(zExInc == 1 || zExInc == 0)
        (zExInc.toInt, zManW1)
      }

      // 8 ops in 5 steps:
      //
      // y-+-> y^2 -+-> y^4       -> pi^4y^4/120 -+-> 1-pi^2y^2/6+pi^4y^4/120 -+->
      //   |        +-> pi^2y^2/6 -> 1-pi^2y^2/6 -+                            |
      //   +-> piy  -----------------------------------------------------------+
      //

      // y^2
      val (ySqExInc, ySqManW1) = multiply(manW, (SafeLong(1)<<manW) + yman,
                                          manW, (SafeLong(1)<<manW) + yman)
      val ySqEx = yex + yex - exBias + ySqExInc
      assert(bit(fracW, ySqManW1) == 1)
//       println(f"sim:ySqEx    = ${ySqEx}")
//       println(f"sim:ySqManW1 = ${ySqManW1.toLong.toBinaryString}%25s")

//       println("piy")
      // pi*y
      val (piyExInc, piyManW1) = multiply(manW, (SafeLong(1)<<manW) + yman,
                                          fracW, coef1ManW1)
      val piyEx = yex + coef1Ex - exBias + piyExInc
      assert(bit(fracW, piyManW1) == 1)
//       println(f"sim:piyEx    = ${piyEx}")
//       println(f"sim:piyManW1 = ${piyManW1.toLong.toBinaryString}%25s")

      if(taylorOrder <= 2) {
        // sin(x) = x = piy
        val taylorEx  = piyEx
        val taylorMan = slice(0, manW, piyManW1)

        return new RealGeneric(x.spec, zSgn, taylorEx.toInt, taylorMan)
      }

      // ---------------- 3rd order --------------------

      // pi^2y^2/6
      val (c3ExInc, c3ManW1) = multiply(fracW, ySqManW1, fracW, coef3ManW1) // XXX
      val c3Ex = ySqEx + coef3Ex - exBias + c3ExInc
      assert(bit(fracW, c3ManW1) == 1)
//       println(f"sim:c3Ex    = ${c3Ex}")
//       println(f"sim:c3ManW1 = ${c3ManW1.toLong.toBinaryString}%25s")

      // 1 - pi^2y^2/6
      assert(c3Ex < exBias)
      val c3Shift = exBias - c3Ex
      val c3Aligned = if(c3Shift > 63) { SafeLong(0) } else {
        (c3ManW1 >> c3Shift) + bit(c3Shift-1, c3ManW1)
      }
      val oneMinusC3 = (SafeLong(1)<<fracW) - c3Aligned
//       println(f"sim:1-c3 = ${oneMinusC3.toLong.toBinaryString}")

      if(taylorOrder <= 4) {
//         println("taylorOrder <= 4. return 1 - c3")
        // piy * (1 - pi^2y^2/6)
        //
        val oneMinusC3MoreThan1 = bit(fracW, oneMinusC3) // == (c3Aligned == 0)
        val oneMinusC3ManW1 = oneMinusC3 << (1 - oneMinusC3MoreThan1)
        val oneMinusC3Ex = exBias - 1 + oneMinusC3MoreThan1

        val (taylorExInc, taylorManW1) = multiply(fracW, piyManW1, fracW, oneMinusC3ManW1)
        val taylorManW1Rounded = (taylorManW1 >> coefPad) + bit(coefPad-1, taylorManW1)
        val taylorManW1MoreThan2AfterRound = bit(2+manW-1, taylorManW1Rounded)

        val taylorEx  = piyEx + oneMinusC3Ex - exBias + taylorExInc + taylorManW1MoreThan2AfterRound
        val taylorMan = slice(0, manW, taylorManW1Rounded)

        return new RealGeneric(x.spec, zSgn, taylorEx.toInt, taylorMan)
      }

      assert(taylorOrder == 5, "taylorOrder should be <= 5")

      // ---------------- 5th order --------------------

      // y^4
      val (yQdExInc, yQdManW1) = multiply(fracW, ySqManW1, fracW, ySqManW1)
      val yQdEx = ySqEx + ySqEx - exBias + yQdExInc
      assert(bit(fracW, yQdManW1) == 1)
//       println(f"sim:yQdEx    = ${yQdEx}")
//       println(f"sim:yQdManW1 = ${yQdManW1.toLong.toBinaryString}%25s")

      // pi^4y^4/120
      val (c5ExInc, c5ManW1) = multiply(fracW, yQdManW1, fracW, coef5ManW1)
      val c5Ex = yQdEx + coef5Ex - exBias + c5ExInc
      assert(bit(fracW, c5ManW1) == 1)
//       println(f"sim:c5Ex    = ${c5Ex}")
//       println(f"sim:c5ManW1 = ${c5ManW1.toLong.toBinaryString}%25s")

      // 1 - pi^2y^2/6 + pi^4y^4/120
      // ~ 1 - 1.645y^2 + 0.8117y^4
      assert(c5Ex < exBias)
      val c5Shift = exBias - c5Ex
      val c5Aligned = if(c5Shift > 63) { SafeLong(0) } else {
        (c5ManW1 >> c5Shift) + bit(c5Shift-1, c5ManW1)
      }
      val oneMinusC3PlusC5 = oneMinusC3 + c5Aligned
      assert(c3Aligned >= c5Aligned)
      assert(oneMinusC3PlusC5 <= (1<<fracW)) // 1 - pi^2y^2/6 + pi^4y^4/120 <= 1

//       println(f"sim:1-c3+c5 = ${oneMinusC3PlusC5.toLong.toBinaryString}")

      val oneMinusC3PlusC5MoreThan1 = bit(fracW, oneMinusC3PlusC5)
      val oneMinusC3PlusC5ManW1 = oneMinusC3PlusC5 << (1 - oneMinusC3PlusC5MoreThan1)
      val oneMinusC3PlusC5Ex = exBias - 1 + oneMinusC3PlusC5MoreThan1

//       println(f"sim:1-c3+c5 Ex    = ${oneMinusC3PlusC5Ex}")
//       println(f"sim:1-c3+c5 ManW1 = ${oneMinusC3PlusC5ManW1.toLong.toBinaryString}")

      // piy * (1 - pi^2y^2/6 + pi^4y^4/120)
      val (taylorExInc, taylorManW1) = multiply(fracW, piyManW1, fracW, oneMinusC3PlusC5ManW1)
      val taylorManW1Rounded = (taylorManW1 >> coefPad) + bit(coefPad-1, taylorManW1)
      val taylorManW1MoreThan2AfterRound = bit(2+manW-1, taylorManW1Rounded)

      val taylorEx  = piyEx + oneMinusC3PlusC5Ex - exBias + taylorExInc + taylorManW1MoreThan2AfterRound
      val taylorMan = slice(0, manW, taylorManW1Rounded)

//       println(f"sim:taylorEx    = ${taylorEx}")
//       println(f"sim:taylorManW1 = ${taylorMan.toLong.toBinaryString}")

      return new RealGeneric(x.spec, zSgn, taylorEx.toInt, taylorMan)

    } else { // table interpolation
      assert(taylorThreshold <= yex - exBias)
//       println("using table")

      // table interpolate for x in [0, 1/2) (x.ex = -2, -3, -4, -5 if FP32)
      val exadr = (exBias - yex - 2).toInt
      val t = ts(exadr)

      val adrW      = t.adrW
      val nOrder    = t.nOrder
      val bp        = t.bp
      val extraBits = bp - manW
      val fracW     = manW + extraBits
      val order     = if(adrW >= manW) { 0 } else { nOrder }

      if(manW <= adrW && nOrder != 0) {
        println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
      }

      val (zEx, zman) = if (order == 0) {
        val adr   = yman.toInt
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

        val lessThanHalf = if(bit(fracW-1, res) == 0) { 1 } else { 0 }
        val ex    = yex+2-lessThanHalf
        val man   = (res << (1+lessThanHalf)).toLong - (SafeLong(1) << fracW)

        (ex.toInt, man)

      } else {
        val dxbp = manW-adrW-1
        val d    = slice(0, manW-adrW, yman) - (SafeLong(1)<<dxbp)
        val adr  = slice(manW-adrW, adrW, yman).toInt

//         println(f"sim:dx  = ${d.toLong.toBinaryString}")
//         println(f"sim:adr = ${exadr.toBinaryString}|${adr.toBinaryString}")

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
//         println(f"sim:zres = ${res.toLong.toBinaryString}")
        val lessThanHalf = if(bit(fracW-1, res) == 0) { 1 } else { 0 }
        ((yex+2-lessThanHalf).toInt, (res << (1+lessThanHalf)) - (SafeLong(1)<<fracW))
      }

      val zmanRound = if (extraBits>0) {(zman>>extraBits) + bit(extraBits-1, zman)} else {zman}
      val zMan = slice(0, manW, zmanRound)
      val zManMoreThan2 = bit(manW, zmanRound).toInt
//       println(f"sim:zmanRound = ${zmanRound.toLong.toBinaryString}")

      new RealGeneric(x.spec, zSgn, zEx + zManMoreThan2, zMan)
    }
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
