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

import rial.math.SqrtSim

object MathFuncACosSim {

  // assuming x, y, and z have the same spec.
  // since acos uses 2 different high-order series expansion,
  // this kind of helper function is useful.
  def multiply(spec: RealSpec,
    xexNobias: Int, xmanW1: Long, yexNobias: Int, ymanW1: Long
    ): (Int, Long) = {
    val manW = spec.manW

    val zProd      = xmanW1 * ymanW1
    val zMoreThan2 = bit(((manW+1)*2-1).toInt, zProd)
    val zRounded   = (zProd >> (manW+zMoreThan2)) +
                     bit((manW+zMoreThan2-1).toInt, zProd)
    val zMoreThan2AfterRound = bit(manW+2-1, zRounded)
    val zExNobias  = xexNobias + yexNobias + zMoreThan2 + zMoreThan2AfterRound
    val zManW1     = if(zMoreThan2AfterRound == 1) {1<<manW} else {zRounded}

    assert(bit(manW, zManW1) == 1)

    (zExNobias.toInt, zManW1)
  }

  // assuming x is in range [-1, 1].
  // otherwise, 0 or pi for positive and negative x, respectively. (to avoid
  // numerical error having a value something like: 1.00000001)
  //
  //   acos( x) = pi/2 - pi/2 + acos(x) = pi/2 - [pi/2 - acos(x)]
  //   acos(-x) = pi - acos(x)          = pi/2 + [pi/2 - acos(x)]
  //
  // for small x, pi/2 - acos(x) = x + x^3/6 + 3x^5/40 + O(x^7)
  // for x < 0.5, pi/2 - acos(x) is approximated by polynomial.
  // for x > 0.5, acos(x) is approximated by polynomial.
  // for x close to 1, use puiseux series:
  //   acos(1-x) = sqrt(2x) * (1 + x/12 + 3x^2/160 + 5x^3/896 + 35x^4/18432 + O(x^5))
  //
  // In case of Taylor series, the condition where pi/2 - acos(x) has enough
  // precision is:
  //   3x^5/40 < 2^-23
  //       x^5 < 2^-20
  //       x   < 2^-4.
  // So if x < 2^-4, that means that x.ex < 2^-5, use Taylor series.
  //
  // In case of Puiseux series, the condition is:
  //   35x^4/18432 < 2^-23
  //     x^4/526.628... < 2^-23
  //     x^4 < 2^-14
  //     x   < 2^-4
  // So if 1-2^-4 < x, that means that x.man(22,20).andR === 1, use Puiseux series.
  //
  def acosSimGeneric(
      ts: Seq[FuncTableInt], tSqrt: FuncTableInt, x: RealGeneric
    ) : RealGeneric = {

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

//     println(f"x = ${x.sgn}|${x.ex}(${ex})|${man.toLong.toBinaryString}")
//     println(f"x = ${x.toDouble}, acos(x) = ${acos(x.toDouble)}")

    // ------------------------------------------------------------------------
    // special value handling

    if(x.isNaN) {
      return RealGeneric.nan(spec)
    }
    if(x.isInfinite) {
      return RealGeneric.nan(spec)
    }
    if(exBias <= x.ex) { // 1 <= |x|
      if (x.sgn == 0) {
        return new RealGeneric(spec, 0.0)
      } else {
        return new RealGeneric(spec, Pi)
      }
    }

    val zSgn = 0 // returns z in [0, 2pi), so always z > 0.

    // ------------------------------------------------------------------------
    // now, ex < 0.

    val constThreshold  = -manW                     // -23, if FP32
    val taylorThreshold = calcTaylorThreshold(manW) // 4, if FP32

//     println(f"constant threshold = ${constThreshold}")
//     println(f"Taylor   threshold = ${taylorThreshold}")

    if ((x.ex - exBias) < constThreshold) {

      // here the x has no effect because |x| << 1.
      return new RealGeneric(spec, Pi * 0.5)

    } else if ((x.ex - exBias) < taylorThreshold) {

      // for small x, pi/2 - acos(x) = x + x^3/6 + O(x^5)
      //                             = x(1 + x^2/6)

      val xexNobias = x.ex - exBias
      val xmanW1    = x.manW1.toLong
      val (xSqExNobias, xSqManW1) =
        MathFuncACosSim.multiply(spec, xexNobias, xmanW1, xexNobias, xmanW1)

      val c1over6ExNobias  = -3
      val c1over6ManW1     = math.round(1.0/6.0 * (1<<(manW+(-c1over6ExNobias))))

      val (xSq6thExNobias, xSq6thManW1) =
        MathFuncACosSim.multiply(spec, xSqExNobias, xSqManW1, c1over6ExNobias, c1over6ManW1)

      assert(xSq6thExNobias < 0)

      val xSq6thAligned  = xSq6thManW1 >> (-xSq6thExNobias)
      val xSq6thPlusOneExNobias = 0
      val xSq6thPlusOneManW1    = xSq6thAligned + (1<<manW)

      val (taylorExNobias, taylorManW1) =
        MathFuncACosSim.multiply(spec, xexNobias, xmanW1, xSq6thPlusOneExNobias, xSq6thPlusOneManW1)
      assert(taylorExNobias < 0)

      val halfPiExNobias = 0
      val halfPiManW1    = math.round(Pi * 0.5 * (1<<manW)).toLong
      assert(bit(manW, halfPiManW1) == 1)

      val taylorAligned  = taylorManW1 >> (-taylorExNobias)
      assert(taylorAligned < halfPiManW1)

      // Let's say we have FP32. For small x < 2^-4,
      //   pi/2 - acos(x) ~ x + x^3/6 < 2^-4 * (1 + 2^-8 / 6) < 2^-3
      //   1 + 2^-3 = 1.125 < pi/2 ~ 1.5707963 < 2 - 2^-3 = 1.875.
      // we don't need any check. exponent never change.
      val taylorSum = if(x.sgn == 1) {
        halfPiManW1 + taylorAligned
      } else {
        halfPiManW1 - taylorAligned
      }

      val taylorEx  = halfPiExNobias + exBias
      val taylorMan = slice(0, manW, taylorSum)

      return new RealGeneric(spec, zSgn, taylorEx, taylorMan)

    } else if(x.ex == exBias-1 && slice(manW-3, 3, x.man) == 7) {
//       println(f"use Puiseux series: |x| = ${x.toDouble.abs}, y = ${1.0 - x.toDouble.abs}")

      // for x close to 1, use puiseux series:
      //   let y = 1 - x
      //   acos(1-y) = sqrt(2y) * (1 + y/12 + 3y^2/160 + 5y^3/896 + O(y^4))
      //             = sqrt(2y) * ((1 + 1/3 * 2^-2 * y) + 3y^2/160 * (1 + 25/21 * 2^-2 * y))

      val oneMinusX    = (1<<(manW+1)) - x.manW1
      // to use sqrt table, we need to normalize 1-x.
      val oneMinusXLen = oneMinusX.toLong.toBinaryString.length
      val ymanW1    = (oneMinusX << (manW+1 - oneMinusXLen)).toLong
      val yman      = ymanW1 - (1<<manW)
      val yexNobias = -(manW+1 - oneMinusXLen)-1 // -1 because oneMinusX has 1bit wider mantissa
      val yex       = yexNobias + exBias

      assert(bit(manW, ymanW1) == 1)
      assert(yman < (1<<manW))
      assert(0 < yex)

//       println(f"yref = ${1.0 - x.toDouble.abs}")
//       println(f"ysim = ${ymanW1.toDouble / (1<<manW) * pow(2.0, yexNobias)}")

      // ----------------------------------------------------------------------
      // sqrt(2y)
      // sqrt table uses the last bit of exponent.

      val adrW      = tSqrt.adrW
      val extraBits = tSqrt.bp - manW
      val y2man  = slice(0, manW+1, ((yex + 1)<<manW) + yman)

      val dxbp   = (manW+1)-adrW-1
      val d      = slice(0, (manW+1)-adrW, y2man) - (SafeLong(1)<<dxbp)
      val adr    = slice((manW+1)-adrW, adrW, y2man).toInt

      // sqrt table returns the mantissa only. we need to add 1<<fracW
      val sqrt2y = tSqrt.interval(adr).eval(d.toLong, dxbp) + (1<<(tSqrt.bp))

      val sqrt2yRounded  = (sqrt2y >> extraBits) + bit(extraBits-1, sqrt2y)
      val sqrt2yMoreThan2AfterRound = bit(manW+2-1, sqrt2yRounded)
      val sqrt2yExNobias = ((yexNobias + 1) >> 1) + sqrt2yMoreThan2AfterRound.toInt
      val sqrt2yManW1    = if(sqrt2yMoreThan2AfterRound == 1) {1<<manW} else {sqrt2yRounded}

//       println(f"sqrt(2y) = ${sqrt(2.0 * (1.0 - x.toDouble.abs))}")
//       println(f"sqrt2y   = ${(sqrt2yManW1.toDouble / (1<<manW)) * pow(2.0, sqrt2yExNobias)}")
      assert(bit(manW, sqrt2yManW1) == 1)

      // ----------------------------------------------------------------------
      // 2^-5 * 3/5 y^2

      val c3over5 = math.round(3.0/5.0 * (1<<(manW+1))) // +1 for normalize
      assert(bit(manW, c3over5) == 1)

      // in the circuit, the width of oneMinusX should be manW+1 because
      // it is the result of (1<<manW+1) - x.manW1. And it is impossible to
      // change the width dynamically, so here we can use ymanW1 without any
      // additional cost.

      val (ySqExNobias, ySqManW1) =
        MathFuncACosSim.multiply(spec, yexNobias, ymanW1, yexNobias, ymanW1)

      val (ySq3over5ExNobias, ySq3over5ManW1) =
        MathFuncACosSim.multiply(spec, ySqExNobias, ySqManW1, -1, c3over5)

      // ----------------------------------------------------------------------
      // 1 + 25/21 * 2^-2 * y
      //
      // here, y < 2^-4. so 1 + 25/21 * 2^-2 * y never exceeds 2.
      // We don't need to check if it exceeds 2 after addition.

      val c25over21 = math.round(25.0/21.0 * (1<<manW))
      assert(bit(manW, c25over21) == 1)

      val (y25over21ExNobias, y25over21ManW1) =
        MathFuncACosSim.multiply(spec, yexNobias, ymanW1, 0, c25over21)
      assert(y25over21ExNobias < 0)

      val onePlus25over21yExNobias = 0
      val onePlus25over21yManW1    = (1<<manW) +
        (y25over21ManW1 >> (-y25over21ExNobias + 2)) +
        bit(-y25over21ExNobias + 2 - 1, y25over21ManW1)

      // ----------------------------------------------------------------------
      // (3/5 * y^2) * (1 + 25/21 * 2^-2 * y)

      val (secondTermExNobias, secondTermManW1) =
        MathFuncACosSim.multiply(spec, ySq3over5ExNobias, ySq3over5ManW1,
          onePlus25over21yExNobias, onePlus25over21yManW1)
      assert(secondTermExNobias < 0)

      // ----------------------------------------------------------------------
      // 1/3 * y

      val c1over3 = math.round(1.0/3.0 * (1<<(manW+2)))
      assert(bit(manW, c1over3) == 1)

      val (yOver3ExNobias, yOver3ManW1) =
        MathFuncACosSim.multiply(spec, yexNobias, ymanW1, -2, c1over3)
      assert(yOver3ExNobias < 0)

      // ----------------------------------------------------------------------
      // 1 + 2^-2 * (1/3 * y) + 2^-5 * (3/5 * y^2) * (1 + 25/21 * 2^-2 * y) < 2
      //                               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      //                               secondTerm

      val puiseuxTermExNobias = 0
      val puiseuxTermManW1    = (1<<manW) +
        (yOver3ManW1     >> (-yOver3ExNobias+2))     + bit(-yOver3ExNobias+2-1,     yOver3ManW1) +
        (secondTermManW1 >> (-secondTermExNobias+5)) + bit(-secondTermExNobias+5-1, secondTermManW1)

      assert(puiseuxTermManW1 < (1<<(manW+1)))

      // ----------------------------------------------------------------------
      // sqrt(2y) * (1 + 2^-2 * (1/3 * y) + 2^-5 * (3/5 * y^2) * (1 + 25/21 * 2^-2 * y))

      val (puiseuxExNobias, puiseuxManW1) =
        MathFuncACosSim.multiply(spec, sqrt2yExNobias,      sqrt2yManW1,
                                       puiseuxTermExNobias, puiseuxTermManW1)
      assert(puiseuxExNobias < 0)

      // ----------------------------------------------------------------------
      // if x < 0, subtract this result from pi.

      if (x.sgn == 1) {
        val piExNobias = 1
        val piManW1    = math.round(Pi * (1<<(manW-piExNobias))).toLong

        val puiseuxAligned = puiseuxManW1 >> (-puiseuxExNobias + 1)
        val puiseuxSub = piManW1 - puiseuxAligned

        // pi = 2 + 1.14. subtracting a value less than 1 never change exponent.
        // we don't need any check.

        val puiseuxEx  = piExNobias + exBias
        val puiseuxMan = slice(0, manW, puiseuxSub)

        return new RealGeneric(spec, zSgn, puiseuxEx, puiseuxMan)

      } else {
        val puiseuxEx  = puiseuxExNobias + exBias
        val puiseuxMan = slice(0, manW, puiseuxManW1)

        return new RealGeneric(spec, zSgn, puiseuxEx, puiseuxMan)
      }

    } else {

//       println(f"taylorThreshold < x. use table.")

      val sgn    = x.sgn
      val exRaw  = x.ex
      val ex     = exRaw-exBias
      val man    = x.man

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

  def calcTaylorThreshold(manW: Int): Int = {
    //   3x^5/40 < 2^-manW
    //       x^5 < 2^-(manW - 3)
    //       x   < 2^-ceil((manW - 3)/5)
    -(math.ceil((manW - 3) / 5.0).toInt)
  }

  // number of tables depending on the exponent and taylorThreshold
  def calcExAdrW(spec: RealSpec): Int = {
    val taylorThreshold = calcTaylorThreshold(spec.manW)
    // from -1 to taylorThreshold (-8 in FP32), 126 to 119 if biased.
    // 0b01111110 to 0b01110111, 4 bits required
    log2Up(abs(taylorThreshold)+1)
  }

  def sqrtTableGeneration( order: Int, adrW: Int, manW: Int, fracW: Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    SqrtSim.sqrtTableGeneration(
      order, adrW, manW, fracW, calcWidthSetting, cbitSetting)
  }

  def acosTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    val taylorThreshold = calcTaylorThreshold(manW)

    if (order == 0 || adrW >= manW) {
      val maxCalcWidth = (-1 to taylorThreshold by -1).map(i => {
        val tableD = new FuncTableDouble( x => acos(scalb(1.0 + x, i)), 0)
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

      (-1 to taylorThreshold by -1).map( i => {
        val tableD = new FuncTableDouble( x => acos(scalb(1.0 + x, i)), 0)
        tableD.addRange(0.0, 1.0, 1<<manW)
        new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
      })
    } else {
      val maxCalcWidth = (-1 to taylorThreshold by -1).map(i => {
        val tableD = new FuncTableDouble(
          x => ((Pi * 0.5) - acos(scalb(1.0 + x, i))) * 0.5, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

      // ex == -1 corresponds to the range [0.5, 1).
      (-1 to taylorThreshold by -1).map( i => {
        val tableD = new FuncTableDouble(
          x => ((Pi * 0.5) - acos(scalb(1.0 + x, i))) * 0.5, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
      })
    }
  }
}
