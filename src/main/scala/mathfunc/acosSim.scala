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
//       println(f"sim: xSqEx    = ${xSqExNobias + exBias }")
//       println(f"sim: xSqManW1 = ${xSqManW1.toLong.toBinaryString}")

      val c1over6ExNobias  = -3
      val c1over6ManW1     = math.round(1.0/6.0 * (1<<(manW+(-c1over6ExNobias))))

      val (xSq6thExNobias, xSq6thManW1) =
        MathFuncACosSim.multiply(spec, xSqExNobias, xSqManW1, c1over6ExNobias, c1over6ManW1)

//       println(f"sim: xSq6thEx    = ${xSq6thExNobias + exBias }")
//       println(f"sim: xSq6thManW1 = ${xSq6thManW1.toLong.toBinaryString}")

      assert(xSq6thExNobias < 0)

      val xSq6thAligned  = xSq6thManW1 >> (-xSq6thExNobias)
      val xSq6thPlusOneExNobias = 0
      val xSq6thPlusOneManW1    = xSq6thAligned + (1<<manW)

//       println(f"sim: xSq6thAligned      = ${xSq6thAligned     .toLong.toBinaryString}%24s")
//       println(f"sim: xSq6thPlusOneManW1 = ${xSq6thPlusOneManW1.toLong.toBinaryString}%24s")

      val (taylorExNobias, taylorManW1) =
        MathFuncACosSim.multiply(spec, xexNobias, xmanW1, xSq6thPlusOneExNobias, xSq6thPlusOneManW1)
      assert(taylorExNobias < 0)

//       println(f"sim: taylorTermEx    = ${taylorExNobias + exBias }")
//       println(f"sim: taylorTermManW1 = ${taylorManW1.toLong.toBinaryString}")

      // ----------------------------------------------------------------------
      // check x.sgn and calculate acos(x) from pi/2 - acos(x)
      //
      // In the corresponding circuit, this part is shared with polynomial that
      // approximates the same formula. Since the precision of polynomial is
      // `extraBits`-bits larger, we need to extend taylor result to avoid
      // rounding error while comparison.
      // The result of taylor expansion is (currently) rounded to `manW`. It is
      // possible to extend the result of Taylor and Puiseux series expansion to
      // fracW. But, to gain the expected precision (<2ULPs), `manW` precision
      // for Taylor and Puiseux series is enough. Normally, the result of those
      // series expansion has small exponent < 0, it might be possible to reduce
      // the precision less than `manW`. But for clarity, we don't do that
      // unless the area and wiring becomes a problem.

      val fracW = ts(0).bp
      val extraBits = fracW - manW

      val halfPiExNobias = 0
      val halfPiManW1    = math.round(Pi * 0.5 * (1<<fracW)).toLong
      assert(bit(fracW, halfPiManW1) == 1)

      val taylorAligned  = (taylorManW1 << (fracW - manW)) >> (-taylorExNobias)
      assert(taylorAligned < halfPiManW1)
//       println(f"sim: halfPiManW1    = ${halfPiManW1.toBinaryString}%24s")
//       println(f"sim: taylorAligned  = ${taylorAligned.toLong.toBinaryString}%24s")

      // Let's say we have FP32. For small x < 2^-4,
      //   pi/2 - acos(x) ~ x + x^3/6 < 2^-4 * (1 + 2^-8 / 6) < 2^-3
      //   1 + 2^-3 = 1.125 < pi/2 ~ 1.5707963 < 2 - 2^-3 = 1.875.
      // we don't need any check. exponent never change.
      val taylorSum = if(x.sgn == 1) {
        halfPiManW1 + taylorAligned
      } else {
        halfPiManW1 - taylorAligned
      }

      val taylorRounded = (taylorSum >> extraBits) + bit(extraBits-1, taylorSum)

      val taylorEx  = halfPiExNobias + exBias
      val taylorMan = slice(0, manW, taylorRounded)
//       println(f"sim: taylorEx  = ${taylorEx}")
//       println(f"sim: taylorMan = ${taylorMan.toLong.toBinaryString}")

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
//       println(f"sim: yex    = ${yexNobias}")
//       println(f"sim: ymanW1 = ${ymanW1.toLong.toBinaryString}")

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
        (y25over21ManW1 >> (-y25over21ExNobias + 2))

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
      // Here, simple rounding cannot be omitted to achieve error < 3ULPs.

//       println(f"sim:  firstTermAligned = ${(yOver3ManW1     >> (-yOver3ExNobias+2)).toLong.toBinaryString}")
//       println(f"sim: secondTermAligned = ${(secondTermManW1 >> (-secondTermExNobias+5)).toLong.toBinaryString}")
//       println(f"sim:  firstTermRounded = ${bit(-yOver3ExNobias+2-1,     yOver3ManW1)}")
//       println(f"sim: secondTermRounded = ${bit(-secondTermExNobias+5-1, secondTermManW1)}")
      assert(puiseuxTermManW1 < (1<<(manW+1)))

      // ----------------------------------------------------------------------
      // sqrt(2y) * (1 + 2^-2 * (1/3 * y) + 2^-5 * (3/5 * y^2) * (1 + 25/21 * 2^-2 * y))

//       println(f"sim: sqrt2y           = ${sqrt2y.toLong.toBinaryString}")
//       println(f"sim: puiseuxTermManW1 = ${puiseuxTermManW1.toLong.toBinaryString}")
      val sqrt2yExNobias   = (yexNobias + 1) >> 1
      val puiseuxProd      = sqrt2y * puiseuxTermManW1
      val puiseuxMoreThan2 = bit((manW+1+manW+extraBits+1)-1, puiseuxProd)
      val puiseuxRounded   = (puiseuxProd >> (manW+extraBits+puiseuxMoreThan2)) +
                             bit((manW+extraBits+puiseuxMoreThan2-1).toInt, puiseuxProd)
      val puiseuxMoreThan2AfterRound = bit(manW+2-1, puiseuxRounded)
      val puiseuxExNobias  = (sqrt2yExNobias + puiseuxTermExNobias +
                             puiseuxMoreThan2 + puiseuxMoreThan2AfterRound).toInt
      val puiseuxManW1     = if(puiseuxMoreThan2AfterRound == 1) {1<<manW} else {puiseuxRounded}
//       println(f"sim: puiseuxTermEx    = ${sqrt2yExNobias + puiseuxTermExNobias + exBias}")
//       println(f"sim: puiseuxTermManW1 = ${puiseuxTermManW1.toLong.toBinaryString}")

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
      val fracW     = manW + extraBits

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
          val piFixed = (math.round(Pi * (1<<fracW))).toLong
          val res0    = piFixed - t.interval(adr).eval(0L, 0)
          // res0 range is [1.57, 3.14]. exponent should be 0 or 1

          val shift = bit(fracW+1, res0)
          val resShifted = res0 >> shift
          val res = resShifted - (1<<fracW)

          (shift.toInt, res)
        } else {
          val res0  = t.interval(adr).eval(0L, 0)
          val shift      = fracW+1 - res0.toLong.toBinaryString.length
          val resShifted = if(shift > 0) {res0 << shift} else {res0}
          val res = resShifted - (1<<fracW)

          (-shift.toInt, res)
        }

      } else { // non-zero order

        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, man).toInt

        val halfPiFixed = math.round(Pi * 0.5 * (1<<fracW))

        // pi/2 - acos(x)
        val res0 = t.interval(adr).eval(d.toLong, dxbp) << 1
        val res  = if (sgn == 1) {
          halfPiFixed + res0
        } else {
          halfPiFixed - res0
        }
        val shift = fracW+2 - res.toLong.toBinaryString.length
        val resShifted = ((res << shift).toLong) >> 1

        ((-shift+1).toInt, resShifted - (1 << fracW))
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
    val acosRequirements = log2Up(abs(taylorThreshold)+1)
    val sqrtRequirements = 1 // does not depends on the spec width

    return max(sqrtRequirements, acosRequirements)
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
        val tableD = new FuncTableDouble( x => {
          val s = scalb(1.0 + x, i) // scaled
          if (1.0 - pow(2.0, -4) < s) {
            0.0
          } else {
            ((Pi * 0.5) - acos(s)) * 0.5
          }
        }, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })
      println(f"acos table width = ${maxCalcWidth}") // 26, 17, 9

      // ex == -1 corresponds to the range [0.5, 1).
      (-1 to taylorThreshold by -1).map( i => {
        val tableD = new FuncTableDouble( x => {
          val s = scalb(1.0 + x, i) // scaled
          if (1.0 - pow(2.0, -4) < s) {
            0.0
          } else {
            ((Pi * 0.5) - acos(s)) * 0.5
          }
        }, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
      })
    }
  }
}
