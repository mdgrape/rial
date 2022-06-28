//% @file acosSim.scala
//
// Simulators for acos(x) function
// Copyright (C) Toru Niina RIKEN BDR 2022
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

object ACosSim {

  // assuming x, y, and z have the same spec.
  // since acos uses 2 different high-order series expansion,
  // this kind of helper function is useful.
  def multiply(spec: RealSpec,
    xexNobias: Int, xmanW1: SafeLong, yexNobias: Int, ymanW1: SafeLong
    ): (Int, SafeLong) = {
    val manW = spec.manW

    val zProd      = xmanW1 * ymanW1
    val zMoreThan2 = bit(((manW+1)*2-1).toInt, zProd)
    val zRounded   = (zProd >> (manW+zMoreThan2)) +
                     bit((manW+zMoreThan2-1).toInt, zProd)
    val zMoreThan2AfterRound = bit(manW+2-1, zRounded)
    val zExNobias  = xexNobias + yexNobias + zMoreThan2 + zMoreThan2AfterRound
    val zManW1     = if(zMoreThan2AfterRound == 1) {SafeLong(1) << manW} else {zRounded}

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
  // for x < 0.5,  pi/2 - acos(x)      is approximated by polynomial.
  // for x > 0.5, (pi/2 - acos(x)) / 2 is approximated by polynomial.
  // for x close to 1, use puiseux series:
  //   acos(1-x) = sqrt(2x) * (1 + x/12 + 3x^2/160 + 5x^3/896 + 35x^4/18432 + O(x^5))
  //
  // In case of Puiseux series, the condition is:
  //   35x^4/18432 < 2^-23
  //     x^4/526.628... < 2^-23
  //     x^4 < 2^-14
  //     x   < 2^-4
  // So if 1-2^-4 < x, that means that x.man(22,20).andR === 1, use Puiseux series.
  //
  def acosSimGeneric(
      tACos: Seq[FuncTableInt], tSqrt: FuncTableInt, x: RealGeneric,
      // if order == 0, exponent will be calculated by a table.
      exTable: Option[Seq[FuncTableInt]] = None
    ) : RealGeneric = {

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val extraBits = tACos(0).bp - manW
    assert(extraBits == tSqrt.bp - manW)

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

    assert(x.ex < exBias)

    val puiseuxThreshold = calcPuiseuxThreshold(manW, extraBits, 0)
    // 1 for hidden bit, 1 for table exponent
    val puiseuxMSBs = puiseuxThreshold.abs - 1 - 1
    if(tACos(0).nOrder != 0 && x.ex == exBias-1 &&
       slice(manW-puiseuxMSBs, puiseuxMSBs, x.man) == maskSL(puiseuxMSBs)) {
//       println(f"use Puiseux series: |x| = ${x.toDouble.abs}%10f = ${x.manW1.toLong.toBinaryString}, y = ${1.0 - x.toDouble.abs}")

      // for x close to 1, use puiseux series:
      //   let y = 1 - x
      //   acos(1-y) = sqrt(2y) * (1 + y/12 + 3y^2/160 + 5y^3/896 + O(y^4))
      //             = sqrt(2y) * ((1 + 1/3 * 2^-2 * y) + 3y^2/160 * (1 + 25/21 * 2^-2 * y))

      val oneMinusX    = (SafeLong(1)<<(manW+1)) - x.manW1
      // to use sqrt table, we need to normalize 1-x.
      val oneMinusXLen = binaryWidthSL(oneMinusX)
      val ymanW1    = (oneMinusX << (manW+1 - oneMinusXLen))
      val yman      = ymanW1 - (SafeLong(1)<<manW)
      val yexNobias = -(manW+1 - oneMinusXLen)-1 // -1 because oneMinusX has 1bit wider mantissa
      val yex       = yexNobias + exBias

      assert(bit(manW, ymanW1) == 1)
      assert(yman < (SafeLong(1)<<manW))
      assert(0 < yex)
//       println(f"sim: yex    = ${yexNobias}")
//       println(f"sim: ymanW1 = ${ymanW1.toLong.toBinaryString}")

//       println(f"yref = ${1.0 - x.toDouble.abs}")
//       println(f"ysim = ${ymanW1.toDouble / (1<<manW) * pow(2.0, yexNobias)}")

      // ----------------------------------------------------------------------
      // 2^-5 * 3/5 y^2

      val c3over5 = (SafeLong(3)<<(manW+1)) / SafeLong(5) // +1 for normalize
      assert(bit(manW, c3over5) == 1)

      // in the circuit, the width of oneMinusX should be manW+1 because
      // it is the result of (1<<manW+1) - x.manW1. And it is impossible to
      // change the width dynamically, so here we can use ymanW1 without any
      // additional cost.

      val (ySqExNobias, ySqManW1) =
        ACosSim.multiply(spec, yexNobias, ymanW1, yexNobias, ymanW1)

      val (ySq3over5ExNobias, ySq3over5ManW1) =
        ACosSim.multiply(spec, ySqExNobias, ySqManW1, -1, c3over5)

      // ----------------------------------------------------------------------
      // 1 + 25/21 * 2^-2 * y
      //
      // here, y < 2^-4. so 1 + 25/21 * 2^-2 * y never exceeds 2.
      // We don't need to check if it exceeds 2 after addition.

//       val c25over21 = math.round(25.0/21.0 * (SafeLong(1)<<manW))
      val c25over21 = (SafeLong(25)<<manW) / SafeLong(21)
      assert(bit(manW, c25over21) == 1)

      val (y25over21ExNobias, y25over21ManW1) =
        ACosSim.multiply(spec, yexNobias, ymanW1, 0, c25over21)
      assert(y25over21ExNobias < 0)

      val onePlus25over21yExNobias = 0
      val onePlus25over21yManW1    = (SafeLong(1)<<manW) +
        (y25over21ManW1 >> (-y25over21ExNobias + 2))

      // ----------------------------------------------------------------------
      // (3/5 * y^2) * (1 + 25/21 * 2^-2 * y)

      val (secondTermExNobias, secondTermManW1) =
        ACosSim.multiply(spec, ySq3over5ExNobias, ySq3over5ManW1,
          onePlus25over21yExNobias, onePlus25over21yManW1)
      assert(secondTermExNobias < 0)

      // ----------------------------------------------------------------------
      // 1/3 * y

      val c1over3 = (SafeLong(1)<<(manW+2)) / SafeLong(3)
      assert(bit(manW, c1over3) == 1)

      val (yOver3ExNobias, yOver3ManW1) =
        ACosSim.multiply(spec, yexNobias, ymanW1, -2, c1over3)
      assert(yOver3ExNobias < 0)

      // ----------------------------------------------------------------------
      // 1 + 2^-2 * (1/3 * y) + 2^-5 * (3/5 * y^2) * (1 + 25/21 * 2^-2 * y) < 2
      //                               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      //                               secondTerm

      val puiseuxTermExNobias = 0
      val puiseuxTermManW1    = (SafeLong(1)<<manW) +
        (yOver3ManW1     >> (-yOver3ExNobias+2))     + bit(-yOver3ExNobias+2-1,     yOver3ManW1) +
        (secondTermManW1 >> (-secondTermExNobias+5)) + bit(-secondTermExNobias+5-1, secondTermManW1)
      // Here, simple rounding cannot be omitted to achieve error < 3ULPs.

//       println(f"sim:  firstTermAligned = ${(yOver3ManW1     >> (-yOver3ExNobias+2)).toLong.toBinaryString}")
//       println(f"sim: secondTermAligned = ${(secondTermManW1 >> (-secondTermExNobias+5)).toLong.toBinaryString}")
//       println(f"sim:  firstTermRounded = ${bit(-yOver3ExNobias+2-1,     yOver3ManW1)}")
//       println(f"sim: secondTermRounded = ${bit(-secondTermExNobias+5-1, secondTermManW1)}")
      assert(puiseuxTermManW1 < (SafeLong(1)<<(manW+1)))

      // ----------------------------------------------------------------------
      // sqrt(2y)
      // sqrt table uses the last bit of exponent.

      val adrW      = tSqrt.adrW
      val extraBits = tSqrt.bp - manW
      val y2man  = slice(0, manW+1, ((SafeLong(yex) + SafeLong(1))<<manW) + yman)

      val dxbp   = (manW+1)-adrW-1
      val d      = slice(0, (manW+1)-adrW, y2man) - (SafeLong(1)<<dxbp)
      val adr    = slice((manW+1)-adrW, adrW, y2man).toInt

      // sqrt table returns the mantissa only. we need to add 1<<fracW
      val sqrt2y = tSqrt.interval(adr).eval(d.toLong, dxbp) + (SafeLong(1)<<(tSqrt.bp))

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
      val puiseuxManW1     = if(puiseuxMoreThan2AfterRound == 1) {
        SafeLong(1) << manW
      } else {puiseuxRounded}
//       println(f"sim: puiseuxTermEx    = ${sqrt2yExNobias + puiseuxTermExNobias + exBias}")
//       println(f"sim: puiseuxTermManW1 = ${puiseuxTermManW1.toLong.toBinaryString}")

      assert(puiseuxExNobias < 0)

      // ----------------------------------------------------------------------
      // if x < 0, subtract this result from pi.

      if (x.sgn == 1) {
        val piExNobias = 1
        val piManW1    = Real.pi(manW-piExNobias)
//           math.round(Pi * (SafeLong(1)<<(manW-piExNobias))).toLong

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

    } else if (x.ex <= exBias - 2) {
//       println(f"xex < -1. x(${x.toDouble} = ${x.manW1.toLong.toBinaryString}) is in [0.0, 0.5). use acos small x table")

      val xmanW1 = x.manW1
      val dex = exBias - 2 - x.ex + 1
      val xAligned = xmanW1 >> dex

//       println(f"x        = ${x.sgn}|${x.ex}(${x.ex-exBias})|${x.man.toLong.toBinaryString}(${x.toDouble})")
//       println(f"xmanW1   = ${xmanW1.toLong.toBinaryString} (${xmanW1.toDouble / (1<<manW)})")
//       println(f"xAligned = ${xAligned.toLong.toBinaryString} (${xAligned.toDouble / (1<<manW)})")

      val t = tACos(0)

      val adrW      = t.adrW
      val nOrder    = t.nOrder
      val fracW     = t.bp

      val order =
        if (adrW>=manW) {
          if (nOrder != 0)
            println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
          0
        } else {
          nOrder
        }

      val (zEx, zMan) = if (order == 0) {
        val adr = xAligned.toInt

        assert(exTable.isDefined, "if order == 0, exTable should be provided")
        val exts = exTable.get

        val zex0  = exts(0).interval(adr).eval(0L, 0)
        val zman0 = t.interval(adr).eval(0L, 0)

        if(x.sgn == 0) {
          (zex0.toInt - exBias, SafeLong(zman0))
        } else {
          val piEx    = 1
//           val piFixed = math.round(math.Pi * (SafeLong(1)<<(fracW-1))).toSafeLong
          val piFixed = Real.pi(fracW-1)
          val zmanShift = piEx - (zex0 - exBias)
          val zmanAligned = if(zmanShift > fracW) {SafeLong(0)} else {
            (zman0+(SafeLong(1)<<fracW)) >> zmanShift.toInt
          }
          val res0 = piFixed - zmanAligned // > 1.57.., ex = 1 (in [2, 4))
          val res0MoreThan2 = bit(fracW, res0)

          if(res0MoreThan2 == 1) {
            (1, res0 - (SafeLong(1)<<fracW))
          } else {
            (0, (res0 << 1) - (SafeLong(1)<<fracW))
          }
        }
      } else { // non-zero order

        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, xAligned) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, xAligned).toInt

        // pi/2 - acos(x)
        val polynomial = t.interval(adr).eval(d.toLong, dxbp).toSafeLong
        assert(0 <= polynomial && polynomial < (SafeLong(1)<<fracW),
          f"polynomial = ${polynomial} > (1<<fracW) = ${SafeLong(1)<<fracW}, fracW = ${fracW}")

        val halfPiFixed = (Real.pi / Real.two)(fracW)
        // math.round(Pi * 0.5 * (SafeLong(1)<<fracW)).toSafeLong

        val res  = if (x.sgn == 1) {
          halfPiFixed + polynomial
        } else {
          halfPiFixed - polynomial
        }
        val shift = fracW+2 - binaryWidthSL(res)
        val resShifted = (res << shift) >> 1

        assert(resShifted > (SafeLong(1)<<fracW))

        ((-shift+1).toInt, resShifted - (SafeLong(1) << fracW))
      }

      val zmanRound = if (extraBits>0) {
        (zMan >> extraBits) + bit(extraBits-1, zMan)
      } else {zMan}

      val zmanMoreThan2AfterRound = bit(manW, zmanRound)
      val zman = slice(0, manW, zmanRound)
      val zex  = zEx + zmanMoreThan2AfterRound + exBias

      return new RealGeneric(x.spec, zSgn, zex.toInt, zman)

    } else {
//       println(f"x(${x.toDouble}%10f = ${x.manW1.toLong.toBinaryString}) is in [0.5, 1.0]; acos(x) = ${acos(x.toDouble)}. use acos large x table.")

      val t = tACos(1)

      val adrW      = t.adrW
      val nOrder    = t.nOrder
      val fracW     = t.bp

      val order =
        if (adrW>=manW) {
          if (nOrder != 0)
            println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
          0
        } else {
          nOrder
        }

      val (zEx, zMan) = if (order == 0) {
        // XXX if order == 0, all the possible value should be tabulized.
        //     But the volume of memory in the real world is limited.
        //     We expect that there should not be 2^32 entry table.
        val adr = x.man.toInt

        assert(exTable.isDefined, "if order == 0, exTable should be provided")
        val exts = exTable.get

        val zex0  = exts(1).interval(adr).eval(0L, 0)
        val zman0 = t.interval(adr).eval(0L, 0).toSafeLong

        if(x.sgn == 0) {
          (zex0.toInt - exBias, zman0)
        } else {
          val piEx    = 1
//           val piFixed = math.round(math.Pi * (1<<(fracW-1))).toLong
          val piFixed = Real.pi(fracW-1)
          val zmanShift = piEx - (zex0 - exBias)
          val zmanAligned = if(zmanShift > fracW) {SafeLong(0)} else {
            (zman0 + (SafeLong(1)<<fracW)) >> zmanShift.toInt
          }
          val res0 = piFixed - zmanAligned // > 1.57.., ex = 1 (in [2, 4))
          val res0MoreThan2 = bit(fracW, res0)

          if(res0MoreThan2 == 1) {
            (1, res0 - (1<<fracW))
          } else {
            (0, (res0 << 1) - (1<<fracW))
          }
        }
      } else { // non-zero order

        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, x.man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, x.man).toInt

//         val halfPiFixed = math.round(Pi * 0.5 * (1<<fracW))
        val halfPiFixed = (Real.pi / Real.two)(fracW)

        // pi/2 - acos(x)
        val polynomial = t.interval(adr).eval(d.toLong, dxbp).toSafeLong
        assert(0 <= polynomial && polynomial < (SafeLong(1)<<fracW),
          f"polynomial = ${polynomial} > (1<<fracW) = ${SafeLong(1)<<fracW}, fracW = ${fracW}")

        val res0 = polynomial << 1
        val res  = if (x.sgn == 1) {
          halfPiFixed + res0
        } else {
          halfPiFixed - res0
        }
        val shift = fracW+2 - binaryWidthSL(res)
        val resShifted = (res << shift) >> 1

        ((-shift+1).toInt, resShifted - (SafeLong(1) << fracW))
      }
//       println(f"zex    = ${zEx}")
//       println(f"zman   = ${zman.toLong.toBinaryString}")

      val zmanRound = if (extraBits>0) {(zMan >> extraBits) + bit(extraBits-1, zMan)} else {zMan}

//       println(f"zmanR  = ${zmanRound.toLong.toBinaryString}")

      val zmanMoreThan2AfterRound = bit(manW, zmanRound)
      val zman = slice(0, manW, zmanRound)
      val zex  = zEx + zmanMoreThan2AfterRound + exBias
//       println(f"Sim: zman = ${z.toBinaryString}")

      return new RealGeneric(x.spec, zSgn, zex.toInt, zman)
    }
  }

  def calcExAdrW(spec: RealSpec): Int = {
    return 1
  }

  // The puiseux series should fill the gap between acos(1-x) to acos(1) where
  // the corresponding table cannot calculate with sufficient precision.
  // The table for x in [0.5, 1) can calculate acos(x) with manW + extraBits-1
  // bits. The reason why the bit precision is subtracted by 1 because
  // acos(0.5) > 1 so the resulting value is divided by 2. So, if the table
  // result becomes less than 2^-(extraBits-1), then the precision becomes
  // insufficient.
  //                                                    prec
  // cos(1.0)  = 0.54030.. <=> acos(0.54030..) = 1.0    manW + extraBits - 1
  // cos(0.5)  = 0.87758.. <=> acos(0.87758..) = 0.5    manW + extraBits - 2
  // cos(0.25) = 0.96891.. <=> acos(0.96891..) = 0.25   manW + extraBits - 3
  //
  // and
  //
  // 1 - 1/2^1 = 0.1_(2)     = 0.5      < 0.54030.. = cos(1.0)
  // 1 - 1/2^3 = 0.111_(2)   = 0.875    < 0.87758.. = cos(0.5)
  // 1 - 1/2^5 = 0.11111_(2) = 0.96875  < 0.96891.. = cos(0.25)
  //
  // So, if
  //   extraBits == 3, table can fill [0.5, 0.96891..).
  //   extraBits == 2, table can fill [0.5, 0.87758..).
  //   extraBits == 1, table can fill [0.5, 0.54030..).
  //   extraBits == 0, table can not fill anywhere with precision >= manW.
  // generally,
  //   extraBits == n, table can fill [0.5, cos 2^-(n-1) ).
  //
  def calcPuiseuxThreshold(manW: Int, extraBits: Int, tolerance: Int = 0): Int = {
    val log2 = (a:Double) => {log(a) / log(2.0)}
    val threshold = 1.0 - cos(pow(2.0, -(extraBits - 1 + tolerance)))
    math.ceil(log2(threshold)).toInt
  }

  def sqrtTableGeneration( order: Int, adrW: Int, manW: Int, fracW: Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    SqrtSim.sqrtTableGeneration(
      order, adrW, manW, fracW, calcWidthSetting, cbitSetting)
  }

  def acosTableGeneration(spec: RealSpec,
      order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {

    assert(spec.manW == manW)
    val puiseuxThreshold = calcPuiseuxThreshold(manW, fracW-manW, 0)

    if (order == 0 || adrW >= manW) {
      // 0 ~ 0.5 (ex <= -2)
      val f = (x: Double) => {
        val z = new RealGeneric(spec, acos(x * 0.5))
        z.man.toDouble / (1<<manW)
      }
      // 0.5 ~ 1 (ex == -1)
      val g = (x:Double) => {
        val s = (1.0 + x) * 0.5
        val z = new RealGeneric(spec, acos(s))
        z.man.toDouble / (1<<manW)
      }

      val maxCalcWidth = Seq(f, g).map(func => {
        val tableD = new FuncTableDouble( func, 0 )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

      Seq(f, g).map(func => {
        val tableD = new FuncTableDouble( func, 0 )
        tableD.addRange(0.0, 1.0, 1<<manW)
        new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
      })
    } else {

      // 0 ~ 0.5 (ex <= -2).
      val f = (x: Double) => {
        // map arg in [0, 1] into pi/2 - [acos(0), acos(0.5)]
        math.Pi * 0.5 - acos(x * 0.5)
      }
      // 0.5 ~ 1 (ex == -1)
      val g = (x:Double) => {
        val s = (1.0 + x) * pow(2.0, -1) // scaled 0~1 to 0.5~1
        if (1.0 - pow(2.0, puiseuxThreshold) < s) { // puiseux threshold
          0.0 // to make >2nd order coefficient smaller
        } else {
          (math.Pi * 0.5 - acos(s)) * 0.5
        }
      }

      val maxCalcWidth = Seq(f, g).map(func => {
        val tableD = new FuncTableDouble( func, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

      // ex == -1 corresponds to the range [0.5, 1).
      Seq(f, g).map( func => {
        val tableD = new FuncTableDouble( func, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
      })
    }
  }

  def acosExTableGeneration(spec: RealSpec,
      order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
      assert(order == 0)

      val exW = spec.exW

      // 0 ~ 0.5 (ex <= -2).
      val f = (x: Double) => {
        val z = new RealGeneric(spec, acos(x * 0.5))
        z.ex.toDouble / (1 << exW)
      }
      // 0.5 ~ 1.0 (ex == -1).
      val g = (x: Double) => {
        val s = (1.0 + x) * 0.5 // 0~1 => 0.5~1
        val z = new RealGeneric(spec, acos(s))
        z.ex.toDouble / (1 << exW)
      }

      val maxCalcWidth = Seq(f, g).map(func => {
        val tableD = new FuncTableDouble( func, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

      Seq(f, g).map( func => {
        val tableD = new FuncTableDouble( func, order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
      })
  }
}
