//% @file atan2Stage2Sim.scala
//
// Simulators for atan2(y, x) stage 2, calculating atan2 from min(x,y)/max(x,y)
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

object ATan2Stage2Sim {

  def atan2Stage2SimGeneric(
    ts : Seq[FuncTableInt], x : RealGeneric, status: Int, special: Int, ysgn: Int
    ): RealGeneric = {

//     println("==================================================")

    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val order     = ts(0).nOrder
    val adrW      = ts(0).adrW
    val fracW     = ts(0).bp
    val extraBits = fracW - manW

//     println(f"x = ${x.sgn}|${x.ex}(${x.ex-exBias})|${x.man.toLong.toBinaryString}")

    // ------------------------------------------------------------------------
    // check special values

    val ysgnUnit = if(ysgn == 0) {1} else {-1}
    if(special == 1) {
      return RealGeneric.nan(x.spec)
    } else if (special == 2) { // zero
      return new RealGeneric(x.spec, ysgnUnit * 0)
    } else if (special == 3) { // pi
      return new RealGeneric(x.spec, ysgnUnit * Pi)
    } else if (special == 4) { // pi/2
      return new RealGeneric(x.spec, ysgnUnit * Pi * 0.5)
    } else if (special == 5) { // pi/4
      return new RealGeneric(x.spec, ysgnUnit * Pi * 0.25)
    } else if (special == 6) { // 3pi/4
      return new RealGeneric(x.spec, ysgnUnit * Pi * 0.75)
    }

    // ------------------------------------------------------------------------
    // atan table
    val linearThreshold = calcLinearThreshold(manW)

//     println(f"sim: xex = ${x.ex}")
//     println(f"sim: linearThreshold = ${linearThreshold}")

    val (atanEx, atanMan) = if(x.ex == 0) {
//       println("atan2Stage2Sim: Less than zero-threshold")
      (x.ex, SafeLong(0))
    } else if(x.ex < linearThreshold + exBias) {
//       println("atan2Stage2Sim: Less than linear-threshold")
      // linear approx
      (x.ex, x.man)
    } else {
      // table interpolation
      val exadr = (-(x.ex - exBias) - 1).toInt
      val t     = ts(exadr)

      val (zEx0, zMan0) = if (t.nOrder == 0) {

        val adr   = x.man
        val res0  = t.interval(adr.toInt).eval(0L, 0)

        val scaling = (-(x.ex - exBias) - 1).toInt
        val shift   = fracW+1 - binaryWidthL(res0)
        val res     = res0 << shift

        (-shift - scaling, SafeLong(res))

      } else {
        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, x.man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, x.man).toInt
        val res0 = t.interval(adr).eval(d.toLong, dxbp)
        assert((SafeLong(1)<<(fracW-2)) <= res0 && res0 < (SafeLong(1)<<fracW))

        val res0MoreThanHalf = bit(fracW-1, res0) == 1
        val shift  = if(res0MoreThanHalf) {1} else {2}

//         println(f"atan2Stage2Sim:              2         1         0")
//         println(f"atan2Stage2Sim:          4321098765432109876543210")
//         println(f"atan2Stage2Sim: res0   = ${res0.toLong.toBinaryString}")
//         println(f"atan2Stage2Sim: res>1/2= ${res0MoreThanHalf}")
//         println(f"atan2Stage2Sim: fracW  = ${fracW}")
//         println(f"atan2Stage2Sim: shift  = ${shift}")
//         println(f"atan2Stage2Sim: x.ex   = ${x.ex - exBias}")
//         println(f"atan2Stage2Sim: -x.ex-1= ${-(x.ex - exBias) - 1}")
//         println(f"atan2Stage2Sim: zex0   = ${-(x.ex - exBias) - 1 - shift}")

        val xexNobias = x.ex - exBias
        val zex0 = xexNobias + 1 - shift
        val res  = SafeLong(res0) << shift // normalize

//         println(f"atan2Stage2Sim:               2         1         0")
//         println(f"atan2Stage2Sim:          54321098765432109876543210")
//         println(f"atan2Stage2Sim: res    = ${res.toLong.toBinaryString}")

        assert((SafeLong(1)<<fracW) <= res && res <= (SafeLong(1)<<(fracW+1)))
        (zex0, res)
      }

      // round the result from table
      val zMan = if (extraBits > 0) {
        (zMan0 >> extraBits) + bit(extraBits-1, zMan0)
      } else {
        zMan0
      }

//       println(f"atan2Stage2Sim:             2         1         0")
//       println(f"atan2Stage2Sim:          321098765432109876543210")
//       println(f"atan2Stage2Sim: zMan   = ${zMan.toLong.toBinaryString}")

      assert((1<<manW) <= zMan)

      (zEx0+exBias, zMan - (SafeLong(1)<<manW))
    }
    assert(0 <= atanMan && atanMan <= (SafeLong(1)<<manW))

//     println(f"atan2Stage2Sim: atanEx  = ${atanEx }(${atanEx - exBias})")
//     println(f"atan2Stage2Sim: atanMan = ${atanMan.toLong.toBinaryString}(${atanMan})")
//     println(f"atan2Stage2Sim: atan(x) = sim(${new RealGeneric(x.spec, ysgn, atanEx, atanMan).toDouble}) == ref(${atan(x.toDouble)})")
//     println(f"atan2Stage2Sim: status  = ${status}")

    // z correction by:
    //   ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y| : 0
    //   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y| : 1
    //   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y| : 2
    //   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y| : 3

    val zsgn = ysgn

    val piManW1     = Real.pi(manW+2)
    val halfpiManW1 = (Real.pi / Real.two)(manW+2)
    val atanManW1   = (atanMan + (SafeLong(1)<<manW)) << 2

    // atan(min(|x|,|y|)/max(|x|,|y|)) < pi/4
    assert(0 <= atanEx && atanEx < exBias)
    val atanShift   = exBias - atanEx
    val atanAligned = atanManW1 >> atanShift

//     println(f"atan2Stage2Sim: piManW1     = ${piManW1    .toLong.toBinaryString}")
//     println(f"atan2Stage2Sim: halfpiManW1 =  ${halfpiManW1.toLong.toBinaryString}")
//     println(f"atan2Stage2Sim: atanManW1   =  ${atanManW1  .toLong.toBinaryString}")
//     println(f"atan2Stage2Sim: atanShift   = ${atanShift  }")
//     println(f"atan2Stage2Sim: atanAligned = ${atanAligned.toLong.toBinaryString}")

//     println(f"atanShift = ${atanShift}")

    if (status == 0) {//   ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|

      return new RealGeneric(x.spec, zsgn, atanEx, atanMan)

    } else if (status == 1) {//   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y| : 1

      // since atanAligned is in [0, pi/4]  pi-atan is in [3/4pi, pi].
      // 3/4pi ~ 2.35, exponent of pi - atan should always be 1.
      val zman0 = piManW1 - atanAligned
      val zman = roundBySpec(RoundSpec.round, 3, zman0)

//       println(f"atan2Stage2Sim: pi  .man = ${piManW1.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: atan.man = ${atanAligned.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: sub .man = ${zman0.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: z   .man = ${zman.toLong.toBinaryString}")

      return new RealGeneric(x.spec, zsgn, 1 + exBias, zman - (SafeLong(1)<<manW))

    } else if (status == 2) {//   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y| : 2

      // since atanAligned is in [0, pi/4]  pi/2-atan is in [pi/4, pi/2].
      // z is in [0.78 ~ 1.57).
      val zman0 = halfpiManW1 - atanAligned

      val zman0LessThan1 = if(bit(manW+2, zman0) == 0) { 1 } else { 0 }
      val zmanRound = roundBySpec(RoundSpec.round, 2-zman0LessThan1, zman0)
      val zman = zmanRound - (SafeLong(1)<<manW)
      val zex  = exBias - zman0LessThan1

//       println(f"atan2Stage2Sim: pi/2.man = ${halfpiManW1.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: atan.man = ${atanAligned.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: sub .man = ${zman0.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: z   .man = ${zman.toLong.toBinaryString}")

      return new RealGeneric(x.spec, zsgn, zex, zman)

    } else                  {//   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y| : 3
      assert(status == 3)
      // since atanAligned is in [0, pi/4]  pi/2-atan is in [pi/2, 3pi/4].
      // z is in [1.57, 2.35).
      val zman0 = halfpiManW1 + atanAligned

      val zman0MoreThan2 = bit(manW+1+3-1, zman0)
      val zmanRound = roundBySpec(RoundSpec.round, 2+zman0MoreThan2, zman0)
      val zman = zmanRound - (SafeLong(1)<<manW)
      val zex  = zman0MoreThan2 + exBias

//       println(f"atan2Stage2Sim:                 2         1         0")
//       println(f"atan2Stage2Sim:            54321098765432109876543210")
//       println(f"atan2Stage2Sim: pi/2.man = ${halfpiManW1.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: atan.man = ${atanAligned.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: add .man = ${zman0.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: rounded  = ${zmanRound.toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: 1<<manW  = ${(1<<manW).toLong.toBinaryString}")
//       println(f"atan2Stage2Sim: z   .man = ${zman.toLong.toBinaryString}") // ?

      return new RealGeneric(x.spec, zsgn, zex, zman)
    }
  }

  def calcLinearThreshold(manW: Int): Int = {
    -math.round(math.ceil(manW / 2.0 + 1.0)).toInt
  }

  // number of tables depending on the exponent and linearThreshold
  def calcExAdrW(spec: RealSpec): Int = {
    val linearThreshold = calcLinearThreshold(spec.manW)
    log2Up(abs(linearThreshold)+1)
  }

  def atanTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {

    // atan(x) = x - x^3/3 + x^5/5 + O(x^7)
    //
    // If x is in [0, 1],
    //                  x/2 < atan(x)          < x
    //         2^ex-1 * 1.m < atan(x)          < 2^ex * 1.m
    //   1/4 < 2^-2   * 1.m < 2^(-ex-1)atan(x) < 1.m / 2 < 1
    //        0.01 ==_2 1/4 < 2^(-ex-1)atan(x) < 1
    //                        ~~~~~~~~~~~~~~~~
    //                        this is calculated in polynomial
    //
    // So the result may take 0.010000 ~ 0.111111.

    val linearThreshold = calcLinearThreshold(manW)

    val nOrder = if (adrW >= manW) { 0 } else { order }

    val maxCalcWidth = (-1 to linearThreshold by -1).map(ex => {
        val tableD = new FuncTableDouble(
          x => scalb(atan(scalb(1.0 + x, ex)), -ex-1),
          nOrder)
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

    (-1 to linearThreshold by -1).map( ex => {
      val tableD = new FuncTableDouble(
        x => scalb(atan(scalb(1.0 + x, ex)), -ex-1),
        nOrder)
      tableD.addRange(0.0, 1.0, 1<<adrW)
      new FuncTableInt( tableD, fracW, Some(maxCalcWidth), cbitSetting )
    })
  }

}
