//% @file atan2Phase2Sim.scala
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

private[rial] object ATan2Phase2Sim {

  // Since atan2-stage1 calculates min(x,y)/max(x,y), so it assumes x <= 1.
  def atan2Phase2SimGeneric(t : FuncTableInt, x : RealGeneric): RealGeneric = {

//     println("==================================================")

    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val order     = t.nOrder
    val adrW      = t.adrW
    val fracW     = t.bp
    val extraBits = fracW - manW

//     println(f"x = ${x.sgn}|${x.ex}(${x.ex-exBias})|${x.man.toLong.toBinaryString}")

    // ------------------------------------------------------------------------
    // check special values

    if(x.ex == maskI(exW)) { // NaN Boxing
      val sgnCoef = if(x.sgn == 0) {1} else {-1}

      if(x.man == 1) { // z == 0
        return new RealGeneric(x.spec, x.sgn, 0, 0)
      } else if(x.man == 2) { // z == pi
        return new RealGeneric(x.spec, sgnCoef * Pi)
      } else if(x.man == 3) { // z == pi/2
        return new RealGeneric(x.spec, sgnCoef * Pi * 0.5)
      } else if(x.man == 4) { // z == pi/4
        return new RealGeneric(x.spec, sgnCoef * Pi * 0.25)
      } else if(x.man == 5) { // z == 3pi/4
        return new RealGeneric(x.spec, sgnCoef * Pi * 0.75)
      } else {
        return RealGeneric.nan(x.spec)
      }
    }

    // ------------------------------------------------------------------------
    // decode states

    val xsgn = 0
    val xex  = x.ex  & ((1 << (exW-1)) - 1) // remove the msb
    val xman = x.man & ((1 <<    manW) - 2) // ignore the last bit (set zero)
    val xmanW1 = xman + (1 << manW)

    val status = (((x.ex >> (exW-1)) & 1) << 1) | (x.man & 1).toInt
    val ysgn   = x.sgn

    assert((new RealGeneric(x.spec, xsgn, xex, xman)).toDouble <= 1.0)

    // ------------------------------------------------------------------------
    // atan(x)/x table

    val xShift0 = (exBias - xex)
    val xShift  = Seq(xShift0, manW+1).min
    val xFixed  = ( xmanW1 | 1 ) >> xShift

    // The case where x == 1, that means that x == y in atan2(x,y), corresponds
    // to the 5- or 6-th special case.
    assert(0 < xShift0, f"xShift = ${xShift}")

    val (resManW1, resEx) = if(xex <= exBias + calcLinearThreshold(manW)) {

      // without this check, polynomial become larger than (not equal to) 1...
      (SafeLong(1) << fracW, 0)

    } else {

      val res0 = if (t.nOrder == 0) {
        val adr   = xFixed
        t.interval(adr.toInt).eval(0L, 0)
      } else {
        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, xFixed) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, xFixed).toInt
        t.interval(adr).eval(d.toLong, dxbp)
      }

      assert(res0 <= (SafeLong(1) << fracW),
             f"xFixed = ${xFixed}, res0 = ${res0} > (1<<fracW) = ${SafeLong(1)<<fracW}")

      (SafeLong(res0) << 1, -1)
    }

    // ------------------------------------------------------------------------
    // calc x * atan(x)/x
    //          ^^^^^^^^^ this term is calculated by polynomial

    val atanProd = resManW1 * xmanW1
    val atanProdMoreThan2 = bit((fracW+1)+(manW+1)-1, atanProd).toInt
    val atanRoundBit = fracW + atanProdMoreThan2
    val atanRound = roundBySpec(RoundSpec.roundToEven, atanRoundBit, atanProd)
    val atanRoundMoreThan2 = bit(manW+1, atanRound).toInt

    val atanMan = slice(0, manW, atanRound >> atanRoundMoreThan2)
    val atanEx  = xex + resEx + atanProdMoreThan2 + atanRoundMoreThan2

//     println(f"atan2Phase2Sim: atanEx  = ${atanEx }(${atanEx - exBias})")
//     println(f"atan2Phase2Sim: atanMan = ${atanMan.toLong.toBinaryString}(${atanMan})")
//     println(f"atan2Phase2Sim: atan(x) = sim(${new RealGeneric(x.spec, ysgn, atanEx, atanMan).toDouble}) == ref(${atan(x.toDouble)})")
//     println(f"atan2Phase2Sim: status  = ${status}")

    // ========================================================================
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

//     println(f"atan2Phase2Sim: piManW1     = ${piManW1    .toLong.toBinaryString}")
//     println(f"atan2Phase2Sim: halfpiManW1 =  ${halfpiManW1.toLong.toBinaryString}")
//     println(f"atan2Phase2Sim: atanManW1   =  ${atanManW1  .toLong.toBinaryString}")
//     println(f"atan2Phase2Sim: atanShift   = ${atanShift  }")
//     println(f"atan2Phase2Sim: atanAligned = ${atanAligned.toLong.toBinaryString}")

//     println(f"atanShift = ${atanShift}")

    if (status == 0) {//   ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|

      return new RealGeneric(x.spec, zsgn, atanEx, atanMan)

    } else if (status == 1) {//   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y| : 1

      // since atanAligned is in [0, pi/4]  pi-atan is in [3/4pi, pi].
      // 3/4pi ~ 2.35, exponent of pi - atan should always be 1.
      val zman0 = piManW1 - atanAligned
      val zman = roundBySpec(RoundSpec.round, 3, zman0)

//       println(f"atan2Phase2Sim: pi  .man = ${piManW1.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: atan.man = ${atanAligned.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: sub .man = ${zman0.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: z   .man = ${zman.toLong.toBinaryString}")

      return new RealGeneric(x.spec, zsgn, 1 + exBias, zman - (SafeLong(1)<<manW))

    } else if (status == 2) {//   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y| : 2

      // since atanAligned is in [0, pi/4]  pi/2-atan is in [pi/4, pi/2].
      // z is in [0.78 ~ 1.57).
      val zman0 = halfpiManW1 - atanAligned

      val zman0LessThan1 = if(bit(manW+2, zman0) == 0) { 1 } else { 0 }
      val zmanRound = roundBySpec(RoundSpec.round, 2-zman0LessThan1, zman0)
      val zman = zmanRound - (SafeLong(1)<<manW)
      val zex  = exBias - zman0LessThan1

//       println(f"atan2Phase2Sim: pi/2.man = ${halfpiManW1.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: atan.man = ${atanAligned.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: sub .man = ${zman0.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: z   .man = ${zman.toLong.toBinaryString}")

      return new RealGeneric(x.spec, zsgn, zex, zman)

    } else                  {//   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y| : 3
      assert(status == 3)

      // since atanAligned is in [0, pi/4]  pi/2-atan is in [pi/2, 3pi/4].
      // z is in [1.57, 2.35).
      val zman0 = halfpiManW1 + atanAligned
      val zmanRound = roundBySpec(RoundSpec.round, 2, zman0)
      val zmanMoreThan2 = bit(manW+1, zmanRound)
      val zman = slice(zmanMoreThan2, manW, zmanRound)
      val zex  = zmanMoreThan2 + exBias

//       println(f"atan2Phase2Sim: pi/2.man = ${halfpiManW1.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: atan.man = ${atanAligned.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: add .man = ${zman0.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: rounded  = ${zmanRound.toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: 1<<manW  = ${(1<<manW).toLong.toBinaryString}")
//       println(f"atan2Phase2Sim: z   .man = ${zman.toLong.toBinaryString}") // ?

      return new RealGeneric(x.spec, zsgn, zex, zman)
    }
  }

  def calcLinearThreshold(manW: Int): Int = {
    -math.ceil(manW / 2.0 + 1.0).toInt
  }

  def atanTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    val tableD = new FuncTableDouble( x => {
      val eps = pow(2.0, -manW)
      val res = if(x == 0) {
        1.0 - eps
      } else {
        Seq(atan(x) / x, 1.0 - eps).min
      }
      assert(0.0 < res && res < 1) // [pi/4, 1)
      res
    }, order)
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt(tableD, fracW, calcWidthSetting, cbitSetting)
  }

}
