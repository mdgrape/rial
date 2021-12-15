//% @file atan2Sim.scala
//
// Simulators for atan2(y, x) stage 1, calculating min(x,y)/max(x,y)
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

object ATan2Stage1Sim {

  def atan2Stage1SimGeneric( t_rec : FuncTableInt, y : RealGeneric, x : RealGeneric ) : RealGeneric = {
//     println("==================================================")
    assert(x.spec == y.spec)

    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val fracW     = t_rec.bp
    val extraBits = fracW - manW

    val xsgn   = x.sgn
    val ysgn   = y.sgn

//     println(f"x = ${x.sgn}|${x.ex}(${xex})|${xman.toLong.toBinaryString}")
//     println(f"y = ${y.sgn}|${y.ex}(${yex})|${yman.toLong.toBinaryString}")
//     println(f"x = ${x.toDouble}, y = ${y.toDouble}, atan2(y, x) = ${atan2(y.toDouble, x.toDouble)}")

    val yIsLarger = slice(0, x.spec.W-1, x.value) < slice(0, x.spec.W-1, y.value)

    val minxy = if (yIsLarger) { x } else { y }
    val maxxy = if (yIsLarger) { y } else { x }


    // ------------------------------------------------------------------------
    // reciprocal table

    val recMan = if (t_rec.nOrder==0) {
      val adr = maxxy.man.toInt
      if (maxxy.man==0) {
        0
      } else {
        t_rec.interval(adr).eval(0L, 0)
      }
    } else {
      val adrW = t_rec.adrW
      val dxbp = manW-adrW-1
      val d    = (SafeLong(1)<<dxbp) - slice(0, manW-adrW, maxxy.man)-1
      val adr  = maskI(adrW)-slice(manW-adrW, adrW, maxxy.man).toInt
      if (maxxy.man==0) {
        0
      } else {
        t_rec.interval(adr).eval(d.toLong, dxbp)
      }
    }

    // ------------------------------------------------------------------------
    // atan2 stage1 postprocess (minxy * rec(maxxy))

    val denomW1 = (1<<fracW) + recMan
    val numerW1 = (1<<manW) + minxy.man

    val zProd = denomW1 * numerW1
    val bp    = fracW + manW

    val zProdMoreThan2 = bit(fracW+1 + manW+1 - 1, zProd)
    val roundBits = bp - manW + zProdMoreThan2
    val zProdRounded = roundBySpec(RoundSpec.roundToEven, roundBits, zProd)
    val zProdMoreThan2AfterRound = bit(manW+1, zProdRounded)
    val zMan = zProdRounded >> zProdMoreThan2AfterRound

    val zSgn = 0
    val zEx0 = minxy.ex - maxxy.ex - 1 + exBias + zProdMoreThan2 + zProdMoreThan2AfterRound
    val zEx = if(zEx0 < 0) {0} else if((1<<exW) <= zEx0) {maskI(exW)} else {zEx0}

    new RealGeneric(x.spec, zSgn, zEx, zMan)
  }
}
