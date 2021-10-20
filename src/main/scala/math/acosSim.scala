//% @file acosSim.scala
//
// Simulators for acos(x) function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

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

object ACosSim {

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
  def acosSimGeneric( ts : Seq[FuncTableIntFixedWidth], x: RealGeneric ) : RealGeneric = {
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

    println(f"x = ${x.sgn}|${x.ex}(${ex})|${man.toLong.toBinaryString}")
    println(f"x = ${x.toDouble}, acos(x) = ${acos(x.toDouble)}")

    if(x.isNaN) {
      return RealGeneric.nan(x.spec)
    }
    if(x.isInfinite) {
      return RealGeneric.nan(x.spec)
    }
    if(0 <= ex) { // 1 <= |x|
      if (ex == 0 && man == 0) { // |x| == 1
        if (sgn == 0) {
          return new RealGeneric(x.spec, 0.0)
        } else {
          return new RealGeneric(x.spec, Pi)
        }
      } else {
        return RealGeneric.nan(x.spec)
      }
    }

    val zSgn = 0 // returns [0, 2pi)

    // here ex < 0.

    val constThreshold  = -manW                   // -23, if FP32
    val linearThreshold = -math.ceil((manW - 1) / 3.0) //  -7
    val halfPi = new RealGeneric(x.spec, Pi * 0.5)

    println(f"constant threshold = ${constThreshold}")
    println(f"linear   threshold = ${linearThreshold}")

    // here the x.sgn has no effect because |x| << 1.
    if (ex < constThreshold) {
      println(f"x < constantThreshold(${pow(2.0, -constThreshold)})")
      return halfPi
    } else if (ex < linearThreshold) { // Pi/2 - x
      println(f"x < linearThreshold(${pow(2.0, -linearThreshold)})")
      // Pi / 2 ~ 1.57... and x < 2^(linearThreshold). we don't need to check MoreThan2.

//       val simz = new RealGeneric(x.spec, Pi * 0.5 - x.toDouble)
//       val refz = new RealGeneric(x.spec, acos(x.toDouble))
//       println(f"pi/2 - x = ${simz.toDouble}(${simz.sgn}|${simz.ex}|${simz.man.toLong.toBinaryString})")
//       println(f"acos(x)  = ${refz.toDouble}(${refz.sgn}|${refz.ex}|${refz.man.toLong.toBinaryString})")

      val shiftOffset = -linearThreshold.toInt
      val calcW       = (manW + 1) + shiftOffset
      val halfPiMan   = ((1 << manW) + halfPi.man) << shiftOffset
      val xManAligned = ((1 << manW) + x.man) >> (linearThreshold.toInt - ex)
      val xManSigned  = if (sgn == 0) { ~(xManAligned.toLong) + 1 } else { xManAligned.toLong }
      val zMan0       = slice(0, calcW, halfPiMan + xManSigned) // rm the leading 1

      println(f"shiftOffset = ${shiftOffset}")
      println(f"calcW       = ${calcW      }")
      println(f"halfPiMan   = ${halfPiMan  .toLong.toBinaryString}")
      println(f"xManAligned = ${xManAligned.toLong.toBinaryString}")
      println(f"xManSigned  = ${xManSigned .toLong.toBinaryString}")
      println(f"zMan0       = ${zMan0      .toLong.toBinaryString}")

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
    } else { // use table.
      println(f"linearThreshold < x. use table.")

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
        val res0  = t.interval(adr).eval(0L, 0)
        // fix range < (1<<bp)
        val resClamp = if(res0 >= (1 << calcW)) {maskL(bp)} else {res0}
        val shift    = calcW+1 - resClamp.toLong.toBinaryString.length
        val res      = (resClamp << shift).toLong - (1 << calcW)

        (-shift, res)

      } else {
        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, man).toInt

        println(f"man = ${man.toLong.toBinaryString}")
        println(f"d   = ${d  .toLong.toBinaryString}")
        println(f"adr = ${adr.toLong.toBinaryString}(${adr})")

        val halfPiFixed = math.round(Pi * 0.5 * (1<<calcW))

        // pi/2 - acos(x)
        val res0 = t.interval(adr).eval(d.toLong, dxbp)
        val res  = if (sgn == 1) {
          halfPiFixed + res0
        } else {
          halfPiFixed - res0
        }
        val shift = calcW+2 - res.toLong.toBinaryString.length
        val resShifted = ((res << shift).toLong) >> 1

        println(f"res0   = ${res0.toLong.toBinaryString}")
        println(f"res0   = ${res0.toDouble / (1<<calcW)}, pi/2-acos(x) = ${Pi*0.5 - acos(x.toDouble)}")
        println(f"halfPi = ${halfPiFixed.toLong.toBinaryString}")
        println(f"res0-hp= ${(-res0 + halfPiFixed).toDouble / (1<<calcW)}, acos(x) = ${acos(x.toDouble)}")
        println(f"res    = ${res .toLong.toBinaryString}")
        println(f"shift  = ${shift}")
        println(f"resS   = ${resShifted.toLong.toBinaryString}")

        (-shift+1, resShifted - (1 << calcW))
      }
      println(f"zman   = ${zman.toLong.toBinaryString}")

      val zmanRound = if (extraBits>0) {roundBySpec(RoundSpec.roundToEven, extraBits, zman)} else {zman}

      println(f"zmanR  = ${zmanRound.toLong.toBinaryString}")

      val z = if (zman<0) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
        0L
      } else if (zmanRound >= (1L<<manW)) {
        println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
        maskL(manW)
      } else {
        zmanRound
      }
      println(f"Sim: zman = ${z.toBinaryString}")

      return new RealGeneric(x.spec, zSgn, zEx + exBias, SafeLong(z))
    }
  }

  def acosTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    val linearThreshold = -math.ceil((manW - 1) / 3.0).toInt

    if (order == 0 || adrW >= manW) {
      val maxCalcWidth = (-2 to linearThreshold by -1).map(i => {
        val tableD = new FuncTableDouble( x => acos(scalb(1.0 + x, i)), 0)
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW ) // convert float table into int
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

      (-1 to linearThreshold by -1).map( i => {
        val tableD = new FuncTableDouble( x => acos(scalb(1.0 + x, i)), 0)
        tableD.addRange(0.0, 1.0, 1<<manW)
        new FuncTableIntFixedWidth( tableD, fracW, maxCalcWidth )
      })
    } else {
      val maxCalcWidth = (-1 to linearThreshold by -1).map(i => {
        val tableD = new FuncTableDouble( x => (Pi * 0.5) - acos(scalb(1.0 + x, i)), order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW ) // convert float table into int
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

      // ex == -1 corresponds to the range [0.5, 1).
      (-1 to linearThreshold by -1).map( i => {
        val tableD = new FuncTableDouble( x => (Pi * 0.5) - acos(scalb(1.0 + x, i)), order )
        tableD.addRange(0.0, 1.0, 1<<adrW)
        new FuncTableIntFixedWidth( tableD, fracW, maxCalcWidth )
      })
    }
  }

  val acosF32TableI = ACosSim.acosTableGeneration( 2, 8, 23, 23+2 )
  val acosF32Sim    = acosSimGeneric(acosF32TableI, _ )

  val acosBF16TableI = ACosSim.acosTableGeneration( 0, 7, 7, 7 )
  val acosBF16Sim    = acosSimGeneric(acosBF16TableI, _ )
}
