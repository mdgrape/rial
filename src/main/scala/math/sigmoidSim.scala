//% @file sigmoidSim.scala
//
// Simulators for sigmoid function
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

object SigmoidSim {

  def sigmoidSimGeneric(
    t : FuncTableInt, x: RealGeneric
  ) : RealGeneric = {

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val adrW      = t.adrW
    val fracW     = t.bp
    val extraBits = fracW - manW
    val nOrder    = t.nOrder

//     println("==============================================================")
//     println(f"SigmoidSim: x    = ${x.toDouble}(${x.sgn}|${x.ex}(${x.ex-exBias})|${x.man.toLong.toBinaryString})")

    if (x.isNaN)                    { return RealGeneric.nan (spec)     }
    if (x.isInfinite && x.sgn == 0) { return new RealGeneric(spec, 1.0) }
    if (x.isInfinite && x.sgn == 1) { return RealGeneric.zero(spec)     }

    val zSgn = 0

    if (adrW >= manW && nOrder != 0) {
      println("WARNING: table address width >= mantissa width, but polynomial" +
        " order is not zero. Polynomial order is overwritten by zero.")
    }

    val order = if (adrW>=manW) {0} else {nOrder}
    val calcW = manW + extraBits

    val rangeMaxLog2 = tableRangeMaxLog2(manW)
    val rangeMaxEx = exBias + rangeMaxLog2 - 1
    if(rangeMaxEx < x.ex) {
      if(x.sgn == 0) {
        return new RealGeneric(spec, 1.0)
      } else {
        return RealGeneric.zero(spec)
      }
    }
    assert(x.ex <= rangeMaxEx)

//     val xPadded = (x.manW1 << rangeMaxLog2)
    val xScaled = x.manW1 >> (rangeMaxEx - x.ex)

    val zScaled0 = if (order==0) {
      if (adrW<manW) {
        println("WARNING: table address width < mantissa width, for polynomial order is zero. address width set to mantissa width.")
      }
      val adr = slice(rangeMaxLog2, adrW, xScaled).toInt
      t.interval(adr).eval(0L, 0)
    } else {
      val dxbp = (manW+1)-adrW-1
      val d    = slice(0, (manW+1)-adrW, xScaled) - (SafeLong(1)<<dxbp)
      val adr  = slice((manW+1)-adrW, adrW, xScaled).toInt
      t.interval(adr).eval(d.toLong, dxbp).toLong
    }

    val zScaled = if(zScaled0 < 0) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial value negative at x=$x%h")
      SafeLong(0)
    } else if (zScaled0 >= (SafeLong(1)<<calcW)) {
      println(f"WARNING (${this.getClass.getName}) : Polynomial range overflow at x=$x%h")
      maskSL(calcW)
    } else {
      SafeLong(zScaled0)
    }

//     println(f"zScaled0 = ${zScaled0.toDouble / (1 << calcW).toDouble}")
//     println(f"zScaled  = ${zScaled .toDouble / (1 << calcW).toDouble}")

    // zScaled = 2 * [1 - 1 / (1 + exp(-x)]
    // if x > 0,  return 1 - zScaled / 2
    // otherwise, return     zScaled / 2

    val zShifted = if(x.sgn == 0) { (SafeLong(1) << (calcW+1)) - zScaled } else { zScaled }
    val zShift   = (calcW+1) - binaryWidthSL(zShifted)
    val zmanW1   = zShifted << zShift
    val zmanW1Rounded = (zmanW1 >> extraBits) + bit(extraBits-1, zmanW1)

//     println(f"zShifted = ${zShifted.toDouble / (1 << calcW).toDouble}")
//     println(f"zShift   = ${zShift}")
//     println(f"zmanW1   = ${zmanW1.toDouble / (1 << calcW).toDouble}")
//     println(f"zmanW1R  = ${zmanW1Rounded.toDouble / (1 << manW).toDouble}")

    val zMan = if(zmanW1Rounded >= (SafeLong(1)<<(manW+1))) {
      println(f"WARNING (${this.getClass.getName}) : rounding overflow at x=$x%h")
      maskSL(manW)
    } else {
      slice(0, manW, zmanW1Rounded)
    }

    val zEx = exBias - (zShift+1)

    val z = new RealGeneric(x.spec, zSgn, zEx, zMan)
//     println(f"SigmoidSim: z    = ${z.toDouble}(${z.sgn}|${z.ex}(${z.ex-exBias})|${z.man.toLong.toBinaryString})")
    z
  }

  def tableRangeMaxLog2( manW: Int ): Int = {
    val log2 = (x: Double) => { log(x)/log(2.0) }

    val delta = pow(2.0, -manW)
//     2.0 - 2.0 / (1.0 + exp(-x)) < delta
//     2.0 - delta < 2.0 / (1.0 + exp(-x))
//     1.0 + exp(-x) < 2.0 / (2.0 - delta)
//     -x < log(2.0 / (2.0 - delta) - 1.0)
    val threshold = -1.0 * log(2.0 / (2.0 - delta) - 1.0)

    ceil(log2(threshold)).toInt
  }

  def tableGeneration(
    order : Int, adrW : Int, manW : Int, fracW : Int,
    calcWidthSetting: Option[Seq[Int]] = None,
    cbitSetting: Option[Seq[Int]] = None
  ) = {
    val rangeMaxLog2 = tableRangeMaxLog2(manW)
    val rangeMax     = 1 << rangeMaxLog2
    val f = ( x01: Double ) => {
      val x = x01 * rangeMax // [0, 1) => [0, rangeMax)
      val z = 2.0 - (2.0 / (1.0 + exp(-x)))
      assert(z <= 1.0, f"x = ${x}, z = ${z}")
      z
    }
    val tableD = new FuncTableDouble( f, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }

  val sigmoidF32TableI = SigmoidSim.tableGeneration( 2, 8, 23, 23+2 )
  val sigmoidF32Sim    = sigmoidSimGeneric(sigmoidF32TableI, _ )

  val sigmoidBF16TableI = SigmoidSim.tableGeneration( 0, 7, 7, 7 )
  val sigmoidBF16Sim    = sigmoidSimGeneric(sigmoidBF16TableI, _ )
}
