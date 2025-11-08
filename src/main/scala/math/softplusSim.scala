//% @file softplusSim.scala
//
// Simulators for softplus function
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

private[rial] object SoftPlusSim {

  def softplusSimGeneric(
    ts : Seq[FuncTableInt], x: RealGeneric
  ) : RealGeneric = {

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val tpos = ts(0)
    val tneg = ts(1)

    val adrW      = tpos.adrW
    val fracW     = tpos.bp
    val extraBits = fracW - manW
    val nOrder    = tpos.nOrder
    assert(tneg.adrW   == tpos.adrW  )
    assert(tneg.bp     == tpos.bp    )
    assert(tneg.nOrder == tpos.nOrder)

//     println("==============================================================")
//     println(f"SoftPlusSim: x    = ${x.toDouble}(${x.sgn}|${x.ex}(${x.ex-exBias})|${x.man.toLong.toBinaryString})")

    if (x.isNaN)                    { return RealGeneric.nan (spec)    }
    if (x.isInfinite && x.sgn == 0) { return RealGeneric.inf (spec, 0) }
    if (x.isInfinite && x.sgn == 1) { return RealGeneric.zero(spec)    }

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
        return x
      } else {
        return RealGeneric.zero(spec)
      }
    }
    assert(x.ex <= rangeMaxEx)

    // +1 for hidden bit
    val xScaled = x.manW1 >> (rangeMaxEx - x.ex + 1)

    val zScaled0 = if (order==0) {
      if (adrW<manW) {
        println("WARNING: table address width < mantissa width, for polynomial order is zero. address width set to mantissa width.")
      }
      val adr = xScaled.toInt
      if(x.sgn == 0) {
        tpos.interval(adr).eval(0L, 0)
      } else {
        tneg.interval(adr).eval(0L, 0)
      }
    } else {
      val dxbp = manW-adrW-1
      val d    = slice(0, manW-adrW, xScaled) - (SafeLong(1)<<dxbp)
      val adr  = slice(manW-adrW, adrW, xScaled).toInt

//       println(f"SoftPlusSim: dx  = ${d.toLong.toBinaryString}")
//       println(f"SoftPlusSim: adr = ${adr.toLong.toBinaryString}")

      if(x.sgn == 0) {
        tpos.interval(adr).eval(d.toLong, dxbp).toLong
      } else {
        tneg.interval(adr).eval(d.toLong, dxbp).toLong
      }
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

//     println(f"SoftPlusSim: zScaled0 = ${zScaled0.toDouble / (1 << calcW).toDouble}")
//     println(f"SoftPlusSim: zScaled  = ${zScaled .toDouble / (1 << calcW).toDouble}")
//     println(f"SoftPlusSim: zres     = ${zScaled .toLong.toBinaryString}")

    val zShift = (calcW+1) - binaryWidthSL(zScaled)
    val zmanW1 = zScaled << zShift
    val zEx = if(x.sgn == 0) {
      exBias - zShift + (rangeMaxLog2 + 1)
    } else {
      exBias - zShift
    }
    val zmanW1Rounded = (zmanW1 >> extraBits) + bit(extraBits-1, zmanW1)

//     println(f"SoftPlusSim: zShift   = ${zShift}")
//     println(f"SoftPlusSim: zmanW1   = ${zmanW1.toLong.toBinaryString}(${zmanW1.toDouble / (1 << calcW).toDouble})")
//     println(f"SoftPlusSim: zmanW1R  = ${zmanW1Rounded.toLong.toBinaryString}(${zmanW1Rounded.toDouble / (1 << manW).toDouble})")

    val zMan = if(zmanW1Rounded >= (SafeLong(1)<<(manW+1))) {
      println(f"WARNING (${this.getClass.getName}) : rounding overflow at x=$x%h")
      maskSL(manW)
    } else {
      slice(0, manW, zmanW1Rounded)
    }

    val z = new RealGeneric(x.spec, zSgn, zEx, zMan)
//     println(f"SoftPlusSim: z    = ${z.toDouble}(${z.sgn}|${z.ex}(${z.ex-exBias})|${z.man.toLong.toBinaryString})")
    z
  }

  def tableRangeMaxLog2( manW: Int ): Int = {
    val log2 = (x: Double) => { log(x)/log(2.0) }

    val delta = pow(2.0, -manW)
    // log(1 + exp(x)) - x   < delta
    // log(1 + exp(x))       < delta + x
    //     1 + exp(x)        < exp(delta + x)
    //     1                 < exp(delta + x) - exp(x)
    //     1                 < (exp(delta) - 1) * exp(x)
    // 1 / (exp(delta)-1)    < exp(x)
    // log(1/(exp(delta)-1)) < x

    // the threshold for the negative part results in the same formula.
    //
    // log(1 + exp(x))    < delta
    //     1 + exp(x)     < exp(delta)
    //     1 - exp(delta) < -exp(x)
    //    -1 + exp(delta) > exp(x)
    //  log(exp(delta)-1) > x
    // -log(exp(delta)-1) < -x

    val threshold = log(1.0 / (exp(delta) - 1))

    ceil(log2(threshold)).toInt
  }

  def tableGeneration(
    order : Int, adrW : Int, manW : Int, fracW : Int,
    calcWidthSetting: Option[Seq[Int]] = None,
    cbitSetting: Option[Seq[Int]] = None
  ) = {
    val rangeMaxLog2 = tableRangeMaxLog2(manW)
    val rangeMax     = 1 << rangeMaxLog2

    val fpos = ( x01: Double ) => {
      val x = x01 * rangeMax // [0, 1) => [0, rangeMax)
      // XXX: not so bit-efficient
      val z = log(1.0 + exp(x)) / (rangeMax*2).toDouble
      assert(z <= 1.0, f"fpos: x = ${x}, z = ${z}")
      z
    }
    val tableDpos = new FuncTableDouble( fpos, order )
    tableDpos.addRange(0.0, 1.0, 1<<adrW)

    val fneg = ( x01: Double ) => {
      val x = -1 * x01 * rangeMax // [0, 1) => [0, rangeMax)
      val z = log(1.0 + exp(x)) // 0 ~ ln(2)
      assert(z <= 1.0, f"fneg: x = ${x}, z = ${z}")
      z
    }
    val tableDneg = new FuncTableDouble( fneg, order )
    tableDneg.addRange(0.0, 1.0, 1<<adrW)

    Seq(
      new FuncTableInt( tableDpos, fracW, calcWidthSetting, cbitSetting ),
      new FuncTableInt( tableDneg, fracW, calcWidthSetting, cbitSetting ),
    )
  }

  val softplusF32TableI = SoftPlusSim.tableGeneration( 2, 8, 23, 23+2 )
  val softplusF32Sim    = softplusSimGeneric(softplusF32TableI, _ )

  val softplusBF16TableI = SoftPlusSim.tableGeneration( 0, 7, 7, 7 )
  val softplusBF16Sim    = softplusSimGeneric(softplusBF16TableI, _ )
}
