//% @file acosStage1Sim.scala
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

// 1. x -> sqrt(1 - |x|)
// 2. x -> acos(1 - x^2)
//
// acos(1 - sqrt(1-|x|)^2) = acos(1-(1-|x|)) = acos(|x|)
//
object ACosStage1Sim {
  def acosStage1SimGeneric(
    tSqrt: FuncTableInt, x: RealGeneric
  ) : (RealGeneric, Boolean, Int) = { // neg flag, special flag (0: x==0, 1: |x|==1, 2: otherwise)
    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val extraBits = tSqrt.bp - manW
    assert(extraBits == tSqrt.bp - manW)

    val negFlag = x.sgn == 1
    val specialFlag =
      if (x.ex == 0) {0} else if(x.ex == exBias && x.man == 0) { 1 } else { 2 }

//     println(f"x = ${x.sgn}|${x.ex}(${ex})|${man.toLong.toBinaryString}")
//     println(f"x = ${x.toDouble}, acos(x) = ${acos(x.toDouble)}")

    // ------------------------------------------------------------------------
    // special value handling

    if(x.isNaN) {
      return (RealGeneric.nan(spec), negFlag, specialFlag)
    }
    if(x.isInfinite) {
      return (RealGeneric.nan(spec), negFlag, specialFlag)
    }
    if(exBias <= x.ex) { // 1 <= |x|
      return (new RealGeneric(spec, 0.0), negFlag, specialFlag)
    }

    assert(x.ex < exBias)
    assert(-1.0 < x.toDouble && x.toDouble < 1.0)

    // ------------------------------------------------------------------------
    // calc 1 - |x|
    //
    //           +1
    //      manW |+1  manW
    //    .-----.|| .-----.
    //  1.00....000 00....0
    // -)         1.00....0
    // --------------------
    //  0.11....111 00....0
    //  1.11....110 00....0 x2^-1 : normalize
    //  1.00....0 x 2^0           : rounding
    //

    val xShift   = Seq(manW+2, exBias - x.ex).min // > 0 always
    val xSubtracted = (SafeLong(1) << (manW + xShift)) - x.manW1

    val xSubLargerThan1 = xSubtracted > (SafeLong(1) << manW)

    // normalize
    val xNormalizeShift = if(xSubLargerThan1) { binaryWidthSL(xSubtracted) - (manW+1) } else { (manW+1) - binaryWidthSL(xSubtracted) }
    val xNormalized     = if(xSubLargerThan1) { xSubtracted >> xNormalizeShift        } else { xSubtracted << xNormalizeShift        }
    val xNormalizedEx   = if(xSubLargerThan1) { exBias + xNormalizeShift - xShift     } else { exBias - xNormalizeShift - xShift     }
    val xNormalizedMan  = slice(0, manW, xNormalized)

//     println(f"x      = 1.${x.man.toLong.toBinaryString} x 2^${x.ex-exBias}")
//     println(f"xShift = ${xShift}")
//     println(f"     1 = ${(SafeLong(1) << (manW + xShift)).toLong.toBinaryString}%50s")
//     println(f"     x = ${x.manW1.toLong.toBinaryString}%50s")
//     println(f"xSub   = ${xSubtracted.toLong.toBinaryString}%50s")
//     println(f"xNrmSft= ${xNormalizeShift}")
//     println(f"xNrmEx = ${xNormalizedEx-exBias}")
//     println(f"xNrm   = ${xNormalized.toLong.toBinaryString}%50s")
//     println(f"xNrmMan= ${xNormalizedMan.toLong.toBinaryString}%50s")

    // xNormalized never be zero because if 1 <= |x| we already returned. |x|<1
    assert(bit(manW, xNormalized) == 1 && xNormalized < (SafeLong(1) << (manW+1)))

    // ------------------------------------------------------------------------
    // calc sqrt(1-|x|)
    //
    // 1 - |x| is in (0, 1].

    if (xNormalizedEx == exBias) {
      assert(xNormalizedMan == 0)
      return (new RealGeneric(spec, 1.0), negFlag, specialFlag)
    }

    val y = new RealGeneric(x.spec, 0, xNormalizedEx, xNormalizedMan)

//     println(f"sim: xNormalizedEx  = ${xNormalizedEx}")
//     println(f"sim: xNormalizedMan = ${xNormalizedMan.toLong.toBinaryString}")
//     println(f"y = 1 - |x| = ${y.toDouble} should be ${1.0 - abs(x.toDouble)}")

    // use the same calc/postproc phase as sqrt.
    val sqrt = SqrtSim.sqrtSimGeneric(tSqrt, y)

    return (sqrt, negFlag, specialFlag)
  }

  def calcExAdrW(spec: RealSpec): Int = {
    return 1
  }

  def sqrtTableGeneration( order: Int, adrW: Int, manW: Int, fracW: Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
    ) = {
    SqrtSim.sqrtTableGeneration(
      order, adrW, manW, fracW, calcWidthSetting, cbitSetting)
  }
}
