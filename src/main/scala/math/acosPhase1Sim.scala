//% @file acosPhase1Sim.scala
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
private[rial] object ACosPhase1Sim {
  def acosPhase1SimGeneric(
    tSqrt: FuncTableInt, x: RealGeneric
  ) : RealGeneric = {

    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val extraBits = tSqrt.bp - manW
    assert(extraBits == tSqrt.bp - manW)

    val zsgn = x.sgn

//     println(f"x = ${x.sgn}|${x.ex}(${x.ex - exBias})|${x.man.toLong.toBinaryString}")
//     println(f"x = ${x.toDouble}, sqrt(1-|x|) = ${scala.math.sqrt(1.0 - scala.math.abs(x.toDouble)) }")

    // ------------------------------------------------------------------------
    // special value handling

    // special values (represented by NaN boxing).
    // if the MSB of the mantissa is 1, it is NaN. if 0, it is special values.
    // since it uses NaN boxing, flag=0 is used by inf.
    // - 1: x ==  0, y == pi/2
    // - 2: x >=  1, y == +0
    // - 3: x <= -1, y == pi
    //
    if(x.isNaN) { // acos(nan) = nan
      return new RealGeneric(spec, zsgn, maskI(spec.exW), (1 << (spec.manW-1)))
    }
    if(x.isInfinite) { // acos(inf) = nan
      return new RealGeneric(spec, zsgn, maskI(spec.exW), (1 << (spec.manW-1)))
    }
    if(x.ex == 0) { // acos(0) = pi/2
      return new RealGeneric(spec, zsgn, maskI(spec.exW), 1)
    }
    if(exBias <= x.ex) { // 1 <= |x|. round to x == +/- 1. retval will be 0 or pi.
      if(x.sgn == 0) {
        return new RealGeneric(spec, zsgn, maskI(spec.exW), 2) // acos(1) = 0
      } else {
        return new RealGeneric(spec, zsgn, maskI(spec.exW), 3) // acos(-1) = pi
      }
    }

    assert(x.ex < exBias)
    assert(-1.0 < x.toDouble && x.toDouble < 1.0)

    // ------------------------------------------------------------------------
    // calc 1 - |x|
    //

    val xShift = exBias - x.ex

    val xShifted = if (xShift > 1+manW+2) {SafeLong(0)} else {(x.manW1 << 2) >> xShift}
    val xSubtracted = (SafeLong(1) << (manW+2)) - xShifted

    val xNormalizeShift = (1+manW+2) - binaryWidthSL(xSubtracted)
    val xNormalized = (xSubtracted << xNormalizeShift)
    assert(bit(manW+2, xNormalized) == 1 && xNormalized < (SafeLong(1) << (manW+3)))

    val xRounded = (xNormalized >> 2) + bit(1, xNormalized)
    val xRoundedMoreThan2 = bit(manW+1, xRounded)
    val xConvertedEx = exBias - xNormalizeShift + xRoundedMoreThan2
    val xConvertedMan = slice(0, manW, xRounded)

//     println(f"x      = 1.${x.man.toLong.toBinaryString} x 2^${x.ex-exBias}")
//     println(f"xShift = ${xShift}")
//     println(f"     1 = ${(SafeLong(1) << (manW + xShift)).toLong.toBinaryString}%50s")
//     println(f"     x = ${x.manW1.toLong.toBinaryString}%50s")
//     println(f"xSub   = ${xSubtracted.toLong.toBinaryString}%50s")
//     println(f"xNrmSft= ${xNormalizeShift}")
//     println(f"xNrmEx = ${xConvertedEx-exBias}")
//     println(f"xNrm   = ${xNormalized.toLong.toBinaryString}%50s")
//     println(f"xNrmMan= ${xConvertedMan.toLong.toBinaryString}%50s")

    // xNormalized never be zero because if 1 <= |x| we already returned.

    // ------------------------------------------------------------------------
    // calc sqrt(1-|x|)
    //
    // 1 - |x| is in (0, 1].

    if (xConvertedEx == exBias) {
      // 1 - |x| == 1. |x| ~ 0
      assert(xConvertedMan == 0, f"xConvertedMan = ${xConvertedMan.toLong.toBinaryString}")

      return new RealGeneric(spec, zsgn, maskI(spec.exW), 1) // acos(0) = pi/2
    }

    val y = new RealGeneric(x.spec, 0, xConvertedEx, xConvertedMan)

//     println(f"sim: xConvertedEx  = ${xConvertedEx}")
//     println(f"sim: xConvertedMan = ${xConvertedMan.toLong.toBinaryString}")
//     println(f"y = 1 - |x| = ${y.toDouble} should be ${1.0 - abs(x.toDouble)}")

    // use the same calc/postproc phase as sqrt.
    val sqrt = SqrtSim.sqrtSimGeneric(tSqrt, y)

    return new RealGeneric(spec, zsgn, sqrt.ex, sqrt.man)
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
