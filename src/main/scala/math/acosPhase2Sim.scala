//% @file acosPhase2Sim.scala
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

// 1. x -> sqrt(1 - |x|) = y
// 2. y -> acos(1 - y^2)
//
// acos(1 - sqrt(1-|x|)^2) = acos(1-(1-|x|)) = acos(|x|)
//
private[rial] object ACosPhase2Sim {
  def acosPhase2SimGeneric(
    tACos: FuncTableInt, x: RealGeneric
  ): RealGeneric = {
    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val adrW   = tACos.adrW
    val extraBits = tACos.bp - manW
    val calcW = manW + extraBits

    // ------------------------------------------------------------------------
    // special value handling

    if(x.ex == maskI(spec.exW)) { //
      if(x.man == 1) { // x == 0, y == pi/2
        return new RealGeneric(x.spec, Pi * 0.5)
      } else if (x.man == 2) { // x >= +1, y == +0
        return new RealGeneric(x.spec, 0.0)
      } else if (x.man == 3) { // x <= -1, y == pi
        return new RealGeneric(x.spec, Pi)
      } else { // normal nan.
        return RealGeneric.nan(spec)
      }
    }

//     println(f"x = ${x.sgn}|${x.ex}(${ex})|${man.toLong.toBinaryString}")
//     println(f"x = ${x.toDouble}, acos(x) = ${acos(x.toDouble)}")

    assert(0.0 <= abs(x.toDouble) && abs(x.toDouble) <= 1.0, f"0.0 <= abs(x(${x.toDouble})) <= 1.0")

    // ------------------------------------------------------------------------
    // calc `acos(1-x^2)/x - 1` by table. then multiply x at the postprocess.
    // - acos(1-x^2)/x --(x->+0)--> sqrt(2) ~ 1.414.
    // - acos(1-x^2)/x --(x-> 1)--> pi/2 ~ 1.57.

    val xAligned = if(x.ex == 0) { SafeLong(0) } else { x.manW1 >> (exBias - x.ex) }

    val dxbp = manW-adrW-1
    val d    = slice(0, dxbp+1, xAligned) - (SafeLong(1) << dxbp)
    val adr  = slice(dxbp+1, adrW, xAligned)

    val polynomial = tACos.interval(adr.toInt).eval(d.toLong, dxbp).toSafeLong

    // ------------------------------------------------------------------------
    // multiply zfrac ( = acos(1-x^2)/x ) and x to get acos(1-x^2).

    val lhsManW1 = polynomial + (SafeLong(1) << calcW) // acos(1-x^2)/x
    val rhsManW1 = x.manW1

    val zProd = lhsManW1 * rhsManW1
    val zMoreThan2 = bit((calcW+1)+(manW+1)-1, zProd)
    val zRoundBit = calcW + zMoreThan2
    val zRounded = roundBySpec(RoundSpec.roundToEven, zRoundBit, zProd)
    val zMoreThan2AfterRound = bit(manW+1, zRounded)
    val zManW1 = zRounded >> zMoreThan2AfterRound

    val lhsEx = exBias
    val rhsEx = x.ex
    val zEx = lhsEx + rhsEx - exBias + zMoreThan2 + zMoreThan2AfterRound

    if (x.sgn == 0) {
      return new RealGeneric(x.spec, 0, zEx, slice(0, manW, zManW1))
    }

    // ------------------------------------------------------------------------
    // subtract z from pi if x was negative

    val pi = new RealGeneric(spec, Pi)
    val zShift   = Seq(1+manW+2, exBias - zEx).min
    assert(zShift >= 0)
    val piAligned = (pi.manW1 << 3)
    val zAligned  = (zManW1 << 2) >> zShift
    val negZ  = piAligned - zAligned

    val negZNormalizeShift = (1+manW+3) - binaryWidthSL(negZ)
    val negZNormalized     = (negZ << negZNormalizeShift)
    assert(bit(manW+3, negZNormalized) == 1)

    val negZRounded          = slice(3, manW, negZNormalized) + bit(2, negZNormalized)
    val negZRoundedMoreThan2 = bit(manW, negZRounded)
    val negZNormalizedEx   = (exBias+1 - negZNormalizeShift) + negZRoundedMoreThan2
    val negZNormalizedMan  = slice(0, manW, negZRounded)

    return new RealGeneric(x.spec, 0, negZNormalizedEx, negZNormalizedMan)
  }

  def acosTableGeneration( order: Int, adrW: Int, manW: Int, fracW: Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
  ) = {
    val f = (x: Double) => {
      if (x == 0.0) {
        // acos(1 - |x|^2) / x -> sqrt(2) when x -> 0
        sqrt(2.0) - 1.0
      } else {
        acos(1.0 - x * x) / x - 1.0
      }
    }
    val tableD = new FuncTableDouble( f, order )

    tableD.addRange(0.0, 1.0, 1<<adrW) // this makes resulting table adrW+1
    new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }
}
