//% @file acosStage2Sim.scala
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
object ACosStage2Sim {
  def acosStage2SimGeneric(
    tACos: FuncTableInt, x: RealGeneric, xneg: Boolean, special: Int
  ): RealGeneric = { // special flag (0: x==0, 1: |x|==1, 2: otherwise)
    val spec   = x.spec
    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val adrW   = tACos.adrW
    val extraBits = tACos.bp - manW
    val calcW = manW + extraBits

    // ------------------------------------------------------------------------
    // special value handling

    if(x.isNaN) {
      return RealGeneric.nan(spec)
    }
    if(x.isInfinite) {
      return RealGeneric.nan(spec)
    }

    if (special == 0) { // if x == 0 {return pi/2}
      return new RealGeneric(x.spec, Pi * 0.5)
    } else if (special == 1) {
      if (xneg) { // acos(-1)
        return new RealGeneric(x.spec, Pi)
      } else {
        return new RealGeneric(x.spec, 0.0)
      }
    }

//     println(f"x = ${x.sgn}|${x.ex}(${ex})|${man.toLong.toBinaryString}")
//     println(f"x = ${x.toDouble}, acos(x) = ${acos(x.toDouble)}")

    val zSgn = 0 // returns acos returns z in [0, pi], so always z >= 0.

    assert(0.0 <= x.toDouble && x.toDouble <= 1.0, f"0.0 <= x(${x.toDouble}) <= 1.0")

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
    val zRounded = (zProd >> (calcW+zMoreThan2)) + bit(calcW+zMoreThan2-1, zProd)
    val zMoreThan2AfterRound = bit(manW+1, zRounded)
    val zManW1 = if(zMoreThan2AfterRound == 1) {SafeLong(1) << manW} else {zRounded}

    val lhsEx = exBias
    val rhsEx = x.ex
    val zEx = lhsEx + rhsEx - exBias + zMoreThan2 + zMoreThan2AfterRound

    // ------------------------------------------------------------------------
    // subtract from pi if x was negative

    // TODO

    return new RealGeneric(x.spec, zSgn, zEx, slice(0, manW, zManW1))
  }

  def calcExAdrW(spec: RealSpec): Int = {
    return 1
  }

  def acosTableGeneration( order: Int, adrW: Int, manW: Int, fracW: Int,
      calcWidthSetting: Option[Seq[Int]] = None,
      cbitSetting: Option[Seq[Int]] = None
  ) = {
    val f = (x: Double) => {
      if (x == 0.0) {
        // acos(1 - |x|^2) / x -> sqrt(2) when x -> 0
        sqrt(2.0)
      } else {
        acos(1.0 - x * x) / x - 1.0
      }
    }
    val tableD = new FuncTableDouble( f, order )

    tableD.addRange(0.0, 1.0, 1<<adrW) // this makes resulting table adrW+1
    new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }
}
