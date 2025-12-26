//% @file scaleMixtureGaussianSim.scala
//
// A Simulator of the differentiated log-likelyhood of scaleMixtureGaussian
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._

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

object ScaleMixtureGaussianSim {

  // let
  //   1/sgm'^2 = (sgmA^2 - sgmB^2) / sgmA^2sgmB^2
  //   g(x) = exp(-x^2 / sgm'^2)
  // where
  //   f(x) = -x/sgm'^2 [ 1 / (sgmB/(sgmA*g(x))+1) + sgm'^2/sgmA^2 ]
  //
  //        = -x/sgmA^2 [ (sgmA^2/sgm'^2) / (sgmB/(sgmA*g(x))+1) + 1 ]
  //                      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  //                      this part will be approximated
  //
  def scaleMixtureGaussianSimGeneric(
    t: FuncTableInt, x: RealGeneric, sgmAd: Double = exp(-1.0), sgmBd: Double = exp(-6.0),
    debugDump: Boolean = false
  ): RealGeneric = {

    if(debugDump) {
      println(f"-------------------------------------------------------------")
      println(f"x    = ${x.toDouble}(${x.sgn}|${x.ex}(${x.ex-x.spec.exBias})|${x.man.toLong.toBinaryString})")
      println(f"xman = ${x.man.toLong.toBinaryString}")
    }

    val spec = x.spec

    assert(sgmBd / sgmAd < 1.0, "sgmB / sgmA << 1")

    val adrW   = t.adrW
    val nOrder = t.nOrder
    val bp     = t.bp

    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias
    val extraBits = bp - manW

    val zSgn = if(x.sgn == 0) { 1 } else { 0 }

    if (x.isNaN)      {return RealGeneric.nan (x.spec)}
    if (x.isInfinite) {return RealGeneric.inf (x.spec, zSgn)}
    if (x.isZero)     {return RealGeneric.zero(x.spec)}

    // ------------------------------------------------------------------------
    // scale x to table domain range

    val xTableMaxDigit = tableDomainDigit(manW, sgmAd, sgmBd)
    val xTableShift0   = xTableMaxDigit - (x.ex - exBias)
    val xTableShift    = abs(xTableShift0)
    val xTableMan      = x.manW1 >> xTableShift

    if(debugDump) {
      val xTableMax = pow(2.0, xTableMaxDigit)
      println(f"xtable: bits = ${xTableMan.toLong.toBinaryString}")
      println(f"xtable: ref = ${x.toDouble / xTableMax}, cir = ${xTableMan.toDouble * pow(2.0, -manW)}")
    }

    // ------------------------------------------------------------------------
    // use table (if x < tableDomain)

    val zTable0 = if (nOrder==0) {
      if (adrW<manW) {
        println("WARNING: table address width < mantissa width, " +
          "for polynomial order is zero. address width set to mantissa width.")
      }
      val adr = slice(0, manW, xTableMan.toInt)
      t.interval(adr).eval(0L, 0)
    } else {
      val dxbp = manW-adrW-1
      val d    = slice(0, manW-adrW, xTableMan) - (SafeLong(1)<<dxbp)
      val adr  = slice(manW-adrW, adrW, xTableMan).toInt
      t.interval(adr).eval(d.toLong, dxbp).toLong
    }

    // if x is large enough, then the first term is zero.
    val zTable = if(xTableShift0 > 0) { zTable0 } else { 0L }

    // ------------------------------------------------------------------------
    // add 1 to the table result

    if(debugDump) {
      val sgmA2d = sgmAd * sgmAd
      val sgmB2d = sgmBd * sgmBd
      val sgmPrime2d = (sgmA2d * sgmB2d) / (sgmA2d - sgmB2d)

      val maxBitDigit = tableMaxBitDigit(sgmAd, sgmBd)
      val zTableD = zTable * pow(2.0, maxBitDigit-t.bp)
      val xd = x.toDouble

      val g     = exp(-xd*xd / (2*sgmPrime2d))
      val numer = sgmA2d / sgmPrime2d
      val denom = (sgmBd / (sgmAd * g)) + 1.0

      println(f"ztable:                 2         1         ")
      println(f"ztable:           65432109876543210987654321")
      println(f"ztable: zTable0 = ${zTable0.toBinaryString}")
      println(f"ztable: zTable  = ${zTable.toBinaryString}")
      println(f"ztable: ref = ${numer / denom}, cir = ${zTableD}")
    }

    val zTableScaleDigit = tableMaxBitDigit(sgmAd, sgmBd)
    val zTableScaled     = SafeLong(zTable) << zTableScaleDigit

    val zman0        = zTableScaled + (SafeLong(1) << t.bp)
    val zman0W       = binaryWidthSL(zman0)
    val zman0Shift   = zman0W - (manW+1+extraBits)
    val zman0Rounded = (zman0 >> zman0Shift) + bit(zman0Shift-1, zman0)
    val zman0MoreThan2AfterRound = bit(manW+1+extraBits+1, zman0Rounded)

    // z before multiplication
    val z0ManW1 = if(zman0MoreThan2AfterRound == 1) {
      zman0Rounded >> 1
    } else {zman0Rounded}

    val z0ExNoBias = zman0W-1 - t.bp + zman0MoreThan2AfterRound

    if(debugDump) {
      println(f"z0: zTableScaled = ${zTableScaled.toLong.toBinaryString}")
      println(f"z0: 1            = ${(SafeLong(1) << t.bp).toLong.toBinaryString}")
      println(f"z0: zman0      = ${zman0.toLong.toBinaryString}")
      println(f"z0: zman0W     = ${zman0W}")
      println(f"z0: zman0Shift = ${zman0Shift}")

      val z0 = new RealGeneric(spec, zSgn, z0ExNoBias + exBias, slice(extraBits, manW, z0ManW1))

      val sgmA2d = sgmAd * sgmAd
      val sgmB2d = sgmBd * sgmBd
      val sgmPrime2d  = (sgmA2d * sgmB2d) / (sgmA2d - sgmB2d)
      val sgmA2overP2 = sgmA2d / sgmPrime2d

      val xd    = x.toDouble
      val g     = exp(-xd*xd / (2*sgmPrime2d))
      val term1 = (sgmA2d / sgmPrime2d) / ((sgmBd / (sgmAd * g)) + 1.0)
      val term2 = 1.0

      val zr = new RealGeneric(spec, term1 + term2)

      println(f"z0: ref = ${-(term1 + term2)}, cir = ${z0.toDouble}")
      println(f"z0: ref.ex = ${zr.ex - exBias}(${zr.ex}) ref.man = ${zr.manW1.toLong.toBinaryString}")
      println(f"z0: cir.ex = ${z0.ex - exBias}(${z0.ex}) cir.man = ${z0.manW1.toLong.toBinaryString}")
    }

    // ------------------------------------------------------------------------
    // multiply x/sgmA2

    // x / sgmA2
    val rsgmA2 = new RealGeneric(spec, 1.0 / (sgmAd * sgmAd))
    val xOverSgmA2Prod = x.manW1 * rsgmA2.manW1
    val xOverSgmA2ProdMoreThan2 = bit((1+manW) + (1+manW) - 1, xOverSgmA2Prod)
    val xOverSgmA2ProdRoundedW1 = Rounding.roundToEven(
      manW + xOverSgmA2ProdMoreThan2, xOverSgmA2Prod)
    val xOverSgmA2ProdMoreThan2AfterRounding = bit(manW+1, xOverSgmA2ProdRoundedW1)
    val xOverSgmA2ProdEx = x.ex + rsgmA2.ex - rsgmA2.spec.exBias +
      xOverSgmA2ProdMoreThan2 + xOverSgmA2ProdMoreThan2AfterRounding

    if(debugDump) {
      val sgmA2d = sgmAd * sgmAd
      val xsgmA2 = new RealGeneric(spec, x.sgn, xOverSgmA2ProdEx, slice(0, manW, xOverSgmA2ProdRoundedW1))
      println(f"x/sgmA2: ref = ${x.toDouble / sgmA2d}, cir = ${xsgmA2.toDouble}")

      println(f"x/sgmA2: ex = ${xsgmA2.ex}, man = ${xsgmA2.manW1.toLong.toBinaryString}")
    }

    // (x / sgmA2) * z0

    val zProd = xOverSgmA2ProdRoundedW1 * z0ManW1
    val zProdMoreThan2 = bit((1+manW)+(1+manW+extraBits)-1, zProd)
    val zProdRoundedW1 = Rounding.roundToEven(
      (manW + extraBits) + zProdMoreThan2, zProd)
    val zProdMoreThan2AfterRound = bit(manW+1, zProdRoundedW1)
    val zProdEx = xOverSgmA2ProdEx + z0ExNoBias +
      zProdMoreThan2 + zProdMoreThan2AfterRound

    if(debugDump) {
      println(f"lhs = ${z0ManW1.toLong.toBinaryString}")
      println(f"rhs = ${xOverSgmA2ProdRoundedW1.toLong.toBinaryString}")
      println(f"zProd = ${zProd}")
      println(f"zProdRoundedW1 = ${zProdRoundedW1.toLong.toBinaryString}")
      println(f"xOverSgmA2ProdEx = ${xOverSgmA2ProdEx}, zProdEx = ${zProdEx}")
    }

    val zProdMan = slice(0, manW, zProdRoundedW1)

    if(zProdEx >= (1 << spec.exW) - 1) {
      return new RealGeneric(spec, zSgn, (1 << spec.exW) - 1, 0)
    } else {
      return new RealGeneric(spec, zSgn, zProdEx, zProdMan)
    }
  }

  def tableDomainDigit(manW: Int, sgmA: Double, sgmB: Double): Int = {
    val log2 = (x: Double) => {log(x) / log(2.0)}
    // Domain starts from 0 to the point where the first term is enough smaller
    // compared to the second term.
    //
    // f(x) = -x/sgm'^2 [ 1 / (sgmB/(sgmA*g)+1) + sgm'^2/sgmA^2 ]
    //                    ^^^^^^^^^^^^^^^^^^^^^
    //                     -> 0 when g -> 0
    //
    // The preprocessor passes x in this range [0, domainMax) converted to [0, 1).
    // the table should consider this conversion while calculating the value.

    val sgmA2 = sgmA * sgmA
    val sgmB2 = sgmB * sgmB
    val sgmPrime2 = (sgmA2 * sgmB2) / (sgmA2 - sgmB2)

    val boundary = sqrt(2.0 * sgmPrime2 * log(
      sgmA / sgmB * (pow(2.0, manW) * sgmA2 / sgmPrime2 - 1.0)
    ))

    scala.math.ceil(log2(boundary)).toInt
  }

  def tableMaxValue(sgmA: Double, sgmB: Double): Double = {
    // table takes the max value at x = 0, because g(0) = 1 and g(inf) -> 0.
    //
    // f(x) = -x/sgm'^2 [ 1 / (sgmB/(sgmA*g)+1) + sgm'^2/sgmA^2 ]
    //      = -x/sgmA^2 [ sgmA^2/sgm'^2 / (sgmB/(sgmA*g)+1)  + 1]
    //
    // P(x) approximates [sgmA^2/sgm'^2 / (sgmB/(sgmA*g)+1)] / max(P).
    // P(inf) = 0
    // P(0) = 1 / (sgmB/sgmA+1) * sgmA^2/sgm'^2
    //      = 1 / (sgmB/sgmA+1) * sgmA^2(sgmA^2 - sgmB^2)/sgmA^2sgmB^2
    //      = 1 / (sgmB/sgmA+1) * (sgmA^2 - sgmB^2)/sgmB^2
    //      = sgmA / (sgmB+sgmA) * (sgmA^2 - sgmB^2)/sgmB^2
    //      = sgmA / (sgmB+sgmA) * (sgmA + sgmB)(sgmA - sgmB)/sgmB^2
    //      = sgmA (sgmA - sgmB)/sgmB^2
    //      = sgmA / sgmB * [(sgmA - sgmB) / sgmB]
    //
    (sgmA / sgmB) * ((sgmA - sgmB) / sgmB)
  }

  def tableMaxBitDigit(sgmA: Double, sgmB: Double): Int = {
    val log2 = (x: Double) => {log(x) / log(2.0)}

    val maxDigit = scala.math.ceil(log2(tableMaxValue(sgmA, sgmB))).toInt
    assert(maxDigit > 0)
    maxDigit
  }

  def tableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int,
    sgmA: Double, sgmB: Double,
    calcWidthSetting: Option[Seq[Int]] = None,
    cbitSetting: Option[Seq[Int]] = None
  ) = {

    val maxDigit = tableMaxBitDigit(sgmA, sgmB)
    // f(x) = 1 / (sgmB/(sgmA*g)+1) * sgmA^2/sgm'^2, g = exp(-x^2/sgm'^2)

    val sgmBoverA = sgmB / sgmA
    val sgmA2 = sgmA * sgmA
    val sgmB2 = sgmB * sgmB
    val sgmPrime2 = (sgmA2 * sgmB2) / (sgmA2 - sgmB2)
    val sgmA2overP2 = sgmA2 / sgmPrime2

    val f = (x: Double) => {
      val g = exp(-x * x / (2*sgmPrime2))
      sgmA2overP2 / (sgmBoverA / g + 1.0)
    }

    val domainMax = pow(2.0, tableDomainDigit(manW, sgmA, sgmB))

    val tableD = new FuncTableDouble( (x0: Double) => {
      val x = domainMax * x0 // x0 in [0, 1) -> x in [0, domainMax)
      f(x) * pow(2.0, -maxDigit) // scales it into [0, 1)
    }, order )

    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW, calcWidthSetting, cbitSetting )
  }

  val smg16F32TableI = ScaleMixtureGaussianSim.tableGeneration( 2, 8, 23, 23+2, exp(-1.0), exp(-6.0) )
  val smg16F32Sim = scaleMixtureGaussianSimGeneric(smg16F32TableI, _, exp(-1.0), exp(-6.0), false )

  val smg16BF16TableI = ScaleMixtureGaussianSim.tableGeneration( 0, 7, 7, 7+2, exp(-1.0), exp(-6.0) )
  val smg16BF16Sim = scaleMixtureGaussianSimGeneric(smg16BF16TableI, _, exp(-1.0), exp(-6.0), false )
}
