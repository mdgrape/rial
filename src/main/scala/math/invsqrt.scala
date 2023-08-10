//% @file invsqrt.scala
//
// 1 / square root function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith.RealSpec
import rial.arith.FloatChiselUtil

import rial.math._

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

// the same as sqrt.

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class InvSqrtTableCoeff(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
) extends Module {

  val manW   = spec.manW
  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input  (UInt((1+adrW).W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {
    val tbl = VecInit( (0L to 1L<<(adrW+1)).map(
      n => {
        val x = if (n < (1L<<adrW)) {
          (n.toDouble / (1L<<(adrW+1))) * 4.0 + 2.0 // 0.0~0.499 -> 2.0~3.999
        } else {
          (n.toDouble / (1L<<(adrW+1))) * 2.0       // 0.5~0.999 -> 1.0~1.999
        }
        val y = round((2.0 / math.sqrt(x)-1.0) * (1L<<fracW))
        if (y >= (1L<<fracW)) {
          println("WARNING: mantissa reaches to 2 while table generation. replaced by 0xFFFF")
          maskL(fracW).U(fracW.W)
        } else if(y < 0.0) { // not used, actually
          0.U(fracW.W)
        } else {
          y.U(fracW.W)
        }
      })
    )
    assert(maxCbit(0) == fracW)

    io.cs.cs(0) := enableIf(io.en, tbl(io.adr(adrW, 0)))

  } else {
    val tableI = InvSqrtSim.invsqrtTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(0 <= diffWidth)
      if(diffWidth != 0) {
        val ci  = coeff(i)
        val msb = ci(cbit(i)-1)
        coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffs.cs(i) := coeff(i)
      }
    }
    io.cs := enableIf(io.en, coeffs)
  }
}
object InvSqrtTableCoeff {
  def getCBits(
    spec:     RealSpec,
    polySpec: PolynomialSpec
  ): Seq[Int] = {

    val order     = polySpec.order
    val adrW      = polySpec.adrW
    val extraBits = polySpec.extraBits
    val fracW     = polySpec.fracW

    if(order == 0) {
      return Seq(fracW)
    } else {
      return InvSqrtSim.invsqrtTableGeneration( order, adrW, spec.manW, fracW ).cbit
    }
  }
  def getCalcW(
    spec:     RealSpec,
    polySpec: PolynomialSpec
  ): Seq[Int] = {

    val order     = polySpec.order
    val adrW      = polySpec.adrW
    val extraBits = polySpec.extraBits
    val fracW     = polySpec.fracW

    if(order == 0) {
      return Seq(fracW)
    } else {
      return InvSqrtSim.invsqrtTableGeneration( order, adrW, spec.manW, fracW ).calcWidth
    }
  }
}

// -------------------------------------------------------------------------
//                        _        _     _                    _   _
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_|
//                                               |_|
// -------------------------------------------------------------------------

// No pathway other than table interpolation. just calculate ex and sgn.
class InvSqrtOtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val x      = Flipped(new DecomposedRealOutput(spec))
    val zother = new RoundingNonTableOutput(spec)
  })

  val xneg = if(spec.disableSign) {false.B} else {io.x.sgn === 1.U(1.W)}

  val znan  = io.x.nan
  val zinf  = io.x.zero || xneg
  val zzero = io.x.inf
  val zman0 = (io.x.man === 0.U) && (!io.x.zero) && (io.x.ex(0) === (exBias % 2).U(1.W))
  // if x.man == 0 && x.ex == 2N, that means x = 2^2N. then z = 2^-N, so zman = 0.

  val zIsNonTable = znan || zinf || zzero || zman0
  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)
  io.zother.znan        := ShiftRegister(znan,        nStage)

  val xExNobias  = io.x.ex - exBias.U
  val xExHalf    = xExNobias(exW-1) ## xExNobias(exW-1, 1) // arith right shift

  val zEx0 = ~xExHalf // -(xex>>1)-1 = ~(xex>>1)+1-1 = ~(xex>>1)
  val zEx  = Mux(zinf || znan, maskU(exW),
             Mux(zzero, 0.U(exW.W), zEx0 + zman0.asUInt + exBias.U))
  val zSgn = 0.U(1.W) // always positive.

  io.zother.zex  := ShiftRegister(zEx , nStage)
  io.zother.zsgn := ShiftRegister(zSgn, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

// use the default RoundingPostProcess.
