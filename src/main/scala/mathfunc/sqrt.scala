//% @file sqrt.scala
//
// square root function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

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

import rial.math.SqrtSim
import rial.mathfunc._

// An implementation of sqrt using (normally 2nd-order) polynomial interpolation
// Since FP is composed of exponent and mantissa, the result of sqrt will be
//     sqrt(2^ex * 1.m) = 2^(ex/2) * sqrt(1.m) .
// The fractional part of exponent appears when ex is an odd number, so
//     sqrt(2^ex * 1.m) = 2^(floor(ex/2)) * sqrt(1.m)   if ex is even
//                        2^(floor(ex/2)) * sqrt(2*1.m) otherwise.
// It requires 2 different tables for sqrt(1.m) and sqrt(2*1.m).
// To distinguish them, we use LSB of ex to detect whether ex is even or odd.

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

// sqrt(x): floating => floating
// - if x < 0, returns 0.
class SqrtPreProcess(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val manW = spec.manW

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val io = IO(new Bundle {
    val en  = Input (UInt(1.W))
    val x   = Input (UInt(spec.W.W))
    val adr = Output(UInt((1+adrW).W)) // we use LSB of x.ex
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  val adr0 = io.x(manW, dxW)
  val adr  = adr0 & Fill(adr0.getWidth, io.en)
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx0  = Cat(~io.x(dxW-1), io.x(dxW-2, 0))
    val dx   = dx0 & Fill(dx0.getWidth, io.en)
    io.dx.get := ShiftRegister(dx, nStage)
  }
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class SqrtTableCoeff(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val maxAdrW  : Int,      // max address width among all math funcs
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
  val stage    : PipelineStageConfig,
) extends Module {

  val manW   = spec.manW
  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order
  val nStage = stage.total

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input(UInt((1+adrW).W))
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
        val y = round((math.sqrt(x)-1.0) * (1L<<manW))
        if (y >= (1L<<manW)) {
          maskL(manW).U(manW.W)
        } else if (y <= 0.0) {
          0.U(manW.W)
        } else {
          y.U(manW.W)
        }
      })
    )
    assert(maxCbit(0) == fracW)

    val c0 = tbl(io.adr(adrW, 0))            // here we use LSB of ex
    val c  = c0 & Fill(c0.getWidth, io.en)
    io.cs.cs(0) := ShiftRegister(c, nStage) // width should be manW + extraBits

  } else {
    val tableI = SqrtSim.sqrtTableGeneration( order, adrW, manW, fracW )
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
    val cs = coeffs.asUInt & Fill(coeffs.asUInt.getWidth, io.en)
    io.cs := ShiftRegister(cs.asTypeOf(new TableCoeffInput(maxCbit)), nStage)
  }
}
object SqrtTableCoeff {
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
      return SqrtSim.sqrtTableGeneration( order, adrW, spec.manW, fracW ).cbit
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

class SqrtNonTableOutput(val spec: RealSpec) extends Bundle {
  val zsgn  = Output(UInt(1.W))
  val zex   = Output(UInt(spec.exW.W))
  val zman  = Output(UInt(spec.manW.W))
  val zIsNonTable = Output(Bool())
}

// No pathway other than table interpolation. just calculate ex and sgn.
class SqrtOtherPath(
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
    val zother = new SqrtNonTableOutput(spec)
  })

  val xneg = if(spec.disableSign) {false.B} else {io.x.sgn === 1.U(1.W)}

  val znan  = io.x.nan
  val zinf  = io.x.inf
  val zzero = io.x.zero || xneg

  val zMan  = Cat(znan, 0.U((manW-1).W)) // in case of znan, zinf, or zzero
  val zIsNonTable = znan || zinf || zzero
  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)
  io.zother.zman        := ShiftRegister(zMan,        nStage)

  val xExNobias  = io.x.ex - exBias.U
  val zExShifted = xExNobias(exW-1) ## xExNobias(exW-1, 1) // arith right shift

  val zex = Mux(zinf || znan, maskU(exW),
            Mux(zzero,        0.U(exW.W),
                              zExShifted + exBias.U))
  val zsgn = 0.U(1.W)

  io.zother.zex  := ShiftRegister(zex , nStage)
  io.zother.zsgn := ShiftRegister(zsgn, nStage) // always positive.
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class SqrtPostProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val nStage = stage.total
  def getStage() = nStage

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order
  val extraBits = polySpec.extraBits

  val io = IO(new Bundle {
    val en = Input(UInt(1.W))
    // ex and some flags
    val zother = Flipped(new SqrtNonTableOutput(spec))
    // table interpolation results
    val zres   = Input(UInt(fracW.W))
    // output
    val z      = Output(UInt(spec.W.W))
  })

  val zsgn = io.zother.zsgn
  val zex  = io.zother.zex
  val zmanNonTable = io.zother.zman
  val zIsNonTable  = io.zother.zIsNonTable

  val zman0 = dropLSB(extraBits, io.zres) +& io.zres(extraBits-1)
  val polynomialOvf = zman0(manW)
  val zmanRounded   = Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW-1,0))
  val zman          = Mux(zIsNonTable, zmanNonTable, zmanRounded)

  val z0 = Cat(zsgn, zex, zman)
  val z = z0 & Fill(z0.getWidth, io.en)

  io.z   := ShiftRegister(z, nStage)
}
