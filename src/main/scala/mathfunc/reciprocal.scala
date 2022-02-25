//% @file reciprocal.scala
//
// Reciprocal function
// Copyright (C) Toru Niina RIKEN BDR 2020
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

import rial.math.ReciprocalSim
import rial.mathfunc._

// 1/x, floating input, floating output
// x   = 2^e * 1.m
// 1/x = 2^(-e-1) * 2/1.m

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class ReciprocalPreProcess(
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
    val adr = Output(UInt(adrW.W))
    val dx  = if (order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  // invert x to make all the polynomial coefficients positive
  val adr0 = ~io.x(manW-1, dxW)
  val adr  = adr0 & Fill(adr0.getWidth, io.en)
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx0  = Cat(io.x(dxW-1), ~io.x(dxW-2, 0))
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

class ReciprocalTableCoeff(
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
    val adr = Input(UInt(adrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {
    val tbl = VecInit( (0L to (1L<<adrW)-1L).map(
      n => {
        val x = 1.0+n.toDouble/(1L<<adrW)
        val y = round((2.0/x-1.0)*(1L<<manW))
        //println(f"$n $x $y")
        if (n==0) {
          0.U(manW.W)
        } else if (y>=(1L<<manW)) {
          println("WARNING: mantissa reaches to 2")
          maskL(manW).U(manW.W)
        } else {
          y.U(manW.W)
        }
      }
    ) )

    assert(maxCbit(0) == fracW)

    val c0 = tbl(~io.adr) // address is inverted in preprocess
    val c  = c0 & Fill(c0.getWidth, io.en)
    io.cs.cs(0) := ShiftRegister(c, nStage) // width should be manW + extraBits

  } else {
    val tableI = ReciprocalSim.reciprocalTableGeneration( order, adrW, manW, fracW )
    val cbit   = tableI.cbit

    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      val ci  = coeff(i)
      val msb = ci(cbit(i)-1)
      coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
    }
    val cs = coeffs.asUInt & Fill(coeffs.asUInt.getWidth, io.en)
    io.cs := ShiftRegister(cs.asTypeOf(new TableCoeffInput(maxCbit)), nStage)
  }
}
object ReciprocalTableCoeff {
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
      return ReciprocalSim.reciprocalTableGeneration( order, adrW, spec.manW, fracW ).cbit
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

class ReciprocalNonTableOutput(val spec: RealSpec) extends Bundle {
  val zsgn  = Output(UInt(1.W))
  val zex   = Output(UInt(spec.exW.W))
  val zman  = Output(UInt(spec.manW.W))
  val zIsNonTable = Output(Bool())
}

// No pathway other than table interpolation. just calculate ex and sgn.
class ReciprocalOtherPath(
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
    val zother = new ReciprocalNonTableOutput(spec)
  })

  val xmanNonZero = io.x.man.orR.asUInt
  val zex0        = (exBias<<1).U((exW+1).W) - io.x.ex - xmanNonZero
  val zexMSB      = zex0(exW)
  // add 1 bit to check carry

  // check overflow/underflow
  // 1/x = 2^(-e-1) * 2/1.m
  val invExMax = -spec.exMin - 1
  val invExMin = -spec.exMax - 1
  val (zinf, zzero) = if(spec.exMax < invExMax) { // check overflow
    (zexMSB || io.x.zero, io.x.inf || zex0 === 0.U)
  } else if (invExMin < spec.exMin) { // check underflow
    (io.x.zero, zexMSB || io.x.inf || zex0 === 0.U)
  } else { // nothing happens.
    (io.x.zero, io.x.inf || zex0 === 0.U)
  }
  val znan  = io.x.nan

  val zsgn = Mux(znan || zzero || zinf, 0.U, io.x.sgn) // TODO sign
  val zex  = Mux(znan || zinf, maskU(exW),
             Mux(zzero, 0.U(exW.W), zex0(exW-1,0)))
  val zman = Cat(znan, 0.U((manW-1).W))
  val zIsNonTable = znan || zinf || zzero

  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)
  io.zother.zsgn := ShiftRegister(zsgn, nStage)
  io.zother.zex  := ShiftRegister(zex , nStage)
  io.zother.zman := ShiftRegister(zman, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class ReciprocalPostProcess(
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
    val zother = Flipped(new ReciprocalNonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
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
