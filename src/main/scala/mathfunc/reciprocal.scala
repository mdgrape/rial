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
    val x   = Flipped(new DecomposedRealOutput(spec))
    val adr = Output(UInt(adrW.W))
    val dx  = if (order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  // invert x to make all the polynomial coefficients positive
  val adr  = enable(io.en, ~io.x.man(manW-1, dxW))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx   = enable(io.en, Cat(io.x.man(dxW-1), ~io.x.man(dxW-2, 0)))
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
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
) extends Module {

  val manW   = spec.manW
  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input(UInt(adrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {
    val tbl = VecInit( (0L to (1L<<adrW)-1L).map(
      n => {
        val x = 1.0+n.toDouble/(1L<<adrW)
        val y = round((2.0/x-1.0)*(1L<<fracW))
        //println(f"$n $x $y")
        if (n==0) {
          0.U(fracW.W)
        } else if (y>=(1L<<fracW)) {
          println("WARNING: mantissa reaches to 2")
          maskL(fracW).U(fracW.W)
        } else {
          y.U(fracW.W)
        }
      }
    ) )

    assert(maxCbit(0) == fracW)

    io.cs.cs(0) := enable(io.en, tbl(~io.adr))

  } else {
    val tableI = ReciprocalSim.reciprocalTableGeneration( order, adrW, manW, fracW )
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
    io.cs := enable(io.en, coeffs)
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
      return ReciprocalSim.reciprocalTableGeneration( order, adrW, spec.manW, fracW ).calcWidth
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
  val znan  = Output(UInt(1.W))
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
  val zIsNonTable = znan || zinf || zzero

  io.zother.znan := ShiftRegister(znan, nStage)
  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)

  val zsgn = Mux(znan || zzero || zinf, 0.U, io.x.sgn) // TODO sign
  val zex  = Mux(znan || zinf, maskU(exW),
             Mux(zzero, 0.U(exW.W), zex0(exW-1,0)))

  io.zother.zsgn := ShiftRegister(zsgn, nStage)
  io.zother.zex  := ShiftRegister(zex , nStage)
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
  def getStage = nStage

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
  val znan = io.zother.znan
  val zIsNonTable  = io.zother.zIsNonTable
  val zmanNonTable = Cat(znan, Fill(manW-1, 0.U(1.W)))

  val zmanRounded = Wire(UInt(manW.W))
  if(extraBits == 0) {
    zmanRounded := io.zres
  } else {
    val zman0 = dropLSB(extraBits, io.zres) +& io.zres(extraBits-1)
    val polynomialOvf = zman0(manW)
    zmanRounded := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW-1,0))
  }

  val zman = Mux(zIsNonTable, zmanNonTable, zmanRounded)
  val z = enable(io.en, Cat(zsgn, zex, zman))

  io.z   := ShiftRegister(z, nStage)
}

// -------------------------------------------------------------------------
//                      _     _                _
//   ___ ___  _ __ ___ | |__ (_)_ __   ___  __| |
//  / __/ _ \| '_ ` _ \| '_ \| | '_ \ / _ \/ _` |
// | (_| (_) | | | | | | |_) | | | | |  __/ (_| |
//  \___\___/|_| |_| |_|_.__/|_|_| |_|\___|\__,_|
// -------------------------------------------------------------------------

class ReciprocalGeneric(
  val spec     : RealSpec,
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage    : MathFuncPipelineConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val pcGap = if(stage.preCalcGap ) {1} else {0}
  val cpGap = if(stage.calcPostGap) {1} else {0}

  val nPreStage  = stage.preStage.total
  val nCalcStage = stage.calcStage.total
  val nPostStage = stage.postStage.total

  val nStage   = stage.total
  def getStage = nStage

  val polySpec = new PolynomialSpec(spec, nOrder, adrW, extraBits,
    enableRangeCheck, enablePolynomialRounding)
  val order = polySpec.order

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val cbits = ReciprocalTableCoeff.getCBits(spec, polySpec)
  val calcW = ReciprocalTableCoeff.getCalcW(spec, polySpec)

  def getCbit  = cbits
  def getCalcW = calcW

  val io = IO(new Bundle {
    val en = Input(Bool())
    val x = Input (UInt(spec.W.W))
    val z = Output(UInt(spec.W.W))
  })

  // --------------------------------------------------------------------------

  val xdecomp = Module(new DecomposeReal(spec))
  xdecomp.io.real := io.x

  val enPCReg    = ShiftRegister(io.en,      nPreStage)
  val enPCGapReg = ShiftRegister(enPCReg,    pcGap)
  val enCPReg    = ShiftRegister(enPCGapReg, nCalcStage)
  val enCPGapReg = ShiftRegister(enCPReg,    cpGap)

  val xdecPCReg    = ShiftRegister(xdecomp.io.decomp, nPreStage)
  val xdecPCGapReg = ShiftRegister(xdecPCReg,         pcGap)
  val xdecCPReg    = ShiftRegister(xdecPCGapReg,      nCalcStage)
  val xdecCPGapReg = ShiftRegister(xdecCPReg,         cpGap)

  // --------------------------------------------------------------------------

  val recPre   = Module(new ReciprocalPreProcess (spec, polySpec, stage.preStage))
  val recTab   = Module(new ReciprocalTableCoeff (spec, polySpec, cbits))
  val recOther = Module(new ReciprocalOtherPath  (spec, polySpec, stage.calcStage))
  val recPost  = Module(new ReciprocalPostProcess(spec, polySpec, stage.postStage))

  val recPreAdrPCGapReg = ShiftRegister(recPre.io.adr, pcGap)

  recPre.io.en  := io.en
  recPre.io.x   := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  recTab.io.en  := enPCGapReg
  recTab.io.adr := recPreAdrPCGapReg
  recOther.io.x := xdecPCGapReg

  // after preprocess
  assert(recPre.io.adr === 0.U               || enPCReg)
  assert(recPre.io.dx.getOrElse(0.U) === 0.U || enPCReg)
  assert(recTab.io.cs.asUInt === 0.U         || enPCGapReg)

  // --------------------------------------------------------------------------

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, cbits, stage.calcStage))

  if(order != 0) {
    polynomialEval.io.dx.get := ShiftRegister(recPre.io.dx.get, pcGap)
  }
  polynomialEval.io.coeffs.cs := recTab.io.cs.cs

  val polynomialResultCPGapReg = ShiftRegister(polynomialEval.io.result, cpGap)

  recPost.io.en     := enCPGapReg
  recPost.io.zother := ShiftRegister(recOther.io.zother, cpGap)
  recPost.io.zres   := polynomialResultCPGapReg

  io.z := recPost.io.z
}


