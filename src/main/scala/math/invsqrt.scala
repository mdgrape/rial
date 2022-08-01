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

    io.cs.cs(0) := enable(io.en, tbl(io.adr(adrW, 0)))

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
    io.cs := enable(io.en, coeffs)
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

class InvSqrtNonTableOutput(val spec: RealSpec) extends Bundle {
  val zsgn  = Output(UInt(1.W))
  val zex   = Output(UInt(spec.exW.W))
  val znan  = Output(Bool())
  val zIsNonTable = Output(Bool())
}

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
    val zother = new InvSqrtNonTableOutput(spec)
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

class InvSqrtPostProcess(
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
    val en  = Input(UInt(1.W))
    // ex and some flags
    val zother = Flipped(new InvSqrtNonTableOutput(spec))
    // table interpolation results
    val zres   = Input(UInt(fracW.W))
    // output
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
    zmanRounded := Mux(polynomialOvf, maskU(manW), zman0(manW-1,0))
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

class InvSqrtGeneric(
  val spec     : RealSpec,
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage    : MathFuncPipelineConfig,
  val dxW0 : Option[Int] = None,
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

  val polySpec = new PolynomialSpec(spec, nOrder, adrW, extraBits, dxW0,
    enableRangeCheck, enablePolynomialRounding)
  val order = polySpec.order

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val cbits = InvSqrtTableCoeff.getCBits(spec, polySpec)
  val calcW = InvSqrtTableCoeff.getCalcW(spec, polySpec)

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

  val invsqrtPre   = Module(new SqrtPreProcess    (spec, polySpec, stage.preStage))
  val invsqrtTab   = Module(new InvSqrtTableCoeff (spec, polySpec, cbits))
  val invsqrtOther = Module(new InvSqrtOtherPath  (spec, polySpec, stage.calcStage))
  val invsqrtPost  = Module(new InvSqrtPostProcess(spec, polySpec, stage.postStage))

  val invsqrtPreAdrPCGapReg = ShiftRegister(invsqrtPre.io.adr, pcGap)

  invsqrtPre.io.en  := io.en
  invsqrtPre.io.x   := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  invsqrtTab.io.en  := enPCGapReg
  invsqrtTab.io.adr := invsqrtPreAdrPCGapReg
  invsqrtOther.io.x := xdecPCGapReg

  // after preprocess
  assert(invsqrtPre.io.adr === 0.U               || enPCReg)
  assert(invsqrtPre.io.dx.getOrElse(0.U) === 0.U || enPCReg)
  assert(invsqrtTab.io.cs.asUInt === 0.U         || enPCGapReg)

  // --------------------------------------------------------------------------

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, cbits, stage.calcStage))

  if(order != 0) {
    polynomialEval.io.dx.get := ShiftRegister(invsqrtPre.io.dx.get, pcGap)
  }
  polynomialEval.io.coeffs.cs := invsqrtTab.io.cs.cs

  val polynomialResultCPGapReg = ShiftRegister(polynomialEval.io.result, cpGap)

  invsqrtPost.io.en     := enCPGapReg
  invsqrtPost.io.zother := ShiftRegister(invsqrtOther.io.zother, cpGap)
  invsqrtPost.io.zres   := polynomialResultCPGapReg

  io.z := invsqrtPost.io.z
}
