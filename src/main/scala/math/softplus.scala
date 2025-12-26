//% @file softplus.scala
//
// square root function
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

// An implementation of the standard softplus function.

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class SoftPlusPreProcessOutput(val spec: RealSpec) extends Bundle {
  val xsgn   = Output(UInt(1.W))
  val xlarge = Output(Bool())
  val zzero  = Output(Bool()) // x = -inf
  val zinf   = Output(Bool()) // x = inf
  val znan   = Output(Bool()) // x = nan
}

class SoftPlusPreProcess(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val dxW       = polySpec.dxW
  val order     = polySpec.order

  val io = IO(new Bundle {
    val en  = Input (UInt(1.W))
    val x   = Flipped(new DecomposedRealOutput(spec))
    val adr = Output(UInt(adrW.W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
    val out  = new SoftPlusPreProcessOutput(spec)
  })

  val rangeMaxLog2   = SoftPlusSim.tableRangeMaxLog2(manW)
  val rangeMaxEx     = exBias + rangeMaxLog2 - 1
  val xExLargeEnough = rangeMaxEx.U < io.x.ex

//   printf("softplus   : rangeMaxLog2   = %d\n", rangeMaxLog2.U)
//   printf("softplus   : rangeMaxEx     = %d\n", rangeMaxEx.U  )
//   printf("softplus   : xExLargeEnough = %b\n", xExLargeEnough)

  val xinf  = io.x.inf && !io.x.nan
  val znan  = io.x.nan
  val zinf  = xinf && io.x.sgn === 0.U
  val zzero = (xinf || xExLargeEnough) && io.x.sgn === 1.U

  val xScaleShift = (rangeMaxEx+1).U(exW.W) - io.x.ex

  val xScaled = Cat(1.U, io.x.man) >> xScaleShift

//   printf("softplus   : xScaled = %b\n", xScaled)

  val adr  = enableIf(io.en, xScaled(manW-1, dxW))
  io.adr := ShiftRegister(adr, nStage)
//   printf("softplus   : adr = %b\n", adr)

  if(order != 0) {
    val dx   = enableIf(io.en, Cat(~xScaled(dxW-1), xScaled(dxW-2, 0)))
    io.dx.get := ShiftRegister(dx, nStage)
//     printf("softplus   : dx = %b\n", dx)
  }

  io.out.xsgn   := ShiftRegister(io.x.sgn,       nStage)
  io.out.xlarge := ShiftRegister(xExLargeEnough, nStage)
  io.out.znan   := ShiftRegister(znan,           nStage)
  io.out.zinf   := ShiftRegister(zinf,           nStage)
  io.out.zzero  := ShiftRegister(zzero,          nStage)
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class SoftPlusTableCoeff(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
) extends Module {

  val manW   = spec.manW
  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val en   = Input(UInt(1.W))
    val adr  = Input(UInt(adrW.W))
    val xsgn = Input(UInt(1.W))
    val cs   = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {

    val rangeMaxLog2 = SoftPlusSim.tableRangeMaxLog2(manW)
    val rangeMax     = 1 << rangeMaxLog2

    val fpos = ( x01: Double ) => {
      val x = x01 * rangeMax // [0, 1) => [0, rangeMax)
      val z = log(1.0 + exp(x)) / (rangeMax*2).toDouble
      assert(z <= 1.0, f"fpos: x = ${x}, z = ${z}")
      z
    }
    val fneg = ( x01: Double ) => {
      val x = -1 * x01 * rangeMax // [0, 1) => [0, rangeMax)
      val z = log(1.0 + exp(x)) // 0 ~ ln(2)
      assert(z <= 1.0, f"fneg: x = ${x}, z = ${z}")
      z
    }

    val tblPos = VecInit( (0L to 1L<<adrW).map(
      n => {
        val x = ( n.toDouble / (1L<<adrW) ) // 0~1
        val y = round( fpos(x) * (1L<<fracW) )
        if (y >= (1L<<fracW)) {
          maskL(fracW).U(fracW.W)
        } else if (y <= 0.0) {
          0.U(fracW.W)
        } else {
          y.toLong.U(fracW.W)
        }
      })
    )
    val tblNeg = VecInit( (0L to 1L<<adrW).map(
      n => {
        val x = ( n.toDouble / (1L<<adrW) ) // 0~1
        val y = round( fneg(x) * (1L<<fracW) )
        if (y >= (1L<<fracW)) {
          maskL(fracW).U(fracW.W)
        } else if (y <= 0.0) {
          0.U(fracW.W)
        } else {
          y.toLong.U(fracW.W)
        }
      })
    )

    val coeffPos = tblPos(io.adr)
    val coeffNeg = tblNeg(io.adr)
    val coeff = Mux(io.xsgn === 0.U, coeffPos, coeffNeg)

    io.cs.cs(0) := enableIf(io.en, coeff)

  } else {

    val tableIs = SoftPlusSim.tableGeneration( order, adrW, manW, fracW )
    val tablePos = tableIs(0)
    val tableNeg = tableIs(1)
    val cbitPos  = tablePos.cbit
    val cbitNeg  = tableNeg.cbit

    val (coeffTablePos, coeffWidthPos) = tablePos.getVectorUnified(/*sign mode =*/0)
    val (coeffTableNeg, coeffWidthNeg) = tableNeg.getVectorUnified(/*sign mode =*/0)
    val coeffPos = getSlices(coeffTablePos(io.adr), coeffWidthPos)
    val coeffNeg = getSlices(coeffTableNeg(io.adr), coeffWidthNeg)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidthPos = maxCbit(i) - cbitPos(i)
      val diffWidthNeg = maxCbit(i) - cbitNeg(i)
      assert(0 <= diffWidthPos)
      assert(0 <= diffWidthNeg)

      val ciPos = Wire(UInt(coeffs.cs(i).getWidth.W))
      val ciNeg = Wire(UInt(coeffs.cs(i).getWidth.W))

      if(diffWidthPos != 0) {
        val ci  = coeffPos(i)
        val msb = ci.head(1)
        // sign extend
        ciPos := Cat(Fill(diffWidthPos, msb), ci)
      } else {
        ciPos := coeffPos(i)
      }

      if(diffWidthNeg != 0) {
        val ci  = coeffNeg(i)
        val msb = ci.head(1)
        // sign extend
        ciNeg := Cat(Fill(diffWidthNeg, msb), ci)
      } else {
        ciNeg := coeffNeg(i)
      }

      coeffs.cs(i) := Mux(io.xsgn === 0.U, ciPos, ciNeg)
//       printf(f"softplus   : coeff(${i}) = %%d(%%b)\n", coeffs.cs(i).asSInt,  coeffs.cs(i))
    }
    io.cs := enableIf(io.en, coeffs)
  }
}

object SoftPlusTableCoeff {
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
      return SoftPlusSim.tableGeneration( order, adrW, spec.manW, fracW ).
        map(_.cbit).reduce( (lhs, rhs) => {
          lhs.zip(rhs).map( x => max(x._1, x._2) )
        } )
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
      return SoftPlusSim.tableGeneration( order, adrW, spec.manW, fracW ).
        map(_.calcWidth).reduce( (lhs, rhs) => {
          lhs.zip(rhs).map( x => max(x._1, x._2) )
        } )
    }
  }
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class SoftPlusPostProcess(
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

  val calcW  = manW + extraBits

  val io = IO(new Bundle {
    val en     = Input(UInt(1.W))
    val preout = Flipped(new SoftPlusPreProcessOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val x      = Flipped(new DecomposedRealOutput(spec))
    val z      = Output(UInt(spec.W.W))
  })

//   printf("softplus   : zres = %b\n", io.zres)

  val xsgn  = io.preout.xsgn
  val zsgn  = 0.U(1.W)

  val xlarge= io.preout.xlarge
  val znan  = io.preout.znan
  val zinf  = io.preout.zinf
  val zzero = io.preout.zzero
  assert((znan ^ zinf) || (!znan && !zinf),
         "znan and zinf cannot be true at the same time")

  val zres = io.zres

  val zShift = Mux(zres === 0.U, zres.getWidth.U, PriorityEncoder(Reverse(zres))) + 1.U
  val zmanW1 = io.zres << zShift
  val zmanW1Rounded = (zmanW1 >> extraBits) + zmanW1(extraBits-1)

//   printf("softplus   : zres   = %b\n", zres)
//   printf("softplus   : zShift = %d\n", zShift)
//   printf("softplus   : zmanW1 = %b\n", zmanW1)
//   printf("softplus   : zmanW1Rounded = %b\n", zmanW1Rounded)

  val zmanW1RoundedMoreThan1 = zmanW1Rounded >= (BigInt(1) << (manW+1)).U
  val zman = Mux(zmanW1RoundedMoreThan1, Fill(manW, 1.U(1.W)), zmanW1Rounded(manW-1, 0))

//   printf("softplus   : zman = %b\n", zman)

  val rangeMaxLog2 = SoftPlusSim.tableRangeMaxLog2(manW)
  val zex = exBias.U(exW.W) - zShift + enableIf(xsgn === 0.U, (rangeMaxLog2 + 1).U)

  val z0 =
    Mux(znan || zinf, Cat(zsgn, Fill(exW, 1.U(1.W)), znan.asUInt, 0.U((manW-1).W)),
    Mux(zzero,        0.U(spec.W.W),
    Mux(xlarge,       Cat(0.U(1.W), io.x.ex, io.x.man),
    /* default = */   Cat(zsgn, zex, zman))))

  val z = enableIf(io.en, z0)
  io.z := ShiftRegister(z, nStage)
}
