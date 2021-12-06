//% @file sinPi.scala
//
// x -> sin(pi * x)
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import chisel3._
import chisel3.util._

import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.FloatChiselUtil
import rial.arith._

import rial.math.SinPiSim
import rial.mathfunc._

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class SinPiPreProcessOutput(val spec: RealSpec) extends Bundle {
  val xConvertedEx  = Output(UInt(spec.manW.W))
  val xConvertedMan = Output(UInt(spec.manW.W))
}

class SinPiPreProcess(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias
  val exAdrW = SinPiSim.calcExAdrW(spec)

  val io = IO(new Bundle {
    val x   = Input (UInt(spec.W.W))
    val adr = Output(UInt((exAdrW+adrW).W))
    val dx  = Output(UInt((manW-adrW).W))
    val xConverted = new SinPiPreProcessOutput(spec)
  })

  val (xsgn, xex, xman) = FloatChiselUtil.decompose(spec,  io.x)

  // --------------------------------------------------------------------------
  // convert everything into [0, 0.5) (ex = [-inf to -2])

  // if ex == 0, x is in [1, 2)
  // if ex ==-1, x is in [1/2, 1)
  val manExceedsHalf = xman(manW-1)
  val halfpos  = xman
  val halfneg  = ~xman + 1.U
  val shiftpos = Mux(halfpos === 0.U, 0.U, PriorityEncoder(Reverse(halfpos)) + 1.U)
  val shiftneg = Mux(halfneg === 0.U, 0.U, PriorityEncoder(Reverse(halfneg)) + 1.U)
  assert(shiftpos <= (manW+1).U)
  assert(shiftneg <= (manW+1).U)

  val halfEx0  = Mux(manExceedsHalf, halfneg , halfpos)
  val shiftEx0 = Mux(manExceedsHalf, shiftneg, shiftpos)

  val yman = MuxCase(xman, Seq(
    (xex === (exBias+0).U) -> (halfEx0 << shiftEx0)(manW-1, 0),
    (xex === (exBias-1).U) -> (halfneg << shiftneg)(manW-1, 0)
  ))

  // shifts never exceeds manW. we can assume exBias-1-shift > 0.
  val yex  = MuxCase(xex, Seq(
      (xex === (exBias  ).U) -> Mux(halfEx0 === 0.U, 0.U(exW.W), exBias.U-shiftEx0),
      (xex === (exBias-1).U) -> ((exBias-1).U - shiftneg)
    ))

  val exOfs = (exBias-2).U - yex
  val exAdr = exOfs.asUInt()(exAdrW-1, 0)

  val adr0 = Cat(exAdr, yman(manW-1, manW-adrW)) // concat exAdr + man
  val dr0  = Cat(~yman(manW-adrW-1), yman(manW-adrW-2,0))

  io.adr   := ShiftRegister(adr0,  nStage)
  io.dx    := ShiftRegister(dr0,   nStage)
  io.xConverted.xConvertedEx  := ShiftRegister(yex,  nStage)
  io.xConverted.xConvertedMan := ShiftRegister(yman, nStage)
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class SinPiTableCoeff(

  val spec : RealSpec,
  val nOrder: Int, val adrW : Int, val extraBits : Int,

  val maxAdrW : Int,      // max address width among all math funcs
  val maxCbit : Seq[Int], // max coeff width among all math funcs

  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val manW   = spec.manW
  val calcW  = manW + extraBits
  val order  = if(adrW == manW) {0} else {nOrder}
  val nStage = stage.total

  val linearThreshold = SinPiSim.calcLinearThreshold(manW)

  val io = IO(new Bundle {
    val adr = Input  (UInt(maxAdrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  val exAdrW = SinPiSim.calcExAdrW(spec, /*cubic =*/ false) // TODO allow cubic
  val exAdr  = io.adr(adrW + exAdrW-1, adrW)
  val manAdr = io.adr(adrW-1, 0)

  if(order == 0) {
    val tbl = VecInit( (-2 to linearThreshold by -1).map( exponent => {
      VecInit((0L to (1L << adrW)).map(
        n => {
          val x = n.toDouble / (1L<<adrW)
          val y = round(scalb(math.sin(Pi * scalb(1.0 + x, exponent)), -exponent-3) * (1L<<calcW))
          assert(y < (1L<<calcW))
          y.U(calcW.W)
        }))
      })
    )

    assert(maxCbit(0) == calcW)

    val c0 = tbl(exAdr)(manAdr)
    io.cs.cs(0) := ShiftRegister(c0, nStage) // width should be manW + extraBits

  } else {

    val cbit = SinPiSim.sinPiTableGeneration( nOrder, adrW, manW, calcW )
      .map( t => {t.getCBitWidth(/*sign mode = */0)} )
      .reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

    val tableIs = VecInit(
      SinPiSim.sinPiTableGeneration( nOrder, adrW, manW, calcW ).map(t => {
        t.getVectorWithWidth(cbit, /*sign mode = */ 0)
      })
    )
    val tableI = tableIs(exAdr)

    val coeff = getSlices(tableI(manAdr), cbit)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to nOrder) {
      val diffWidth = maxCbit(i) - cbit(i)
      val ci  = coeff(i)
      val msb = ci(cbit(i)-1)
      if(0 < diffWidth) {
        coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffs.cs(i) := ci
      }
    }
    io.cs := ShiftRegister(coeffs, nStage)
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

class SinPiNonTableOutput(val spec: RealSpec) extends Bundle {
  val zIsNonTable = Output(Bool())
  val zsgn        = Output(UInt(1.W))
  val zex         = Output(UInt(spec.exW.W))
  val zman        = Output(UInt(spec.manW.W))
}

class SinPiOtherPath(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val x          = Input(UInt(spec.W.W))
    val xConverted = Flipped(new SinPiPreProcessOutput(spec))
    val zother     = new SinPiNonTableOutput(spec)
  })

  val (xsgn,  xex,  xman) = FloatChiselUtil.decompose(spec, io.x)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(spec, io.x)

  val yex  = io.xConverted.xConvertedEx
  val yman = io.xConverted.xConvertedMan

  val out_of_bounds = ((exBias+1).U <= xex) && (xman =/= 0.U) // 2 < |x|
  val znan  = xnan || xinf || out_of_bounds
  val zzero = ((exBias.U <= xex)    && xman === 0.U) ||
              (yex === 0.U          && yman === 0.U)
  val zone  = (yex === (exBias-1).U && yman === 0.U)

  val zSgn  = !zzero && ((xsgn === 1.U) || (xex === exBias.U))

  // --------------------------------------------------------------------------
  // linear approximation around zero

  val linearThreshold = SinPiSim.calcLinearThreshold(manW)
  val pi = new RealGeneric(spec, Pi)

  val isLinear = (yex < (linearThreshold + exBias).U)

  // y is in [0, 0.5) so it never overflows
  val linearProdEx        = (pi.ex-exBias).U(exW.W) + yex
  val linearProdMan       = (pi.man + (1<<manW)).toLong.U((manW+1).W) *
                            (yman   + (1<<manW).U((manW+1).W))
  val linearProdbp        = manW + manW
  val linearProdMoreThan2 = linearProdMan(linearProdbp+1)
  val linearProdRoundBits = linearProdbp - manW

  val linearProdSticky    = linearProdMan(linearProdRoundBits-2, 0).orR |
                           (linearProdMoreThan2 & linearProdMan(linearProdRoundBits-1))
  val linearProdRound     = Mux(linearProdMoreThan2, linearProdMan(linearProdRoundBits),
                                                     linearProdMan(linearProdRoundBits-1))
  val linearProdShift     = Mux(linearProdMoreThan2, linearProdMan(linearProdbp, linearProdRoundBits+1),
                                                     linearProdMan(linearProdbp-1, linearProdRoundBits))
  val linearProdLSB       = linearProdShift(0)

  val linearProdInc       = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
                                linearProdLSB, linearProdRound, linearProdSticky)
  val linearProdRounded   = linearProdShift +& linearProdInc
  val linearProdMoreThan2AfterRound = linearProdRounded(manW)

  val zManLinear = linearProdRounded(manW-1, 0)
  val zExLinear  = linearProdEx + (linearProdMoreThan2 | linearProdMoreThan2AfterRound)

  // TODO: cubic

  val zIsNonTable = znan || zzero || zone || isLinear
  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)

  val zEx = Mux(isLinear, zExLinear,
            Mux(znan,     Fill(exW, 1.U(1.W)),
            Mux(zzero,    0.U,
            Mux(zone,     exBias.U,
                          yex))))   // default

  val zeroFlush = znan || zzero || zone
  val zMan = Mux(zeroFlush, Cat(znan, 0.U((manW-1).W)), zManLinear)

  io.zother.zsgn := ShiftRegister(zSgn, nStage)
  io.zother.zex  := ShiftRegister(zEx,  nStage)
  io.zother.zman := ShiftRegister(zMan, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class SinPiPostProcess(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val nStage = stage.total
  def getStage() = nStage

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val bp = manW + extraBits

  val io = IO(new Bundle {
    val zother = Flipped(new SinPiNonTableOutput(spec))
    val zres   = Input(UInt((manW+extraBits).W))
    val z      = Output(UInt(spec.W.W))
  })

  val resLessThanHalf = io.zres(bp-1) === 0.U

  val zEx0Table  = Mux(resLessThanHalf, io.zother.zex + 1.U, io.zother.zex + 2.U)
  val zExTable   = zEx0Table(exW-1, 0)

  val zMan0Table = Mux(resLessThanHalf, Cat(io.zres, 0.U(2.W))(bp-1, 0),
                                        Cat(io.zres, 0.U(1.W))(bp-1, 0))
  val zManTable  = zMan0Table(bp-1, bp-manW) + zMan0Table(bp-manW-1)

  val zSgn = io.zother.zsgn
  val zEx  = Mux(io.zother.zIsNonTable, io.zother.zex,  zExTable)
  val zMan = Mux(io.zother.zIsNonTable, io.zother.zman, zManTable)

  val z0  = Cat(zSgn, zEx, zMan)
  io.z   := ShiftRegister(z0, nStage)
}
