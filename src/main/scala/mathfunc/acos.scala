//% @file acos.scala
//
// ACos function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

import java.lang.Math.scalb
import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith._

import rial.math.ACosSim
import rial.mathfunc._

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class ACosPreProcess(
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
  val exAdrW    = ACosSim.calcExAdrW(spec)

  val io = IO(new Bundle {
    val en  = Input (UInt(1.W))
    val x   = Input (UInt(spec.W.W))
    val adr = Output(UInt((exAdrW+adrW).W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  val (xsgn, xex, xman) = FloatChiselUtil.decompose(spec, io.x & Fill(spec.W, io.en))

  val exAdr0 = (exBias - 1).U(exW.W) - xex
  val exAdr  = Mux(exAdr0(exW-1), 0.U(exAdrW.W), exAdr0(exAdrW-1, 0))

  val adr0 = Cat(exAdr, xman(manW-1, dxW))
  val adr  = adr0 & Fill(adr0.getWidth, io.en)
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx0  = Cat(~xman(dxW-1), xman(dxW-2, 0))
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

class ACosTableCoeff(
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
  val exAdrW = ACosSim.calcExAdrW(spec)
  val nStage = stage.total

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input  (UInt((exAdrW+adrW).W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  val exAdr = io.adr(exAdrW + adrW - 1, adrW)
  val adr   = io.adr(adrW - 1, 0)

  val linearThreshold = ACosSim.calcLinearThreshold(manW)

  if(order == 0) {
    val tbl = VecInit( (-1 to linearThreshold.toInt by -1).map( exponent => {
      VecInit( (0L to (1L<<adrW)-1L).map(
        n => {
          val x = scalb(1.0 + n.toDouble/(1L<<adrW), exponent.toInt)
          val y = round(acos(x)*(1L<<fracW))
          y.U((fracW+1).W)
        } ) )
      } ) )
    assert(maxCbit(0) == fracW)

    val c0 = tbl(exAdr)(adr)
    val c  = c0 & Fill(c0.getWidth, io.en)
    io.cs.cs(0) := ShiftRegister(c, nStage) // width should be manW + extraBits

  } else {

    val cbit = ACosSim.acosTableGeneration( order, adrW, manW, fracW )
      .map( t => {t.getCBitWidth(/*sign mode = */0)} )
      .reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

    val tableIs = VecInit(
      ACosSim.acosTableGeneration( order, adrW, manW, fracW ).map(t => {
        t.getVectorWithWidth(cbit, /*sign mode = */ 0)
      })
    )
    val tableI = tableIs(exAdr)
    val coeff = getSlices(tableI(adr), cbit)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      val ci  = if(diffWidth != 0) {
        Cat(0.U(diffWidth.W), coeff(i))
      } else {
        coeff(i) // no need to extend; this is the largest value in all the tables
      }
      coeffs.cs(i) := ci
    }
    val cs = coeffs.asUInt & Fill(coeffs.asUInt.getWidth, io.en)
    io.cs := ShiftRegister(cs.asTypeOf(new TableCoeffInput(maxCbit)), nStage)
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

class ACosNonTableOutput(val spec: RealSpec) extends Bundle {
  val zsgn  = Output(UInt(1.W))
  val zex   = Output(UInt(spec.exW.W))
  val zman  = Output(UInt(spec.manW.W))
  val zIsNonTable = Output(Bool())
  val xsgn  = Output(UInt(1.W))
}

// No pathway other than table interpolation. just calculate ex and sgn.
class ACosOtherPath(
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
    val zother = new ACosNonTableOutput(spec)
  })

  val xMoreThan1 = exBias.U <= io.x.ex

  val znan = io.x.nan
  val zzero = xMoreThan1 && io.x.sgn === 0.U
  val zpi   = xMoreThan1 && io.x.sgn === 1.U

  // --------------------------------------------------------------------------
  // constant

  val constThreshold = -manW
  val isConstant     = io.x.ex < (exBias + constThreshold).U

  val halfpi       = new RealGeneric(spec, Pi * 0.5)
  val zexConstant  = halfpi.ex .U(exW.W)
  val zmanConstant = halfpi.man.toLong.U(manW.W)

  // --------------------------------------------------------------------------
  // linear

  val linearThreshold = ACosSim.calcLinearThreshold(manW)
  val isLinear        = io.x.ex < (exBias + linearThreshold).U

  val xExNobias = io.x.ex.zext - exBias.S

  val constThresholdDigit = log2Up(abs(constThreshold))
  val xmanW1          = io.x.man + (1<<manW).U((4+manW).W)
  val linearExDiff    = (~xExNobias(constThresholdDigit-1, 0)) - 2.U(constThresholdDigit.W)
  val linearXAligned  = Mux(io.x.sgn === 1.U, (xmanW1 >> linearExDiff), ~(xmanW1 >> linearExDiff) + 1.U)
  val linearManDiff   = (((1<<manW) + halfpi.man.toLong) << 3).U((manW+1+3).W) + linearXAligned
  // since pi/2 = 1.1001... and xEx < linearThreshold, mandiff always larger than 1 and never be larger than 2.

  val linearLSB    = linearManDiff(3)
  val linearRound  = linearManDiff(2)
  val linearSticky = linearManDiff(1) | linearManDiff(0)
  val linearInc    = FloatChiselUtil.roundIncBySpec(
    RoundSpec.roundToEven, linearLSB, linearRound, linearSticky)

  assert(linearManDiff(manW+3-1, 3) =/= maskL(manW).U(manW.W))

  val zexLinear  = exBias.U(exW.W)
  val zmanLinear = linearManDiff(manW+3-1, 3) + linearInc

  // --------------------------------------------------------------------------
  // mux

  val pi = new RealGeneric(spec, Pi)

  val zIsNonTable = znan || zzero || zpi || isConstant || isLinear

  val zsgn = 0.U(1.W)

  val zex  = Mux(znan,  Fill(exW, 1.U(1.W)),
             Mux(zzero, 0.U(exW.W),
             Mux(zpi,   pi.ex.U(exW.W),
             Mux(isConstant, zexConstant, zexLinear))))

  val zman = Mux(znan || zzero, Cat(znan, 0.U((manW-1).W)),
             Mux(zpi,           pi.man.toLong.U(manW.W),
             Mux(isConstant,    zmanConstant, zmanLinear)))

  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)
  io.zother.zman        := ShiftRegister(zman, nStage)
  io.zother.zex         := ShiftRegister(zex , nStage)
  io.zother.zsgn        := ShiftRegister(zsgn, nStage)
  io.zother.xsgn        := ShiftRegister(io.x.sgn, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class ACosPostProcess(
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
    val zother = Flipped(new ACosNonTableOutput(spec))
    // table interpolation results
    val zres   = Input(UInt(fracW.W))
    // output
    val z      = Output(UInt(spec.W.W))
  })

  val zIsNonTable  = io.zother.zIsNonTable
  val zsgn         = io.zother.zsgn
  val zexNonTable  = io.zother.zex
  val zmanNonTable = io.zother.zman
  val xsgn         = io.zother.xsgn

  val halfPiFixed = math.round(Pi * 0.5 * (1 << fracW)).U((fracW+1).W)

  val res0  = Cat(io.zres, 0.U(1.W))
  val res   = halfPiFixed + Mux(xsgn === 1.U(1.W), Cat(0.U(1.W), res0), Cat(1.U(1.W), ~res0 + 1.U))

  val shift = (fracW+2).U - (res.getWidth.U - PriorityEncoder(Reverse(res)))
  val resShifted = (res << shift)(fracW+1, 1) - (1<<fracW).U

  val zexTable  = (exBias+1).U(exW.W) - shift
  val zmanTable = (resShifted >> extraBits) + resShifted(extraBits-1)

  val zex  = Mux(zIsNonTable, zexNonTable,  zexTable(exW-1, 0))
  val zman = Mux(zIsNonTable, zmanNonTable, zmanTable(manW-1, 0))

  val z0 = Cat(zsgn, zex, zman)
  val z  = z0 & Fill(z0.getWidth, io.en)

  io.z := ShiftRegister(z, nStage)
}
