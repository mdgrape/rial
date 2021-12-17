//% @file acos.scala
//
// ATan2 function
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

import rial.math.ATan2Sim
import rial.mathfunc._

// ATan2 Stage1 calculates min(x,y)/max(x,y).
//       Stage2 calculates atan(min(x,y)/max(x,y)) +/- constant.
// Some flags are needed to be saved.

object ATan2Status {
  val W = 2
  val xIsPosIsLarger  = 0.U(W.W) // ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|
  val xIsNegIsLarger  = 1.U(W.W) // ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y|
  val xIsPosIsSmaller = 2.U(W.W) // ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y|
  val xIsNegIsSmaller = 3.U(W.W) // ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y|
}
object ATan2SpecialValue {
  val W = 3
  val zNormal     = 0.U(W.W)
  val zNaN        = 1.U(W.W) // == nan
  val zZero       = 2.U(W.W) // == zero
  val zPi         = 3.U(W.W) // == pi
  val zHalfPi     = 4.U(W.W) // == pi/2
  val zQuarterPi  = 5.U(W.W) // == pi/4
  val z3QuarterPi = 6.U(W.W) // == 3pi/4
}
class ATan2Flags extends Bundle {
  val status  = UInt(ATan2Status.W.W)
  val special = UInt(ATan2SpecialValue.W.W)
  val ysgn    = UInt(1.W)
}

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------
//
// Stage1 preprocess only checks the special cases.
// Stage1 calculates min(x,y)/max(x,y), so ReciprocalPreProcess is re-used.
//
class ATan2Stage1PreProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val dxW    = polySpec.dxW
  val order  = polySpec.order
  val exAdrW = ATan2Sim.calcExAdrW(spec)

  val io = IO(new Bundle {
    val x       = Input (new DecomposedRealOutput(spec))
    val y       = Input (new DecomposedRealOutput(spec))
    val special = Output(UInt(ATan2SpecialValue.W.W))
  })

  val xpos      = io.x.sgn === 0.U
  val xneg      = io.x.sgn === 1.U
  val znan      =  (io.x.nan ||  io.y.nan) || ( io.x.zero &&  io.y.zero)
  val zzero     = ((io.x.inf && !io.y.inf) || (!io.x.zero &&  io.y.zero)) && xpos
  val zpi       = ((io.x.inf && !io.y.inf) || (!io.x.zero &&  io.y.zero)) && xneg
  val zhalfpi   = (!io.x.inf &&  io.y.inf) || ( io.x.zero && !io.y.zero)
  val z1piover4 =  (io.x.inf &&  io.y.inf) &&  xpos
  val z3piover4 =  (io.x.inf &&  io.y.inf) &&  xneg

  val special = MuxCase(ATan2SpecialValue.zNormal, Seq(
    znan      -> ATan2SpecialValue.zNaN,
    zzero     -> ATan2SpecialValue.zZero,
    zpi       -> ATan2SpecialValue.zPi,
    zhalfpi   -> ATan2SpecialValue.zHalfPi,
    z1piover4 -> ATan2SpecialValue.zQuarterPi,
    z3piover4 -> ATan2SpecialValue.z3QuarterPi
  ))

  io.special := ShiftRegister(special, nStage)
}

// Stage2 calculates atan(x), so we need to extract the address value

class ATan2Stage2PreProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val dxW    = polySpec.dxW
  val order  = polySpec.order
  val exAdrW = ATan2Sim.calcExAdrW(spec)

  val io = IO(new Bundle {
    val x   = Input (UInt(spec.W.W))
    val adr = Output(UInt((exAdrW+adrW).W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  val (xsgn, xex, xman) = FloatChiselUtil.decompose(spec, io.x)

  val exAdr0 = (exBias - 1).U(exW.W) - xex
  val exAdr  = exAdr0(exAdrW-1, 0)

  val adr0 = Cat(exAdr, xman(manW-1, dxW))
  io.adr := ShiftRegister(adr0, nStage)

  if(order != 0) {
    val dx0  = Cat(~xman(manW-adrW-1), xman(manW-adrW-2,0))
    io.dx.get := ShiftRegister(dx0, nStage)
  }
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------
//
// Stage1 re-use the reciprocal table. we don't need to implement it for atan2.
//
class ATan2Stage2TableCoeff(
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
  val exAdrW = ATan2Sim.calcExAdrW(spec)
  val nStage = stage.total

  val io = IO(new Bundle {
    val adr = Input  (UInt((exAdrW+adrW).W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  val exAdr = io.adr(exAdrW + adrW - 1, adrW)
  val adr   = io.adr(adrW - 1, 0)

  val linearThreshold = ATan2Sim.calcLinearThreshold(manW)

  if(order == 0) {

    val tbl = VecInit( (-1 to linearThreshold.toInt by -1 ).map( exponent => {
      VecInit((0L to (1L<<adrW)-1L).map( n => {
        // atan(x) < x.
        val x = scalb(1.0 + n.toDouble / (1L<<adrW), exponent.toInt)
        val y = math.round(scalb(atan(x), -exponent-1) * (1L<<fracW))
        assert(y < (1L << fracW))
        y.U((fracW+1).W)
      } ) )
    } ) )

    assert(maxCbit(0) == fracW)

    val c0 = tbl(exAdr)(adr)
    io.cs.cs(0) := ShiftRegister(c0, nStage)

  } else {

    val cbit = ATan2Sim.atanTableGeneration( order, adrW, manW, fracW )
      .map( t => {t.getCBitWidth(/*sign mode = */0)} )
      .reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

    val tableIs = VecInit(
      ATan2Sim.atanTableGeneration( order, adrW, manW, fracW ).map(t => {
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
    io.cs := ShiftRegister(coeffs, nStage)
  }
}

// -------------------------------------------------------------------------
//                        _        _     _                    _   _       _
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__   / |
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \  | |
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | | | |
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_| |_|
//                                               |_|
// -------------------------------------------------------------------------

class ATan2Stage1NonTableOutput(val spec: RealSpec) extends Bundle {
  val zex       = Output(UInt(spec.exW.W)) // sign of min(x,y)/max(x,y)
  val maxXYMan0 = Output(Bool())           // max(x,y).man === 0.U
  val xySameMan = Output(Bool())
}

class ATan2Stage1OtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val yIsLarger = Input(Bool())
    val x       = Flipped(new DecomposedRealOutput(spec))
    val y       = Flipped(new DecomposedRealOutput(spec))

    val zother  = new ATan2Stage1NonTableOutput(spec)
  })

  val maxXYMan0 = Mux(io.yIsLarger, (!io.y.man.orR.asBool), (!io.x.man.orR.asBool))
//   printf("x.man = %b\n", io.x.man)
//   printf("y.man = %b\n", io.y.man)
//   printf("x < y = %b\n", io.yIsLarger)

  val xySameMan = io.x.man === io.y.man

  io.zother.maxXYMan0 := ShiftRegister(maxXYMan0, nStage)
  io.zother.xySameMan := ShiftRegister(xySameMan, nStage)

  // --------------------------------------------------------------------------
  // Here we don't need to check if 1/max(x,y) is a special value because
  // Stage1PreProcess.io.special covers all the cases.

  // 1/x = 2^(-e-1) * 2/1.m
  // ex = -(x.ex - exBias) - 1 + exBias
  //    = -x.ex + exBias - 1 + exBias
  //    = exBias * 2 - 1 - x.ex
  //
  // y/x = 2^(y.e) * y.m * 2^(-x.e-1) * 2/x.m
  //     = 2^(y.e - x.e - 1) * y.m * (2/x.m)
  //
  // y/x.ex = (y.e - exBias - x.e + exBias - 1 + exBias
  //        = (y.e - x.e - 1 + exBias)

  val xexBiased = Mux(io.yIsLarger, io.y.ex, io.x.ex) // if y>x, swap x and y
  val yexBiased = Mux(io.yIsLarger, io.x.ex, io.y.ex)
  val zex0 = Wire(UInt(exW.W))

  // Since y < x, ex is always smaller than exBias. In normal cases, exBias is
  // around a half of 2^exW-1, so we have almost a half of the space of the
  // output port, UInt(exW.W). We re-interpret it as a signed integer to keep
  // information. If we round the negative value to zero, the postprocess
  // might consider the result is a small but non-zero value, though actually
  // that is less than the minimum.
  zex0 := yexBiased - xexBiased + (exBias-1).U(exW.W)

  // exponent of min(x,y)/max(x,y). we will later correct +/- 1 by checking
  // the mantissa of min(x,y) and max(x,y)
  io.zother.zex := ShiftRegister(zex0, nStage)
}

// -----------------------------------------------------------------------------
//                        _        _     _                    _   _       ____
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__   |___ \
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \    __) |
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |  / __/
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_| |_____|
//                                               |_|
// -----------------------------------------------------------------------------

class ATan2Stage2NonTableOutput(val spec: RealSpec) extends Bundle {
  val zsgn        = Output(UInt(1.W))
  val zex         = Output(UInt(spec.exW.W))
  val zman        = Output(UInt(spec.manW.W))
  val zIsNonTable = Output(Bool())
  val correctionNeeded = Output(Bool())
}

class ATan2Stage2OtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val flags   = Input(new ATan2Flags())
    val x       = Flipped(new DecomposedRealOutput(spec))
    val zother  = new ATan2Stage2NonTableOutput(spec)
  })

  val pi         = new RealGeneric(spec, Pi)
  val halfPi     = new RealGeneric(spec, Pi * 0.5)
  val quarterPi  = new RealGeneric(spec, Pi * 0.25)
  val quarter3Pi = new RealGeneric(spec, Pi * 0.75)

  val linearThreshold = (ATan2Sim.calcLinearThreshold(manW) + exBias)
  val isLinear = io.x.ex < linearThreshold.U(exW.W)

  val zex0 = io.x.ex - (exBias - 1).U(exW.W) // later we need to correct it

  val xzero = !io.x.ex.orR

  val defaultEx  = Mux(isLinear, io.x.ex, zex0) // isLinear includes xzero.
  val defaultMan = Mux(xzero,    0.U(exW), io.x.man) // we need to re-set man if x is zero
  // non-linear mantissa is calculated by table.

  io.zother.zIsNonTable      := isLinear || (io.flags.special =/= ATan2SpecialValue.zNormal)
  io.zother.correctionNeeded := isLinear || (io.flags.special === ATan2SpecialValue.zNormal)

  io.zother.zsgn := io.flags.ysgn
  io.zother.zex  := MuxCase(defaultEx, Seq(
    (io.flags.special === ATan2SpecialValue.zNaN)        -> Fill(exW, 1.U(1.W)),
    (io.flags.special === ATan2SpecialValue.zZero)       -> 0.U(exW.W),
    (io.flags.special === ATan2SpecialValue.zPi)         -> pi.ex.U(exW.W),
    (io.flags.special === ATan2SpecialValue.zHalfPi)     -> halfPi.ex.U(exW.W),
    (io.flags.special === ATan2SpecialValue.zQuarterPi)  -> quarterPi.ex.U(exW.W),
    (io.flags.special === ATan2SpecialValue.z3QuarterPi) -> quarter3Pi.ex.U(exW.W)
    ))
  io.zother.zman := MuxCase(defaultMan, Seq(
    (io.flags.special === ATan2SpecialValue.zNaN)        -> Fill(manW, 1.U(1.W)),
    (io.flags.special === ATan2SpecialValue.zZero)       -> 0.U(manW.W),
    (io.flags.special === ATan2SpecialValue.zPi)         -> pi.man.toLong.U(manW.W),
    (io.flags.special === ATan2SpecialValue.zHalfPi)     -> halfPi.man.toLong.U(manW.W),
    (io.flags.special === ATan2SpecialValue.zQuarterPi)  -> quarterPi.man.toLong.U(manW.W),
    (io.flags.special === ATan2SpecialValue.z3QuarterPi) -> quarter3Pi.man.toLong.U(manW.W)
    ))
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

// Stage1: takes 1/max(x, y) and min(x, y), returns min(x, y) / max(x, y)

class ATan2Stage1PostProcess(
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
    val zother = Flipped(new ATan2Stage1NonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val minxy  = Flipped(new DecomposedRealOutput(spec))
    val z      = Output(UInt(spec.W.W))
  })

  val zsgn      = 0.U(1.W)
  val zex0      = io.zother.zex
  val maxXYMan0 = io.zother.maxXYMan0
  val xySameMan = io.zother.xySameMan

//   printf("cir: zex0 = %b\n", zex0)
//   printf("cir: zres = %b\n", io.zres)

  val denomW1 = Cat(1.U(1.W), Mux(maxXYMan0, 0.U, io.zres))
  val numerW1 = Cat(1.U(1.W), io.minxy.man)

//   printf("cir: denomW1 = %b\n", denomW1)
//   printf("cir: numerW1 = %b\n", numerW1)

  val zProd     = denomW1 * numerW1
  val bp        = fracW + manW
  val roundBits = fracW + manW - manW
  val zProdMoreThan2 = zProd((fracW+1)+(manW+1)-1)
  val zProdSticky    = zProd(roundBits-2, 0).orR | (zProdMoreThan2 & zProd(roundBits-1))
  val zProdRound     = Mux(zProdMoreThan2, zProd(roundBits),       zProd(roundBits-1))
  val zProdShifted   = Mux(zProdMoreThan2, zProd(bp, roundBits+1), zProd(bp-1, roundBits))
  assert(zProdShifted.getWidth == manW)
  val zProdLSB       = zProdShifted(0)
  val zProdInc       = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
    zProdLSB, zProdRound, zProdSticky)
  val zProdRounded   = zProdShifted +& zProdInc
  assert(zProdRounded.getWidth == manW+1)
  val zProdMoreThan2AfterRound = zProdRounded(manW)

  val zex = Wire(UInt(exW.W))

  // the result of OtherPath might cause underflow.
  val zex0Inc = zex0 + zProdMoreThan2 + zProdMoreThan2AfterRound + xySameMan.asUInt
  val canUnderflow = (spec.exMin - spec.exMax + exBias < 0)
  if (canUnderflow) {
    zex := Mux(zex0Inc(exW-1), 0.U, zex0Inc)
  } else {
    zex := zex0Inc
  }

  val zman = Mux(~zex.orR || xySameMan, 0.U(manW.W), zProdRounded(manW-1, 0))

  val z0 = Cat(zsgn, zex, zman)

  io.z := ShiftRegister(z0, nStage)
}

// -------------------------------------------------------------------------
//                  _                                       ____
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___  |___ \
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|   __) |
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \  / __/
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/ |_____|
// |_|                 |_|
// -------------------------------------------------------------------------
//
// Stage2: takes status flags, returns corrected result
//

class ATan2Stage1PostProcess(
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
    val zother = Flipped(new ATan2Stage1NonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val flags  = Input(new ATan2Flags())
    val z      = Output(UInt(spec.W.W))
  })

// TODO
}
