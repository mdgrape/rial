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
  val ysgn    = UInt(1.W)
  val special = UInt(ATan2SpecialValue.W.W)
}

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

// TODO: post preprocess

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

// TODO

// -------------------------------------------------------------------------
//                        _        _     _                    _   _
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_|
//                                               |_|
// -------------------------------------------------------------------------

class ATan2Stage1NonTableOutput(val spec: RealSpec) extends Bundle {
  val zex       = Output(UInt(spec.exW.W)) // sign of min(x,y)/max(x,y)
  val maxXYMan0 = Output(Bool())           // max(x,y).man === 0.U
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

    val special = Output(UInt(ATan2SpecialValue.W.W))
    val zother  = new ATan2Stage1NonTableOutput(spec)
  })

  val xpos = io.x.sgn === 0.U
  val xneg = io.x.sgn === 1.U
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
  val maxXYMan0 = Mux(io.yIsLarger, (!io.y.man.orR.asBool), (!io.x.man.orR.asBool))
//   printf("x.man = %b\n", io.x.man)
//   printf("y.man = %b\n", io.y.man)
//   printf("x < y = %b\n", io.yIsLarger)

  io.special := ShiftRegister(special, nStage)
  io.zother.maxXYMan0 := ShiftRegister(maxXYMan0, nStage)

  // --------------------------------------------------------------------------
  // When 1/max(x,y) becomes a special value, the atan2 will also become a
  // special value like:
  // - 1/max(x,y) == nan -> atan2(y,x) == nan
  // - 1/max(x,y) == inf -> atan2(y,x) == nan      (because x == y == 0)
  // - 1/max(x,y) == zero-> atan2(y,x) == +/- pi/2 (because max(x,y) == inf)
  // So here we don't need to check if 1/max(x,y) is a special value because
  // io.special covers all the cases.

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
  val yOverXEx0 = Wire(UInt(exW.W))
  // since y < x, yOverXEx0 is always negative, so we don't need to check overflow
  val canUnderflow = (spec.exMin - spec.exMax + exBias < 0)
  if (canUnderflow) {
    val exDiff = yexBiased +& (exBias-1).U - xexBiased
    yOverXEx0 := Mux(exDiff(exW), 0.U(exW.W), exDiff(exW-1, 0))
  } else {
    yOverXEx0 := yexBiased - xexBiased + (exBias-1).U
  }

  // exponent of min(x,y)/max(x,y). we will later correct +/- 1 by checking
  // the mantissa of min(x,y) and max(x,y)
  io.zother.zex := ShiftRegister(yOverXEx0, nStage)
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

  val zex = zex0 + zProdMoreThan2 + zProdMoreThan2AfterRound
  val zman = Mux(~zex.orR, 0.U(manW.W), zProdRounded(manW-1, 0))

  val z0 = Cat(zsgn, zex, zman)

  io.z := ShiftRegister(z0, nStage)
}

// Post: takes status flags, returns corrected result
