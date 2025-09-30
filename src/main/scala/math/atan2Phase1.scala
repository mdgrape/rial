//% @file atan2.scala
//
// ATan2 function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import java.lang.Math.scalb
import scala.language.reflectiveCalls
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.math.Real
import spire.implicits._

import chisel3._
import chisel3.util._
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith._

// ATan2 Phase1 calculates min(x,y)/max(x,y).
//       Phase2 calculates atan(min(x,y)/max(x,y)) +/- constant.

object ATan2Status {
  val width = 2
  val xIsPosIsLarger  = 0.U(width.W) // ysgn *   atan(|y/x|)      .. x>0, |x|>|y|
  val xIsNegIsLarger  = 1.U(width.W) // ysgn * (-atan(|y/x|)+pi)  .. x<0, |x|>|y|
  val xIsPosIsSmaller = 2.U(width.W) // ysgn * (pi/2-atan(|x/y|)) .. x>0, |x|<|y|
  val xIsNegIsSmaller = 3.U(width.W) // ysgn * (pi/2+atan(|x/y|)) .. x<0, |x|<|y|
}

// encoded by NaN-boxing. sign of NaN encodes ysgn.
//
// we assume that, if x < 1.0, the 2 bit at the MSB side of the exponent is not
// used (the MSB of the exBias is zero, like 0b0111_1111).
object ATan2SpecialValue {
  val width = 3
  val zNormal     = 0.U(width.W) // normal
  val zZero       = 1.U(width.W) // |x| >> |y|, x > 0 : atan2(y,x) = zero
  val zPi         = 2.U(width.W) // |x| >> |y|, x < 0 : atan2(y,x) = pi
  val zHalfPi     = 3.U(width.W) // |y| >> |x|        : atan2(y,x) = (+/-) pi/2
  val zQuarterPi  = 4.U(width.W) // |x| == |y|, x > 0 : atan2(y,x) = (+/-) pi/4
  val z3QuarterPi = 5.U(width.W) // |x| == |y|, x < 0 : atan2(y,x) = (+/-) 3pi/4
  val zNaN        = 7.U(width.W) // x or y is nan.
}
class ATan2Flags extends Bundle {
  val isSpecial = Bool()
  val special   = UInt(ATan2SpecialValue.width.W)
  val status    = UInt(ATan2Status.width.W)
  val ysgn      = UInt(1.W)
}

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------
//
// Phase1 preprocess only checks the special cases.
// Phase1 calculates |min(x,y)|/|max(x,y)|, so ReciprocalPreProcess is re-used.
//
class ATan2Phase1PreProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig
) extends Module {

  assert(spec.exBias == maskI(spec.exW-1), """
    |ATan2Phase1 calculates |min(x,y)| / |max(x,y)| and encodes information
    |like xsgn, ysgn, |x| >? |y| into the sign, msb of exponent, lsb of mantissa.
    |It depends on |min(x,y)| / |max(x,y)| <= 1 and in case of == 1 we use
    |special value flags, so we don't use the msb of the exponent.
    |It means that, if x:float < 1, x.ex.msb must be zero.""".stripMargin)

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val dxW    = polySpec.dxW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val en        = Input(UInt(1.W))
    val x         = Flipped(new DecomposedRealOutput(spec))
    val y         = Flipped(new DecomposedRealOutput(spec))
    val yIsLarger = Input(Bool())

    val flags     = Output(new ATan2Flags)
  })

  val minex     = Mux(io.yIsLarger, io.x.ex, io.y.ex)
  val maxex     = Mux(io.yIsLarger, io.y.ex, io.x.ex)
  val diffexDec = Mux(io.yIsLarger, io.y.man > io.x.man, io.x.man > io.y.man)
  val zeroed    = ( minex +& exBias.U(exW.W) ) <= ( maxex +& diffexDec.asUInt )

  // |min(x,y)| / |max(x,y)| = 0 means atan2(y,x) = (n/2)pi, n=0,1,2,3
  // case |y| << |x| && 0 < x : z = 0
  // case |y| << |x| && x < 0 : z = pi
  // case |x| << |y| && 0 < y : z = pi/2
  // case |x| << |y| && y < 0 : z = 3pi/2

  val tooLargeX = (zeroed && !io.yIsLarger) || ( io.x.inf && !io.y.inf) || (!io.x.zero &&  io.y.zero)
  val tooLargeY = (zeroed &&  io.yIsLarger) || (!io.x.inf &&  io.y.inf) || ( io.x.zero && !io.y.zero)

  val samexy    = io.x.ex === io.y.ex && io.x.man === io.y.man && !io.x.nan && !io.y.nan

  val xpos      = io.x.sgn === 0.U
  val xneg      = io.x.sgn === 1.U
  val zzero     = tooLargeX && xpos
  val zpi       = tooLargeX && xneg
  val zhalfpi   = tooLargeY
  val z1piover4 = samexy && xpos
  val z3piover4 = samexy && xneg
  val znan      = (io.x.nan || io.y.nan) || (io.x.zero && io.y.zero)

  val flags = WireDefault(0.U.asTypeOf(new ATan2Flags))
  flags.isSpecial := znan || zzero || zpi || zhalfpi || z1piover4 || z3piover4
  flags.special   := MuxCase(ATan2SpecialValue.zNormal, Seq(
    zzero     -> ATan2SpecialValue.zZero,
    zpi       -> ATan2SpecialValue.zPi,
    zhalfpi   -> ATan2SpecialValue.zHalfPi,
    z1piover4 -> ATan2SpecialValue.zQuarterPi,
    z3piover4 -> ATan2SpecialValue.z3QuarterPi,
    znan      -> ATan2SpecialValue.zNaN
  ))
  flags.status := Cat(io.yIsLarger, io.x.sgn)
  flags.ysgn   := io.y.sgn

  io.flags := ShiftRegister(enableIf(io.en, flags), nStage)
}

// -------------------------------------------------------------------------
//                        _        _     _                    _   _
//  _ __   ___  _ __     | |_ __ _| |__ | | ___   _ __   __ _| |_| |__
// | '_ \ / _ \| '_ \ ___| __/ _` | '_ \| |/ _ \ | '_ \ / _` | __| '_ \
// | | | | (_) | | | |___| || (_| | |_) | |  __/ | |_) | (_| | |_| | | |
// |_| |_|\___/|_| |_|    \__\__,_|_.__/|_|\___| | .__/ \__,_|\__|_| |_|
//                                               |_|
// -------------------------------------------------------------------------

class ATan2Phase1NonTableOutput(val spec: RealSpec) extends Bundle {
  val zex       = Output(UInt(spec.exW.W)) // sign of min(x,y)/max(x,y)
  val maxXYMan0 = Output(Bool())           // max(x,y).man === 0.U
  val xySameMan = Output(Bool())
}

class ATan2Phase1OtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val x         = Flipped(new DecomposedRealOutput(spec))
    val y         = Flipped(new DecomposedRealOutput(spec))
    val yIsLarger = Input(Bool())

    val zother = new ATan2Phase1NonTableOutput(spec)
  })

  val maxXYMan0 = Mux(io.yIsLarger, (!io.y.man.orR.asBool), (!io.x.man.orR.asBool))
//   printf("x.man = %b\n", io.x.man)
//   printf("y.man = %b\n", io.y.man)
//   printf("x < y = %b\n", io.yIsLarger)

  val xySameMan = io.x.man === io.y.man

  io.zother.maxXYMan0 := ShiftRegister(maxXYMan0, nStage)
  io.zother.xySameMan := ShiftRegister(xySameMan, nStage)

  val exDec  = Mux(Mux(io.yIsLarger, io.x.man < io.y.man, io.y.man < io.x.man),
                   1.U(1.W), 0.U(1.W))
  val maxEx  = Mux(io.yIsLarger, io.y.ex, io.x.ex)
  val minEx  = Mux(io.yIsLarger, io.x.ex, io.y.ex)
  val zeroed = minEx +& exBias.U(exW.W) <= maxEx + exDec.asUInt
  val zex0  = Mux(io.x.inf && io.y.inf, exBias.U(exW.W),
              Mux((io.x.inf && !io.y.inf) || (!io.x.inf && io.y.inf), 0.U(exW.W),
                  ((minEx +& exBias.U) - maxEx) - exDec.asUInt))
  val zex   = Mux(zeroed, 0.U, zex0)

  // exponent of min(x,y)/max(x,y). we will later correct +/- 1 by checking
  // the mantissa of min(x,y) and max(x,y)
  io.zother.zex := ShiftRegister(zex, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

// Phase1: takes 1/max(x, y) and min(x, y), returns min(x, y) / max(x, y)

class ATan2Phase1MultArgs(
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
    val en     = Input(Bool())
    val zother = Flipped(new ATan2Phase1NonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val minxy  = Flipped(new DecomposedRealOutput(spec))

    val lhs    = Output(UInt((1+fracW).W))
    val rhs    = Output(UInt((1+manW).W))
  })

  val zFrac = Mux(io.zother.maxXYMan0, 0.U, io.zres)

  io.lhs := enableIf(io.en, Cat(1.U(1.W), zFrac))
  io.rhs := enableIf(io.en, Cat(1.U(1.W), io.minxy.man))
}

class ATan2Phase1PostProcess(
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
    val en     = Input(Bool())
    val zother = Flipped(new ATan2Phase1NonTableOutput(spec))
    val zman0  = Input(UInt(manW.W))
    val minxy  = Flipped(new DecomposedRealOutput(spec))
    val flags  = Input(new ATan2Flags)

    val z      = Output(UInt(spec.W.W))
  })

  // --------------------------------------------------------------------------
  // calc result

  val zsgn      = 0.U(1.W)
  val zex       = io.zother.zex
  val maxXYMan0 = io.zother.maxXYMan0
  val xySameMan = io.zother.xySameMan
  val zman      = Mux(~zex.orR || xySameMan, 0.U(manW.W),
                  Mux(maxXYMan0, io.minxy.man, io.zman0))

  // --------------------------------------------------------------------------
  // encode state and ysgn

  val zSgnSpecial = io.flags.ysgn
  val zExSpecial  = Fill(exW, 1.U(1.W))
  val zManSpecial = io.flags.special

  assert(zex(exW-1) === 0.U)
  assert(zex(exW-2) === 0.U)

  val zSgnEncoded = io.flags.ysgn
  val zExEncoded  = Cat(io.flags.status(1), zex(exW-2, 0))
  val zManEncoded = Cat(zman(manW-1, 1), io.flags.status(0))

  val zSpecial = Cat(zSgnSpecial, zExSpecial, zManSpecial)
  val zEncoded = Cat(zSgnEncoded, zExEncoded, zManEncoded)

  val z = Mux(io.flags.isSpecial, zSpecial, zEncoded)

  io.z := ShiftRegister(enableIf(io.en, z), nStage)
}

