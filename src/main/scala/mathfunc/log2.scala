//% @file log2.scala
//
// log2 function
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
import rial.arith.RealGeneric
import rial.arith.FloatChiselUtil

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

//
// log2(2^ex * 1.man) = log2(2^ex) + log2(1.man)
//                    = ex + log2(1.man)
//

class Log2PreProcess(
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
    val x   = Input (UInt(spec.W.W))
    val adr = Output(UInt(adrW.W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  val adr0 = io.x(manW-1, dxW)
  io.adr := ShiftRegister(adr0, nStage)

  if(order != 0) {
    val dx0  = Cat(~io.x(dxW-1), io.x(dxW-2, 0))
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

class Log2TableCoeff(
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
    val adr = Input  (UInt((1+adrW).W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  val log2 = (a:Double) => {log(a) / log(2.0)}

  if(order == 0) {
    val tbl = VecInit( (0L to (1L<<adrW)-1L).map(
      n => {
        val x = n.toDouble/(1L<<adrW)
        val y = round( log2(1.0+x) * (1L<<manW) )
        if (y>=(1L<<manW)) {
          println("WARNING: mantissa reaches to 2")
          maskL(manW).U(manW.W)
        } else {
          y.U(manW.W)
        }
      }
    ) )

    assert(maxCbit(0) == fracW)

    val c0 = tbl(io.adr(adrW, 0))            // here we use LSB of ex
    io.cs.cs(0) := ShiftRegister(c0, nStage) // width should be manW + extraBits

  } else {

    val tableD = new FuncTableDouble( x => log2(1.0 + x), order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    val tableI = new FuncTableInt(tableD, fracW)
    val cbit   = tableI.cbit

    // TODO: add specific table for the case of ex == 0 and ex == -1

    // both 1st and 2nd derivative of 2^x is larger than 0
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      val ci  = coeff(i)
      val msb = ci(cbit(i)-1)
      coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
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

class Log2NonTableOutput(val spec: RealSpec) extends Bundle {

  // 8 -> 3
  val log2Ceil = (a:Double) => {ceil(log(a) / log(2.0)).toInt}

  val zsgn  = Output(UInt(1.W))
  val zint  = Output(UInt(spec.exW.W))
  val znan  = Output(Bool())
  val zinf  = Output(Bool())
  val zzero = Output(Bool())
  val zIsNonTable = Output(Bool())
}

class Log2OtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val padding = exBias

  val io = IO(new Bundle {
    val x      = Flipped(new DecomposedRealOutput(spec))
    val zother = new Log2NonTableOutput(spec)
  })

  // --------------------------------------------------------------------------
  // calc zex
  //
  // z = log2(2^ex * 1.man)
  //   = ex + log2(1.man)
  //
  // 0 = log2(1) <= log2(1.man) < log2(2) = 1
  //
  // ex           <= ex + log2(1.man)   < ex + 1
  // ex           <= z                  < ex + 1
  // ex           <= 2^(zex) * 1.zman   < ex + 1
  // log2(ex)     <= zex + log2(1.zman) < log2(ex + 1)
  // log2(ex) - 1 <  zex                < log2(ex + 1)
  //
  val xexPos = io.x.ex - exBias.U // ==  x
  val xexNeg = (exBias-1).U - io.x.ex // == -x == abs(x)
  val xexPosW = PriorityEncoder(xexPos)
  val xexNegW = PriorityEncoder(xexNeg)

  // integer part of z (fractional part is calculated in polynomial module)
  val zint = Mux(io.x.sgn === 0.U, xexPos, xexNeg)

  // TODO: If xexNobias === 0 or -1, we need to count zeros in MSB.
  //       And, to determine `shift` in postprocess, we need to pass it to post.

  io.zother.zint := ShiftRegister(zint, nStage)

  // --------------------------------------------------------------------------
  // check special value
  //
  // log2(nan) -> nan
  // log2(inf) -> inf
  // log2(0)   -> -inf
  // log2(1)   -> 0
  // log2(-|x|) -> nan

  val xmanAllZero = !io.x.man.orR
  val znan  = io.x.nan || io.x.sgn === 1.U
  val zinf  = io.x.inf || io.x.zero
  val zzero = xmanAllZero && io.x.ex === exBias.U

  val zIsNonTable = znan || zinf || zzero

  val zsgn = io.x.ex < exBias.U

  io.zother.zsgn        := ShiftRegister(zsgn,  nStage)
  io.zother.znan        := ShiftRegister(znan,  nStage)
  io.zother.zinf        := ShiftRegister(zinf,  nStage)
  io.zother.zzero       := ShiftRegister(zzero, nStage)
  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class Log2PostProcess(
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

  val padding = extraBits

  val io = IO(new Bundle {
    val x      = Flipped(new DecomposedRealOutput(spec))
    val zother = Flipped(new Log2NonTableOutput(spec))
    val zres   = Input(UInt(fracW.W))
    val z      = Output(UInt(spec.W.W))
  })

  // --------------------------------------------------------------------------
  // z is large enough (xex != 0, -1)

  val zint   = io.zother.zint
  val zfrac0 = io.zres
  val zfrac  = Mux(io.x.ex >= exBias.U, zfrac0, (~zfrac0) + 1.U)

//   printf("c: zfrac0 = %b\n", zfrac0)
//   printf("c: zfrac  = %b\n", zfrac )

  // z = log2(2^xex * 1.xman)
  //   = xex + log2(1.xman)
  // if xex < 0, we need to subtract log2(1.xman) from |xex| and set sgn to 1

  val log2 = (a:Double) => {log(a) / log(2.0)}

  val zFull    = Cat(zint, zfrac)
  val zShiftW  = Mux(zFull(exW+fracW-1) === 1.U, 0.U, PriorityEncoder(Reverse(zFull)))
  val zShifted = zFull << zShiftW
//   printf("c: zfull    = %b\n", zfrac )
//   printf("c: zshiftW  = %b\n", zfrac )
//   printf("c: zshifted = %b\n", zShifted )

  assert(zShifted(fracW+exW-1) === 1.U)

//   printf("c: zman0 = %b + %b\n", zShifted(fracW+exW-2, fracW+exW-1-manW), zShifted(fracW+exW-1-manW-1))
  val zman0    = zShifted(fracW+exW-2, fracW+exW-1-manW) + zShifted(fracW+exW-1-manW-1)
  val zex0     = (exBias + exW - 1).U - zShiftW

  // --------------------------------------------------------------------------

  val zSgn = io.zother.zsgn
  val zEx  = Mux(io.zother.znan || io.zother.zinf, Fill(exW, 1.U(1.W)),
             Mux(io.zother.zzero, 0.U(exW.W), zex0))
  val zMan = Mux(io.zother.zIsNonTable, Cat(io.zother.znan, 0.U((manW-1).W)),
                 zman0(manW-1, 0))
  val z0   = Cat(zSgn, zEx, zMan)

  assert(zEx .getWidth == exW)
  assert(zMan.getWidth == manW)
  assert(z0  .getWidth == spec.W)

  io.z   := ShiftRegister(z0, nStage)
}
