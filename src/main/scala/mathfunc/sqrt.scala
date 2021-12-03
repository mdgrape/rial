//% @file sqrt.scala
//
// square root function
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
import rial.arith.FloatChiselUtil

import rial.math.SqrtSim
import rial.mathfunc._

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

// sqrt(x): floating => floating
// - if x < 0, returns 0.
class SqrtPreProcess(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val manW = spec.manW

  val io = IO(new Bundle {
    val x   = Input (UInt(spec.W.W))
    val adr = Output(UInt((1+adrW).W))
    val dx  = Output(UInt((manW-adrW).W))
  })

  io.adr := io.x(manW, manW-adrW) // include LSB of x.ex
  io.dx  := Cat(~io.x(manW-adrW-1), io.x(manW-adrW-2,0))
}

// -------------------------------------------------------------------------
//  _        _     _                        __  __
// | |_ __ _| |__ | | ___    ___ ___   ___ / _|/ _|
// | __/ _` | '_ \| |/ _ \  / __/ _ \ / _ \ |_| |_
// | || (_| | |_) | |  __/ | (_| (_) |  __/  _|  _|
//  \__\__,_|_.__/|_|\___|  \___\___/ \___|_| |_|
// -------------------------------------------------------------------------

class SqrtTableCoeff(

  val spec : RealSpec,
  val nOrder: Int, val adrW : Int, val extraBits : Int,

  val maxAdrW : Int,      // max address width among all math funcs
  val maxCbit : Seq[Int], // max coeff width among all math funcs

  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val manW  = spec.manW
  val calcW = manW + extraBits
  val order = if(adrW == manW) {0} else {nOrder}

  val io = IO(new Bundle {
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
        val y = round((math.sqrt(x)-1.0) * (1L<<manW))
        if (y >= (1L<<manW)) {
          maskL(manW).U(manW.W)
        } else if (y <= 0.0) {
          0.U(manW.W)
        } else {
          y.U(manW.W)
        }
      })
    )
    assert(maxCbit(0) == calcW)

    val c0 = tbl(io.adr(adrW, 0)) // here we use LSB of ex
    io.cs.cs(0) := c0             // width should be the same, manW + extraBits

  } else {
    val tableI = SqrtSim.sqrtTableGeneration( nOrder, adrW, manW, calcW )
    val cbit   = tableI.cbit

    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    for (i <- 0 to nOrder) {
      val diffWidth = maxCbit(i) - cbit(i)
      val ci  = coeff(i)
      val msb = ci(cbit(i)-1)
      io.cs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
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

class SqrtNonTableOutput(val spec: RealSpec) extends Bundle {
  val zex   = Output(UInt(spec.exW.W))
  val zsgn  = Output(UInt(1.W))
  val znan  = Output(Bool())
  val zinf  = Output(Bool())
  val zzero = Output(Bool())
}

// No pathway other than table interpolation. just calculate ex and sgn.
class SqrtOtherPath(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val x      = Input(UInt(spec.W.W))
    val zother = new SqrtNonTableOutput(spec)
  })

  val (xsgn,  xex,  xman) = FloatChiselUtil.decompose(spec, io.x)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(spec, io.x)

  val xneg = if(spec.disableSign) {false.B} else {xsgn === 1.U(1.W)}

  val znan  = xnan
  val zinf  = xinf
  val zzero = xzero || xneg

  io.zother.znan  := znan
  io.zother.zinf  := zinf
  io.zother.zzero := zzero

  val xExNobias = xex - exBias.U
  val zex0 = xExNobias(exW-1) ## xExNobias(exW-1, 1) // arith right shift

  io.zother.zex  := Mux(zinf || znan, maskU(exW),
                    Mux(zzero,        0.U(exW.W),
                                      zex0 + exBias.U))
  io.zother.zsgn := 0.U(1.W) // always positive.
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class SqrtPostProcess(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val nStage = stage.total
  def getStage() = nStage

  val io = IO(new Bundle {
    // ex and some flags
    val zother = Flipped(new SqrtNonTableOutput(spec))
    // table interpolation results
    val zres   = Input(UInt((spec.manW + extraBits).W))
    // output
    val z      = Output(UInt(spec.W.W))
  })

  val zinf  = io.zother.zinf
  val znan  = io.zother.znan
  val zzero = io.zother.zzero

  val zsgn = io.zother.zsgn
  val zex  = io.zother.zex

  val zman0 = dropLSB(extraBits, io.zres) +& io.zres(extraBits-1)
  val polynomialOvf = zman0(manW)
  val zeroFlush     = zinf || zzero || znan
  val zman          = Mux(zeroFlush,     Cat(znan, 0.U((manW-1).W)),
                      Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman0(manW-1,0)))

  val z0 = Cat(zsgn, zex, zman)

  io.z   := ShiftRegister(z0, nStage)
}
