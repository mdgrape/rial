//% @file exponential.scala
//
// pow2 and exp function
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

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class Pow2PreProcess(
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
    val x    = Input (UInt(spec.W.W))
    val xint = Output(UInt(exW.W)) // integer part of x
    val adr  = Output(UInt((1+adrW).W)) // we use LSB of x.ex
    val dx   = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  val padding = if (adrW>=manW) {
    if (nOrder != 0) {
      println("ERROR: table address width >= mantissa width, but polynomial order is not zero.")
      sys.exit(1)
      0
    } else {
      adrW-manW
    }
  } else {
    extraBits
  }

  val bp       = manW + padding    // for rounding
  val exValidW = log2Up(bp+expW-1) // valid exponent bits

  val (sgn, ex, man) = FloatChiselUtil.decompose(spec, io.x)

  // Float to Integer, keep fractional width bp
  // max exponent is expW-2
  val lsbPadding = expW-2+padding
  val manWith1   = 1.U(1.W) ## man ## 0.U(lsbPadding.W)
  val man1W      = manW + padding + expW - 1
  // manWith1 width: manW+extraBits+expW-1
  //   binary point: manW+extraBits

  // right shift amount :
  //   0  if ex==exBias+expW-2
  //   shift = exBias+expW-2-ex
  //   0 if shift >= bp+expW-1
  //     => ex <= exBias-bp-1
  val shift_base  = (exBias+expW-2) & maskI(exValidW)
  val shiftToZero = ex < (exBias-bp).U(expW.W)
  val shift       = shift_base.U - ex(exValidW-1,0)
  val xshift      = manWith1 >> shift
  // bp includes extraBits bit for rounding / extension
  val xu = 0.U(1.W) ## Mux(shiftToZero, 0.U(man1W.W), xshift)
  val xi = Mux(sgn, -xu, xu) // bp+expW
  val xiInt  = xi(bp+expW-1,bp) // Integral part - for exponent
  val xiFrac = xi(bp-1, 0)
  //printf("ex=%d shiftToZero=%b shift=%d xi=%x\n", ex, shiftToZero, shift, xi)

  val adr0 = xiFrac(bp-1, bp-1-adrW+1)
  io.adr := ShiftRegister(adr0, nStage)

  if(order != 0) {
    val dx0  = Cat(~xiFrac(bp-1-adrW), xiFrac(bp-1-adrW-1, 0))
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

class Pow2TableCoeff(
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

  if(order == 0) {
    val tbl = VecInit( (0L to (1L<<adrW)-1L).map(
      n => {
        val x = n.toDouble/(1L<<adrW)
        val y = round((pow(2.0,x)-1.0)*(1L<<manW))
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
    val tableD = new FuncTableDouble( x => pow(2.0,x)-1.0, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)

    val tableI = new FuncTableInt(tableD, fracW)
    val cbit   = tableI.cbit

    // both 1st and 2nd derivative of 2^x is larger than 0
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/1)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      val ci  = coeff(i)
      coeffs.cs(i) := Cat(Fill(diffWidth, 1.U(1.W)), ci) // zero extension
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

class Pow2NonTableOutput(val spec: RealSpec) extends Bundle {
  //  zsgn  = 0, always
  val zex   = Output(UInt(spec.exW.W))
  val znan  = Output(Bool())
  val zIsNonTable = Output(Bool())
}

// No pathway other than table interpolation. just calculate ex and sgn.
class Pow2OtherPath(
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
    val xint   = Input(UInt(exW.W))
    val zother = new Pow2NonTableOutput(spec)
  })

  val exOvfUdf = io.x.ex >= (exW-1+exBias).U
  val exOvf    = exOvfUdf && !x.io.sgn
  val exUdf    = exOvfUdf &&  x.io.sgn
  val znan     = if (spec.disableNaN) {false.B} else {// iff x is nan, z is nan.
    (io.x.ex === Fill(exW, 1.U(1.W))) && ( io.x.man =/= 0.U )
  }
  val zinf  = exOvf
  val zzero = exUdf

  val zExNegMax = (xint(exW-1)===1.U) && (xint(exW-2,1)===0.U)
  val zEx0 = exBias.U + xiInt
  val zEx = Mux( exOvf || nan , Fill(exW, 1.U(1.W)),
            Mux( exUdf || zExNegMax, 0.U, zEx0))

  val zIsNonTable = znan || zinf || zzero

  io.zother.zex  := ShiftRegister(zEx,  nStage)
  io.zother.znan := ShiftRegister(znan, nStage)
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

class Pow2PostProcess(
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
    // ex and some flags
    val zother = Flipped(new Pow2NonTableOutput(spec))
    // table interpolation results
    val zres   = Input(UInt(fracW.W))
    // output
    val z      = Output(UInt(spec.W.W))
  })

  val zman0  = dropLSB(extraBits, zres) + zres(extraBits-1)
  val polynomialOvf = zman0.head(2) === 1.U // >=1
  val polynomialUdf = zman0.head(1) === 1.U // Negative
  val zeroFlush     = polynomialUdf || io.zother.zIsNonTable
  val zman_checkovf = Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman0(manW-1,0))

  val zMan = Mux(zeroFlush, Cat(io.zother.znan, 0.U((manW-1).W)), zman_checkovf)
  val z0   = Cat(0.U(1.W), io.zother.zex, zMan)

  io.z   := ShiftRegister(z0, nStage)
}
