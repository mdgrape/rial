//% @file scaleMixtureGaussian.scala
//
// scale mixture gaussian function
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

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

// -------------------------------------------------------------------------
//
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

class ScaleMixtureGaussianPreProcess(
  val sgmA:      Double,
  val sgmB:      Double,
  val spec:      RealSpec, // Input / Output floating spec
  val polySpec:  PolynomialSpec,
  val stage:     PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val dxW    = polySpec.dxW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val en  = Input (UInt(1.W))
    val x   = Flipped(new DecomposedRealOutput(spec))
    val adr = Output(UInt(adrW.W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
    val useTable = Output(Bool())
  })

  val xTableMaxDigit = ScaleMixtureGaussianSim.tableDomainDigit(manW, sgmA, sgmB)
  assert(xTableMaxDigit <= spec.exMax)

  val xTableShift = (xTableMaxDigit + exBias).U(exW.W) - io.x.ex
  io.useTable := enable(io.en, (xTableMaxDigit + exBias).U(exW.W) > io.x.ex)

  val xTableMan = Cat(1.U(1.W), io.x.man) >> xTableShift

//   printf("cir: xTableMan = %b\n", xTableMan)

  val adr  = enable(io.en, xTableMan(manW-1, dxW))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx   = enable(io.en, Cat(~xTableMan(manW-adrW-1), xTableMan(manW-adrW-2,0)))
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
//
//
class ScaleMixtureGaussianTableCoeff(
  val sgmA:      Double,
  val sgmB:      Double,
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
    val adr = Input(UInt(adrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  val adr = io.adr

  if(order == 0) {

    val tableI = ScaleMixtureGaussianSim.tableGeneration( order, adrW, manW, fracW, sgmA, sgmB )
    val cbit   = tableI.cbit

    // sign mode 0: always include sign
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode = */1)
    val coeff = getSlices(coeffTable(adr), coeffWidth)

    io.cs.cs(0) := enable(io.en, coeff(0))

  } else {

    val tableI = ScaleMixtureGaussianSim.tableGeneration( order, adrW, manW, fracW, sgmA, sgmB )
    val cbit   = tableI.cbit

    // sign mode 0: always include sign
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode = */0)
    val coeff = getSlices(coeffTable(adr), cbit)

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

object ScaleMixtureGaussianTableCoeff {
  def getCBits(
    sgmA:     Double,
    sgmB:     Double,
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
      return ScaleMixtureGaussianSim.tableGeneration(
        order, adrW, spec.manW, fracW, sgmA, sgmB
      ).getCBitWidth(/*sign mode = */0)
    }
  }

  def getCalcW(
    sgmA:     Double,
    sgmB:     Double,
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
      return ScaleMixtureGaussianSim.tableGeneration(
        order, adrW, spec.manW, fracW, sgmA, sgmB
      ).calcWidth
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
//
// round PreProcMultiplier result (x / sgmA^2).
//

// PreProcMultiplier argument (1/sgmA^2)
object ScaleMixtureGaussianPreMulArgs {
  def lhsW(spec: RealSpec): Int = {
    1 + spec.manW
  }

  def lhs(sgmA: Double, sgmB: Double, spec: RealSpec, w: Int): UInt = {
    val diffW = w - lhsW(spec)
    assert(diffW >= 0, f"PreProcMultiplier lhsW(${w}) should be wider or " +
      f"equal to lhsW(${lhsW(spec)})")

    val coef = new RealGeneric(spec, 1.0 / (sgmA * sgmA))
    if(diffW != 0) {
      Cat(coef.manW1.toBigInt.U(lhsW(spec).W), 0.U(diffW.W))
    } else {
      coef.manW1.toBigInt.U(lhsW(spec).W)
    }
  }

  def prodW(spec: RealSpec): Int = {
    (1 + spec.manW) + lhsW(spec)
  }
  def prod(spec: RealSpec, out: UInt): UInt = {
    val outW = out.getWidth
    out(outW-1, outW - prodW(spec))
  }
}

class ScaleMixtureGaussianOtherPath(
  val sgmA:     Double,
  val sgmB:     Double,
  val spec:     RealSpec, // Input / Output floating spec
  val polySpec: PolynomialSpec,
  val stage:    PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val fracW     = polySpec.fracW
  val extraBits = polySpec.extraBits

  val io = IO(new Bundle {
    val xex        = Input(UInt(exW.W))
    val xsgmA2Prod = Input(UInt(ScaleMixtureGaussianPreMulArgs.prodW(spec).W))
    val xsgmA2Man  = Output(UInt(manW.W))
    val xsgmA2Ex   = Output(UInt(exW.W))
  })

  val rsgmA2 = new RealGeneric(spec, 1.0 / (sgmA * sgmA))

  val prodW = ScaleMixtureGaussianPreMulArgs.prodW(spec)
  val prod  = io.xsgmA2Prod

  val moreThan2 = prod(prodW-1)

  val roundBits = manW
  val sticky    = prod(roundBits-2, 0).orR | (moreThan2 & prod(roundBits-1))
  val round     = Mux(moreThan2, prod(roundBits), prod(roundBits-1))
  val shifted   = Mux(moreThan2, prod(prodW-2, (manW+1)), prod(prodW-3, (manW+1)-1))
  val lsb       = shifted(0)
  val roundInc  = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven, lsb, round, sticky)
  val rounded   = shifted +& roundInc
  val moreThan2AfterRound = rounded(manW)
  val exInc     = moreThan2 | moreThan2AfterRound

  val man = rounded(manW-1, 0)

  assert(rsgmA2.toDouble > 1.0, "sgmA < 1.0, so 1 / sgmA^2 > 1.0")
  val ex0 = io.xex +& (rsgmA2.ex - spec.exBias).U +& exInc
  val ex  = Mux(ex0 > (spec.exMax + exBias).U, (spec.exMax + exBias + 1).U, ex0)

//   printf("cir: xsgmA2.ex = %d, man = %b\n", ex, man)

  io.xsgmA2Ex  := ShiftRegister(ex,  nStage)
  io.xsgmA2Man := ShiftRegister(man, nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------


// at the postprocess, we multiply (x/sgmA^2) and (polynomial+1).
class ScaleMixtureGaussianPostMulArgs(
  val sgmA:     Double,
  val sgmB:     Double,
  val spec:     RealSpec, // Input / Output floating spec
  val polySpec: PolynomialSpec,
  val stage:    PipelineStageConfig,
) extends Module {

  assert(stage.total == 0)

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val fracW  = polySpec.fracW

  val io = IO(new Bundle {
    val en     = Input(UInt(1.W))
    val useTable = Input(Bool())
    val xman   = Input(UInt(manW.W))  // result of the non-table path
    val zres   = Input(UInt(fracW.W)) // result of the table
    val lhs    = Output(UInt((1+fracW).W))
    val rhs    = Output(UInt((1+manW).W))
    val z0ex   = Output(UInt(exW.W))
  })

//   printf("cir: useTable = %b\n", io.useTable)
  val zTable = Mux(io.useTable, io.zres, 0.U(fracW.W))

  val zTableScaleDigit = ScaleMixtureGaussianSim.tableMaxBitDigit(sgmA, sgmB)
  val zTableScaled = Cat(zTable, 0.U(zTableScaleDigit.W))

//   printf("cir: zTableScaled = %b\n", zTableScaled)
//   printf("cir: 1            = %b\n", Cat(1.U(1.W), 0.U(fracW.W)))

  val z0 = zTableScaled + Cat(1.U(1.W), 0.U(fracW.W))
  val z0W = (fracW + zTableScaleDigit).U - PriorityEncoder(Reverse(z0))
  val z0Shift = z0W - (1+fracW).U
//   printf("cir: zman0      = %b\n", z0)
//   printf("cir: zman0W     = %d\n", z0W)
//   printf("cir: zman0Shift = %d\n", z0Shift)

  val z0Rounded0 = Cat(z0, 0.U(1.W)) >> z0Shift
  val z0Rounded  = z0Rounded0(z0Rounded0.getWidth-1, 1) +& z0Rounded0(0)
  val z0MoreThan2AfterRound = z0Rounded.head(1)

  val z0man = Mux(z0MoreThan2AfterRound === 1.U, (z0Rounded >> 1), z0Rounded)

  val z0ex0 = z0W +& (exBias - 1 - fracW).U(exW.W) +& z0MoreThan2AfterRound
//   printf("cir: z0ex0 = %d\n", z0ex0)
  val z0ex  = Mux(z0ex0.head(2).orR, (spec.exMax + 1 + exBias).U(exW.W), z0ex0(exW-1, 0))

//   printf("cir: z0man = %b\n", z0man)
//   printf("cir: z0ex = %d\n", z0ex)

  io.lhs := enable(io.en, z0man)                  // = polynomial + 1.0
  io.rhs := enable(io.en, Cat(1.U(1.W), io.xman)) // = x / sgmA^2

  io.z0ex := z0ex
}

class ScaleMixtureGaussianPostProcess(
  val sgmA:     Double,
  val sgmB:     Double,
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
    val en       = Input(UInt(1.W))
    val xsgn     = Input(UInt(1.W))
    val xsgmA2Ex = Input(UInt(exW.W))
    val z0ex     = Input(UInt(exW.W))
    // postproc mult
    val zman0    = Input(UInt(manW.W))
    val zexInc   = Input(UInt(1.W))
    // result
    val z        = Output(UInt(spec.W.W))
  })

  val zsgn = ~io.xsgn

  val zex00 = io.xsgmA2Ex +& io.z0ex +& io.zexInc
  val zex0  = Mux(zex00 < exBias.U, 0.U, zex00 - exBias.U)
  val zex = Mux(zex0 > (spec.exMax + exBias).U,
    (spec.exMax + 1 + exBias).U(exW.W), zex0(exW-1, 0))
//   printf("cir: zex0 = %d, zex0NoBias = %d\n", zex0, zex0 - exBias.U)

  // TODO isnan
  val zman = Mux(zex === 0.U || zex === (spec.exMax + 1 + exBias).U, 0.U, io.zman0)

  val z0 = Cat(zsgn, zex, zman)
  val z = enable(io.en, z0)

  io.z := ShiftRegister(z, nStage)
}
