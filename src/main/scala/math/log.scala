//% @file log.scala
//
// log2 and ln function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._

import spire.math.SafeLong
import spire.math.Real
import spire.implicits._

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

//
// case 1: x < 0.5       or        2 <= x. calc (ex + log2(1.man)) * ln2
// case 2:     0.5 <= x < 1.               calc       log(x)/(1-x) * (1-x)
// case 3:                1 <= x < 2.      calc       log(x)/(x-1) * (x-1)
//                                                    ^^^^^^^^^^^^   ^^^^^
//                                                    table approx   constant

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------
//
// calculate table addresses.
// (ln2/1-x/x-1) will be calculated in OtherPath, in parallel with polynomial.
//
class LogPreProcess(
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
    val adr = Output(UInt((2+adrW).W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  // 0: others
  // 1: xexNobias == ex - exBias == -1 <=> ex == exBias - 1
  // 2: xexNobias == ex - exBias ==  0 <=> ex == exBias

  // used to distinguish 3 tables
  val exAdr = Mux(io.x.ex === exBias.U, 2.U, (io.x.ex === (exBias - 1).U).asUInt)

  val adr0  = Cat(exAdr, io.x.man(manW-1, dxW))
  val adr   = enable(io.en, adr0)
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx0 = Cat(~io.x.man(dxW-1), io.x.man(dxW-2, 0))
    val dx  = enable(io.en, dx0)
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

class LogTableCoeff(
  val spec     : RealSpec,
  val polySpec : PolynomialSpec,
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
) extends Module {

  val manW      = spec.manW
  val adrW      = polySpec.adrW
  val fracW     = polySpec.fracW
  val order     = polySpec.order
  val extraBits = polySpec.extraBits

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input  (UInt((2+adrW).W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {
    val exadr = io.adr(adrW+2-1, adrW)
    val adr   = io.adr(adrW-1, 0)

    val tblNormal = VecInit((0L to 1L<<adrW).map(
      n => {
        val log2 = (a:Double) => {log(a) / log(2.0)}
        val x = (n.toDouble / (1L<<adrW)) + 1.0
        val res = log2(x)
        assert(0.0 <= res && res <= 1.0, f"log2(x) = ${res}")
        val y = round(res * (1L<<fracW))
        if (y >= (1L<<fracW)) {
          maskL(fracW).U(fracW.W)
        } else if (y <= 0.0) {
          0.U(fracW.W)
        } else {
          y.U(fracW.W)
        }
      })
    )
    val tblSmallNeg = VecInit((0L to 1L<<adrW).map(
      n => {
        val x = ((n.toDouble / (1L<<adrW)) + 1.0) * 0.5 // 0~1 -> 0.5~1.0
        val res = if(x == 1.0) {0.0} else {-log(x) / (1.0 - x) - 1.0}
        assert(0.0 <= res && res <= 1.0, f"-log(x)/(1-x) - 1 = ${res}")
        val y = round(res * (1L<<fracW))

        if (y >= (1L<<fracW)) {
          maskL(fracW).U(fracW.W)
        } else if (y <= 0.0) {
          0.U(fracW.W)
        } else {
          y.U(fracW.W)
        }
      })
    )
    val tblSmallPos = VecInit((0L to 1L<<adrW).map(
      n => {
        val x = (n.toDouble / (1L<<adrW)) + 1.0 // 0~1 -> 1~2
        val res = if(n == 0) {1.0} else {log(x) / (x - 1.0)}
        assert(0.0 <= res && res <= 1.0, f"log(x)/(x-1) = ${res}")
        val y = round(res * (1L<<fracW))

        if (y >= (1L<<fracW)) {
          maskL(fracW).U(fracW.W)
        } else if (y <= 0.0) {
          0.U(fracW.W)
        } else {
          y.U(fracW.W)
        }
      })
    )

    val coeffNormal   = tblNormal  (adr)
    val coeffSmallNeg = tblSmallNeg(adr)
    val coeffSmallPos = tblSmallPos(adr)

    val coeff = Mux(exadr === 0.U, coeffNormal,
                Mux(exadr === 1.U, coeffSmallNeg, coeffSmallPos))

    io.cs.cs(0) := enable(io.en, coeff)

  } else {

    // split address
    val exadr = io.adr(adrW+2-1, adrW) // distinguish 3 tables
    val adr   = io.adr(adrW-1, 0)      // address in a table

//     printf("cir: exadr = %b\n", exadr)
//     printf("cir: adr   = %b\n", adr)

    // -----------------------------------------------------------------------
    // default table
    val tableNormalI  = LogSim.logNormalTableGeneration(spec, order, adrW, extraBits)
    val tableSmallNeg = LogSim.logSmallNegativeTableGeneration(spec, order, adrW, extraBits) // 0.5~1.0
    val tableSmallPos = LogSim.logSmallPositiveTableGeneration(spec, order, adrW, extraBits) // 1.0~2.0

    val cbitNormal   = tableNormalI .getCBitWidth(/*sign mode = */0)
    val cbitSmallNeg = tableSmallNeg.getCBitWidth(/*sign mode = */0)
    val cbitSmallPos = tableSmallPos.getCBitWidth(/*sign mode = */0)

    val tableINormal   = tableNormalI .getVectorWithWidth(cbitNormal  , /*sign mode = */0)
    val tableISmallNeg = tableSmallNeg.getVectorWithWidth(cbitSmallNeg, /*sign mode = */0)
    val tableISmallPos = tableSmallPos.getVectorWithWidth(cbitSmallPos, /*sign mode = */0)

    val coeffNormal   = getSlices(tableINormal  (adr), cbitNormal  )
    val coeffSmallNeg = getSlices(tableISmallNeg(adr), cbitSmallNeg)
    val coeffSmallPos = getSlices(tableISmallPos(adr), cbitSmallPos)

    val coeffsNormal   = Wire(new TableCoeffInput(maxCbit))
    val coeffsSmallNeg = Wire(new TableCoeffInput(maxCbit))
    val coeffsSmallPos = Wire(new TableCoeffInput(maxCbit))

    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbitNormal(i)
      assert(cbitNormal(i) <= maxCbit(i))

      if(0 < diffWidth) {
        val ci  = coeffNormal(i)
        val msb = ci(cbitNormal(i)-1)
        coeffsNormal.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffsNormal.cs(i) := coeffNormal(i)
      }
    }
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbitSmallNeg(i)
      assert(cbitSmallNeg(i) <= maxCbit(i))

      if(0 < diffWidth) {
        val ci  = coeffSmallNeg(i)
        val msb = ci(cbitSmallNeg(i)-1)
        coeffsSmallNeg.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffsSmallNeg.cs(i) := coeffSmallNeg(i)
      }
    }
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbitSmallPos(i)
      assert(cbitSmallPos(i) <= maxCbit(i))

      if(0 < diffWidth) {
        val ci  = coeffSmallPos(i)
        val msb = ci(cbitSmallPos(i)-1)
        coeffsSmallPos.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffsSmallPos.cs(i) := coeffSmallPos(i)
      }
    }

    val coeffs = Mux(exadr === 0.U, coeffsNormal.asUInt,
                 Mux(exadr === 1.U, coeffsSmallNeg.asUInt, coeffsSmallPos.asUInt)
                     ).asTypeOf(new TableCoeffInput(maxCbit))

    io.cs := enable(io.en, coeffs)
  }
}

object LogTableCoeff {
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
      val tableNormalI = LogSim.logNormalTableGeneration(
        spec, order, adrW, extraBits)
      val cbitNormal   = tableNormalI.cbit

      val tableSmallPos = LogSim.logSmallPositiveTableGeneration(
        spec, order, adrW, extraBits)
      val cbitSmallPos = tableSmallPos.cbit

      val tableSmallNeg = LogSim.logSmallNegativeTableGeneration(
        spec, order, adrW, extraBits)
      val cbitSmallNeg = tableSmallNeg.cbit

      return Seq(cbitNormal, cbitSmallPos, cbitSmallNeg).
        reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )
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
      val tableNormalI = LogSim.logNormalTableGeneration(
        spec, order, adrW, extraBits)
      val calcWidthNormal   = tableNormalI.calcWidth

      val tableSmallPos = LogSim.logSmallPositiveTableGeneration(
        spec, order, adrW, extraBits)
      val calcWidthSmallPos = tableSmallPos.calcWidth

      val tableSmallNeg = LogSim.logSmallNegativeTableGeneration(
        spec, order, adrW, extraBits)
      val calcWidthSmallNeg = tableSmallNeg.calcWidth

      return Seq(calcWidthNormal, calcWidthSmallPos, calcWidthSmallNeg).
        reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )
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
// calculate constant and zex. Also checks special cases (e.g. nan/inf, zero).
//
// case 1: x < 0.5       or        2 <= x. calc (ex + log2(1.man)) * ln2
// case 2:     0.5 <= x < 1.               calc       log(x)/(1-x) * (1-x)
// case 3:                1 <= x < 2.      calc       log(x)/(x-1) * (x-1)
//                                                    ^^^^^^^^^^^^   ^^^^^
//                                                    table approx   constant
//
class LogNonTableOutput(val spec: RealSpec, val polySpec: PolynomialSpec) extends Bundle {
  val zsgn  = Output(UInt(1.W))
  val znan  = Output(Bool())
  val zinf  = Output(Bool())
  val zzero = Output(Bool())

  val x0_5to1_0  = Output(Bool())
  val x1_0to2_0  = Output(Bool())
  val xOtherwise = Output(Bool())

  val constant  = Output(UInt((spec.manW+1).W))    // ln2, 1-x or x-1 > 0.
  val zex0      = Output(UInt(spec.exW.W))         // exponent before multiplication
  val zint      = Output(UInt(spec.exW.W))         // for default table
  val zintShift = Output(UInt(log2Up(spec.exW).W)) // for default table
}

class LogOtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val fracW     = polySpec.fracW
  val extraBits = polySpec.extraBits

  val io = IO(new Bundle {
    val x      = Flipped(new DecomposedRealOutput(spec))
    val zother = new LogNonTableOutput(spec, polySpec)
  })

  val x0_5to1_0  = io.x.ex === (exBias-1).U
  val x1_0to2_0  = io.x.ex === exBias.U
  val xOtherwise = !x0_5to1_0 && !x1_0to2_0

  io.zother.x0_5to1_0  := ShiftRegister(x0_5to1_0 , nStage)
  io.zother.x1_0to2_0  := ShiftRegister(x1_0to2_0 , nStage)
  io.zother.xOtherwise := ShiftRegister(xOtherwise, nStage)

  // --------------------------------------------------------------------------
  // constant

  // 0.5 < x < 1.0. c = 1-x.
  val c0_5to1_0 = ~io.x.man +& 1.U
  // 1.0 < x < 2.0. c = x-1.
  val c1_0to2_0 = Cat(0.U(1.W), io.x.man)

  // normalize 1-x or x-1.
  val c0 = Mux(x0_5to1_0, c0_5to1_0, c1_0to2_0)
  val c0Shift = PriorityEncoder(Reverse(c0))
  val c0Shifted = (c0 << c0Shift)(manW, 0)

  assert((BigInt(1) << manW).U <= c0Shifted || c0 === 0.U)

  // exponent of z in case of 0.5 <= x < 2.0.
  val zex0closeto1 = (exBias-1).U(exW.W) - c0Shift

  // constant in case of x < 0.5 or 2.0 <= x.
  val ln2 = Real.log(Real.two)(manW+1).toBigInt.U((manW+1).W)
  assert((BigInt(1)<<manW).U < ln2 && ln2 < (BigInt(1) << (manW+1)).U)

  // select correct one.
  val constant = Mux(xOtherwise, ln2, c0Shifted)
  io.zother.constant := ShiftRegister(constant, nStage)

  // --------------------------------------------------------------------------
  // integer part of z.
  val xexPos = io.x.ex - exBias.U
  val xexNeg = (exBias-1).U - io.x.ex // absolute value. -x.ex - 1
  val zint   = Mux(io.x.ex >= exBias.U, xexPos, xexNeg)
  io.zother.zint := ShiftRegister(zint, nStage)

  val zintShift    = PriorityEncoder(Reverse(zint))
  io.zother.zintShift := ShiftRegister(zintShift, nStage)

  // exponent of z in case of x < 0.5 || 2.0 <= x
  val zex0farfrom1 = (exBias + exW - 2).U - zintShift

  val zex0 = Mux(xOtherwise, zex0farfrom1, zex0closeto1)
  io.zother.zex0 := ShiftRegister(zex0, nStage)

  // --------------------------------------------------------------------------
  // check special value
  //
  // log(nan)  ->  nan
  // log(inf)  ->  inf
  // log(-|x|) ->  nan
  // log(0)    -> -inf
  // log(1)    ->    0

  val znan  = io.x.nan || io.x.sgn === 1.U
  val zinf  = io.x.inf || io.x.zero
  val zzero = !io.x.man.orR && io.x.ex === exBias.U
  val zsgn  = io.x.ex < exBias.U

  io.zother.znan  := ShiftRegister(znan,  nStage)
  io.zother.zinf  := ShiftRegister(zinf,  nStage)
  io.zother.zzero := ShiftRegister(zzero, nStage)
  io.zother.zsgn  := ShiftRegister(zsgn,  nStage)
}

// -------------------------------------------------------------------------
//                  _
//  _ __   ___  ___| |_ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \ / _ \/ __| __| '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | (_) \__ \ |_| |_) | | | (_) | (_|  __/\__ \__ \
// | .__/ \___/|___/\__| .__/|_|  \___/ \___\___||___/___/
// |_|                 |_|
// -------------------------------------------------------------------------

class LogPostProcess(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
  val isAlwaysLn : Boolean = false
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

  val log2 = (a:Double) => {log(a) / log(2.0)}

  val io = IO(new Bundle {
    val en     = Input(UInt(1.W))
    val isln   = if(!isAlwaysLn) {Some(Input(Bool()))} else {None}
    val zother = Flipped(new LogNonTableOutput(spec, polySpec))
    val zres   = Input(UInt(fracW.W)) // polynomial
    val z      = Output(UInt(spec.W.W))
  })

  val x0_5to1_0  = io.zother.x0_5to1_0
  val x1_0to2_0  = io.zother.x1_0to2_0
  val xOtherwise = io.zother.xOtherwise

  val zsgn = io.zother.zsgn

  //
  // case 1: x < 0.5 or 2 <= x. calc (ex + log2(1.man)) * ln2.
  //
  val zfrac0 = Mux(zsgn === 0.U,    io.zres, // means 0 <= xexNobias
               Mux(io.zres === 0.U, Fill(fracW, 1.U(1.W)),
                                    ~io.zres + 1.U))
  val zfull0 = Cat(io.zother.zint, zfrac0)
  val zfullShifted = (zfull0 << io.zother.zintShift)(exW+fracW-1, 0)

//   printf("cir: zres         = %b\n", io.zres     )
//   printf("cir: zfrac0       = %b\n", zfrac0      )
//   printf("cir: zfull0       = %b\n", zfull0      )
//   printf("cir: zfullShifted = %b\n", zfullShifted)

  assert(zfullShifted(exW + fracW-1) === 1.U || !xOtherwise)

  // fracW+1 width
  val ymanW1 = Mux(x0_5to1_0, Cat(1.U(1.W), io.zres), // case 2
               Mux(x1_0to2_0, Cat(io.zres, 0.U(1.W)), // case 3
                              zfullShifted(exW+fracW-1, exW-1))) // case 1
  // manW+1 width
  val cmanW1 = io.zother.constant

//   printf("cir: ymanW1 = %b\n", ymanW1)
//   printf("cir: cmanW1 = %b\n", cmanW1)
//   printf("cir: zex0   = %b\n", io.zother.zex0)

  // --------------------------------------------------------------------------
  // multiply (XXX rounding method is simple)

  assert(ymanW1.getWidth == fracW+1)
  assert(cmanW1.getWidth == manW+1)

  val zmanProd = ymanW1 * cmanW1
  val zmanProdMoreThan2 = zmanProd((fracW+1)+(manW+1)-1)
  val zmanShifted = Mux(zmanProdMoreThan2, zmanProd((fracW+1)+(manW+1)-2, fracW+1),
                                           zmanProd((fracW+1)+(manW+1)-3, fracW  ))
  val zmanRoundInc = Mux(zmanProdMoreThan2, zmanProd(fracW), zmanProd(fracW-1))
  val zmanRound = zmanShifted +& zmanRoundInc
  val zmanRoundMoreThan2 = zmanRound(manW)

  val lnman = zmanRound(manW-1, 0)
  val lnex0 = io.zother.zex0 +& (zmanProdMoreThan2 + zmanRoundMoreThan2)
  val lnex  = Mux(lnex0(exW) === 1.U, Fill(exW, 1.U(1.W)), lnex0(exW-1, 0))

  val zman0 = Wire(UInt(manW.W))
  val zex0  = Wire(UInt(exW.W))

  if(isAlwaysLn) {
    // we don't need to multiply log2(e).
    zman0 := lnman
    zex0  := lnex
  } else {
    assert(io.isln.isDefined)

    // -------------------------------------------------------------------------
    // multiply log2(e); XXX this part can be simplified because when
    // x < 0.5 || 2.0 <= x the polynomial result is already log2. Or, we can
    // just drop log2 functionality because it is not used so frequently.

    // 1 < log2e < 2. log2e.ex == 0
    val log2e = (Real.one / Real.log(Real.two))(manW).toBigInt.U((manW+1).W)

    val log2Prod = Cat(1.U(1.W), lnman) * log2e
    val log2ProdMoreThan2 = log2Prod((manW+1)+(manW+1)-1)
    val log2Shifted = Mux(log2ProdMoreThan2, log2Prod((manW+1)+(manW+1)-2, manW+1),
                                             log2Prod((manW+1)+(manW+1)-3, manW  ))
    val log2RoundInc = Mux(log2ProdMoreThan2, log2Prod(manW  ), log2Prod(manW-1))
    val log2Round = log2Shifted +& log2RoundInc
    val log2RoundMoreThan2 = log2Round(manW)

    val log2man = log2Round(manW-1, 0)
    val log2ex0 = lnex +& (log2ProdMoreThan2 + log2RoundMoreThan2)
    val log2ex  = Mux(log2ex0(exW) === 1.U, Fill(exW, 1.U(1.W)), log2ex0(exW-1, 0))

    zman0 := Mux(io.isln.get, lnman, log2man)
    zex0  := Mux(io.isln.get, lnex,  log2ex)
  }

  val znan  = io.zother.znan
  val zinf  = io.zother.zinf
  val zzero = io.zother.zzero

  val zman = Mux(znan || zinf || zzero, Cat(znan.asUInt, Fill(manW-1, 0.U(1.W))), zman0)
  val zex  = Mux(znan || zinf, Fill(exW, 1.U(1.W)),
                               Mux(zzero, 0.U(exW.W), zex0))

  val z0 = Cat(io.zother.zsgn, zex, zman)
  assert(z0.getWidth == spec.W)

  val z = enable(io.en, z0)
  io.z := ShiftRegister(z, nStage)
}

// -------------------------------------------------------------------------
//                      _     _                _
//   ___ ___  _ __ ___ | |__ (_)_ __   ___  __| |
//  / __/ _ \| '_ ` _ \| '_ \| | '_ \ / _ \/ _` |
// | (_| (_) | | | | | | |_) | | | | |  __/ (_| |
//  \___\___/|_| |_| |_|_.__/|_|_| |_|\___|\__,_|
// -------------------------------------------------------------------------

class LogGeneric(
  val isln: Boolean,
  val spec: RealSpec,
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage: MathFuncPipelineConfig,
  val dxW0 : Option[Int] = None,
  val enableRangeCheck: Boolean = true,
  val enablePolynomialRounding: Boolean = false,
) extends Module {

  val pcGap = if(stage.preCalcGap ) {1} else {0}
  val cpGap = if(stage.calcPostGap) {1} else {0}

  val nPreStage  = stage.preStage.total
  val nCalcStage = stage.calcStage.total
  val nPostStage = stage.postStage.total

  val nStage   = stage.total
  def getStage = nStage

  val polySpec = new PolynomialSpec(spec, nOrder, adrW, extraBits, dxW0,
    enableRangeCheck, enablePolynomialRounding)
  val order = polySpec.order

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val cbits = LogTableCoeff.getCBits(spec, polySpec)
  val calcW = LogTableCoeff.getCalcW(spec, polySpec)

  def getCbit  = cbits
  def getCalcW = calcW

  val io = IO(new Bundle {
    val en = Input(Bool())
    val x = Input (UInt(spec.W.W))
    val z = Output(UInt(spec.W.W))
  })

  // --------------------------------------------------------------------------

  val xdecomp = Module(new DecomposeReal(spec))
  xdecomp.io.real := io.x

  val enPCReg    = ShiftRegister(io.en,      nPreStage)
  val enPCGapReg = ShiftRegister(enPCReg,    pcGap)
  val enCPReg    = ShiftRegister(enPCGapReg, nCalcStage)
  val enCPGapReg = ShiftRegister(enCPReg,    cpGap)

  val xdecPCReg    = ShiftRegister(xdecomp.io.decomp, nPreStage)
  val xdecPCGapReg = ShiftRegister(xdecPCReg,         pcGap)
  val xdecCPReg    = ShiftRegister(xdecPCGapReg,      nCalcStage)
  val xdecCPGapReg = ShiftRegister(xdecCPReg,         cpGap)

  // --------------------------------------------------------------------------

  val logPre   = Module(new LogPreProcess (spec, polySpec, stage.preStage))
  val logTab   = Module(new LogTableCoeff (spec, polySpec, cbits))
  val logOther = Module(new LogOtherPath  (spec, polySpec, stage.calcStage))
  val logPost  = Module(new LogPostProcess(spec, polySpec, stage.postStage, false)) // can be both

  logPre.io.en  := io.en
  logPre.io.x   := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  val logPreExAdr = logPre.io.adr(logPre.io.adr.getWidth-1, logPre.io.adr.getWidth-2)
  logTab.io.en  := enPCGapReg
  logTab.io.adr := ShiftRegister(logPre.io.adr, pcGap)
  logOther.io.x := xdecPCGapReg

  // after preprocess
  assert(logPre.io.adr === 0.U               || enPCReg)
  assert(logPre.io.dx.getOrElse(0.U) === 0.U || enPCReg)
  assert(logTab.io.cs.asUInt === 0.U         || enPCGapReg)

  // --------------------------------------------------------------------------

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, cbits, stage.calcStage))

  if(order != 0) {
    polynomialEval.io.dx.get := ShiftRegister(logPre.io.dx.get, pcGap)
  }
  polynomialEval.io.coeffs.cs := logTab.io.cs.cs

  val polynomialResultCPGapReg = ShiftRegister(polynomialEval.io.result, cpGap)

  logPost.io.en     := enCPGapReg
  logPost.io.zother := ShiftRegister(logOther.io.zother, cpGap)
  logPost.io.zres   := polynomialResultCPGapReg
  logPost.io.isln.get := isln.B

  io.z := logPost.io.z
}


