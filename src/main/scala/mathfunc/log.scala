//% @file log.scala
//
// log2 and ln function
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
// It also checks if x.ex == 0 or -1 and add it to MSB side of the address.
// To distinguish 0, -1, and others, we need 3 states = 2 bits.

class Log2PreProcess(
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

  val ex = io.x.ex
  val exAdr = Mux(ex === exBias.U, 2.U, (ex === (exBias - 1).U).asUInt)

  // if xexNobias == -1, (x is in [1/2,1)), use 1<<manW - x.man.
  // otherwise, use x.man.
  val manPos = io.x.man
  val manNeg = ~(io.x.man) + 1.U // === (1<<manW).U - x.man

  val adr0 = Cat(exAdr,
    Mux(ex === (exBias-1).U, manNeg(manW-1, dxW), manPos(manW-1, dxW)))
  val adr  = enable(io.en, adr0)
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx0 = Mux(ex === (exBias-1).U,
      Cat(~manNeg(dxW-1), manNeg(dxW-2, 0)),
      Cat(~manPos(dxW-1), manPos(dxW-2, 0)))
    val dx = enable(io.en, dx0)
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

class Log2TableCoeff(
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

  val log2 = (a:Double) => {log(a) / log(2.0)}

  if(order == 0) {
    val exadr = io.adr(adrW+2-1, adrW)
    val adr   = io.adr(adrW-1, 0)

    val tableNormalI  = MathFuncLogSim.logNormalTableGeneration(spec, order, adrW, extraBits)
    val tableSmallNeg = MathFuncLogSim.logSmallNegativeTableGeneration(spec, order, adrW, extraBits)
    val tableSmallPos = MathFuncLogSim.logSmallPositiveTableGeneration(spec, order, adrW, extraBits)

    val cbitNormal   = tableNormalI .getCBitWidth(/*sign mode = */1) // 2's complement and no sign bit
    val cbitSmallNeg = tableSmallNeg.getCBitWidth(/*sign mode = */1)
    val cbitSmallPos = tableSmallPos.getCBitWidth(/*sign mode = */1)

    val cbit = Seq(cbitNormal, cbitSmallNeg, cbitSmallPos).reduce(
      (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )
    val tableIs = VecInit(Seq(
        tableNormalI .getVectorWithWidth(cbit, /*sign mode = */1),
        tableSmallNeg.getVectorWithWidth(cbit, /*sign mode = */1),
        tableSmallPos.getVectorWithWidth(cbit, /*sign mode = */1)
      ))
    val tableI = tableIs(exadr)
    val coeff = getSlices(tableI(adr), cbit)

    assert(maxCbit(0) == fracW)
    assert(coeff(0).getWidth == fracW, f"coeff = ${coeff(0).getWidth}, fracW = ${fracW}")

    io.cs.cs(0) := enable(io.en, coeff(0))

  } else {

    // split address
    val exadr = io.adr(adrW+2-1, adrW)
    val adr   = io.adr(adrW-1, 0)

    // -----------------------------------------------------------------------
    // default table
    val tableNormalI  = MathFuncLogSim.logNormalTableGeneration(spec, order, adrW, extraBits)
    val tableSmallNeg = MathFuncLogSim.logSmallNegativeTableGeneration(spec, order, adrW, extraBits)
    val tableSmallPos = MathFuncLogSim.logSmallPositiveTableGeneration(spec, order, adrW, extraBits)

    val cbitNormal   = tableNormalI .getCBitWidth(/*sign mode = */0)
    val cbitSmallNeg = tableSmallNeg.getCBitWidth(/*sign mode = */0)
    val cbitSmallPos = tableSmallPos.getCBitWidth(/*sign mode = */0)

    val cbit = Seq(cbitNormal, cbitSmallNeg, cbitSmallPos).reduce(
      (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )
    val tableIs = VecInit(Seq(
        tableNormalI .getVectorWithWidth(cbit, /*sign mode = */0),
        tableSmallNeg.getVectorWithWidth(cbit, /*sign mode = */0),
        tableSmallPos.getVectorWithWidth(cbit, /*sign mode = */0)
      ))
    val tableI = tableIs(exadr)
    val coeff = getSlices(tableI(adr), cbit)

    val coeffs = Wire(new TableCoeffInput(maxCbit))

    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(cbit(i) <= maxCbit(i))

      if(0 < diffWidth) {
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

object Log2TableCoeff {
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
      val tableNormalI = MathFuncLogSim.logNormalTableGeneration(
        spec, order, adrW, extraBits)
      val cbitNormal   = tableNormalI.cbit

      val tableSmallPos = MathFuncLogSim.logSmallPositiveTableGeneration(
        spec, order, adrW, extraBits)
      val cbitSmallPos = tableSmallPos.cbit

      val tableSmallNeg = MathFuncLogSim.logSmallNegativeTableGeneration(
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
      val tableNormalI = MathFuncLogSim.logNormalTableGeneration(
        spec, order, adrW, extraBits)
      val calcWidthNormal   = tableNormalI.calcWidth

      val tableSmallPos = MathFuncLogSim.logSmallPositiveTableGeneration(
        spec, order, adrW, extraBits)
      val calcWidthSmallPos = tableSmallPos.calcWidth

      val tableSmallNeg = MathFuncLogSim.logSmallNegativeTableGeneration(
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
// This has two taylor pathes.
//

class Log2NonTableOutput(val spec: RealSpec, val polySpec: PolynomialSpec) extends Bundle {
  // always required
  val zsgn  = Output(UInt(1.W))
  val zIsNonTable = Output(Bool())
  val znan  = Output(Bool())
  val zinf  = Output(Bool())
  val zzero = Output(Bool())
  // taylor result & special value result
  val zman  = Output(UInt(polySpec.fracW.W))
  val zex   = Output(UInt(spec.exW.W))
  // default table
  val zint  = Output(UInt(spec.exW.W))
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

  val fracW     = polySpec.fracW
  val extraBits = polySpec.extraBits

  val io = IO(new Bundle {
    val x      = Flipped(new DecomposedRealOutput(spec))
    // exadr = 1: xexNobias == -1, x in [1/2, 1)
    // exadr = 2: xexNobias ==  0, x in [1, 2)
    // exadr = 0: others
    val exadr  = Input(UInt(2.W))
    val zother = new Log2NonTableOutput(spec, polySpec)
    val xmanbp = Output(UInt(log2Up(manW).W))
  })

  val xmanNeg = (1L<<manW).U - io.x.man

  val xmanbpPos =  manW.U - PriorityEncoder(Reverse(io.x.man))
  val xmanbpNeg = (1+manW).U - PriorityEncoder(Reverse(xmanNeg)) // the width changed

  val taylorThreshold = MathFuncLogSim.calcTaylorThreshold(spec)
  val invln2   = math.round((1.0 / log(2.0)) * (1L << fracW)).toLong.U((fracW+1).W)
  val oneThird = math.round((1.0 / 3.0)      * (1L << fracW)).toLong.U( fracW   .W)

  val isTaylorSmallPos = (io.exadr === 2.U) && (io.x.man < (1L << (manW-taylorThreshold)).U)
  val isTaylorSmallNeg = (io.exadr === 1.U) && (xmanNeg  < (1L << (manW-taylorThreshold+1)).U)

//   printf("cir: xmanNeg           = %b\n", xmanNeg)
//   printf("cir: xmanbpPos         = %d\n", xmanbpPos)
//   printf("cir: xmanbpNeg         = %d\n", xmanbpNeg)
//   printf("cir: isTaylorSmallPos  = %d\n", isTaylorSmallPos)
//   printf("cir: isTaylorSmallNeg  = %d\n", isTaylorSmallNeg)

  // here you cannot use isTaylor because this value is used even in the non-taylor path
  val xmanbp0 = Mux(io.exadr === 2.U, xmanbpPos, xmanbpNeg)
  io.xmanbp := ShiftRegister(xmanbp0, nStage)

  // --------------------------------------------------------------------------
  // taylor x in [1, 2)
  //
  // log(1+x) = x/ln(2) - x^2/2ln(2) + x^3/3ln(2) + O(x^4)
  //          = x(1 - x/2 + x^2/3) / ln(2)

  // 1 - x/2 < 1
  val oneMinusHalfx = if(extraBits == 0) {
    (1L << fracW).U - io.x.man(manW-1, 1)
  } else if(extraBits == 1) {
    (1L << fracW).U - io.x.man
  } else {
    (1L << fracW).U - Cat(io.x.man, 0.U((extraBits-1).W))
  }

  // x^2/3
  val xsqPos       = io.x.man * io.x.man
  val xsqThirdPos0 = xsqPos * oneThird
  val xsqThirdPos  = xsqThirdPos0 >> (manW * 2).U

  // 1 - x/2 + x^2/3 < 1, x < 2^-8
  val taylorTermPos = oneMinusHalfx +& xsqThirdPos
  // x < x/ln2 ~ x * 1.44 < 2x
  val convTermPos   = invln2 * io.x.man

  // x < x/ln2 * (1 - x/2 + x^2/3) < 2x
  val resPosProd      = (convTermPos * taylorTermPos) >> xmanbpPos
  val resPosMoreThan2 = resPosProd(fracW*2)
  val resPosShifted   = Mux(resPosMoreThan2,
    resPosProd(resPosProd.getWidth-1, fracW    ),
    resPosProd(resPosProd.getWidth-1, fracW - 1))
  val resPosProdInc   = Mux(resPosMoreThan2,
    resPosProd(fracW - 1),
    resPosProd(fracW - 2))
  val resPosProdRounded = resPosShifted +& resPosProdInc
  val resPosShiftedMoreThan2 = resPosProdRounded(fracW+1)

  val zmanTaylorPos = resPosProdRounded(fracW-1, 0)
  val zexTaylorPos  = (exBias - manW - 1).U + xmanbpPos + resPosMoreThan2 + resPosShiftedMoreThan2

  // --------------------------------------------------------------------------
  // taylor x in [1/2, 1)
  //
  // log(1-x) = -x/ln(2) - x^2/2ln(2) - x^3/3ln(2) - O(x^4)
  //          = -x(1 + x/2 + x^2/3) * (1 / ln(2))

  // 1 + x/2 > 1
  val onePlusHalfx = Wire(UInt((fracW+1).W))
  if(extraBits > 2) {
    onePlusHalfx := (1L << fracW).U + Cat(xmanNeg, 0.U((extraBits-2).W))
  } else if (extraBits == 2) {
    onePlusHalfx := (1L << fracW).U + xmanNeg
  } else {
    onePlusHalfx := (1L << fracW).U + (xmanNeg >> (2-extraBits))
  }
//   printf("cir: onePlusHalfx  = %b\n", onePlusHalfx     )

  // x^2/3
  val xsqNeg      = xmanNeg * xmanNeg
  val xsqThirdNeg = (xsqNeg * oneThird) >> ((manW+1) * 2).U
//   printf("cir: xsqNeg        = %b\n", xsqNeg     )
//   printf("cir: xsqThirdNeg   = %b\n", xsqThirdNeg)

  // 1 + x/2 + x^2/3 > 1
  val taylorTermNeg = onePlusHalfx + xsqThirdNeg
  val convTermNeg   = xmanNeg * invln2

//   printf("cir: taylorTermNeg = %b\n", taylorTermNeg)
//   printf("cir: convTermNeg   = %b\n", convTermNeg  )
//
//   printf("cir: xmanbpNeg   = %d\n", xmanbpNeg  )
  val resNegProd      = (convTermNeg * taylorTermNeg) >> xmanbpNeg
  val resNegMoreThan2 = resNegProd(fracW + fracW)
  val resNegShifted   = Mux(resNegMoreThan2,
    resNegProd(resNegProd.getWidth-1, fracW    ),
    resNegProd(resNegProd.getWidth-1, fracW - 1))
  val resNegProdInc   = Mux(resNegMoreThan2,
    resNegProd(fracW - 1),
    resNegProd(fracW - 2))
  val resNegProdRounded = resNegShifted + resNegProdInc
  val resNegShiftedMoreThan2 = resNegProdRounded(fracW+1)

//   printf("cir: resNegProdRounded  = %b\n", resNegProdRounded     )

  val zmanTaylorNeg = resNegProdRounded(fracW-1, 0)
  val zexTaylorNeg  = (exBias - manW - 2).U + xmanbpNeg + resNegMoreThan2 + resNegShiftedMoreThan2
//   printf("cir: zmanTaylor = %b\n", zmanTaylorNeg)
//   printf("cir: zexTaylor  = %d\n", zexTaylorNeg)

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
  val xexPos = io.x.ex - exBias.U     // ==  exNobias
  val xexNeg = (exBias-1).U - io.x.ex // == -exNobias - 1 == -(ex - exBias) - 1 = -ex + exBias - 1

  // integer part of z (fractional part is calculated in polynomial module)
  val zint = Mux(io.x.ex >= exBias.U, xexPos, xexNeg)
  io.zother.zint := ShiftRegister(zint, nStage)

  // --------------------------------------------------------------------------
  // check special value
  //
  // log2(nan) ->  nan
  // log2(inf) ->  inf
  // log2(0)   -> -inf
  // log2(2)   ->  1
  // log2(1)   ->  0
  // log2(1/2) -> -1
  // log2(-|x|) -> nan

  val xmanAllZero = !io.x.man.orR
  val znan  = io.x.nan || io.x.sgn === 1.U
  val zinf  = io.x.inf || io.x.zero
  val zzero = xmanAllZero && io.x.ex === exBias.U
  val xtwo  = xmanAllZero && io.x.ex === (exBias+1).U
  val xhalf = xmanAllZero && io.x.ex === (exBias-1).U

  io.zother.znan  := ShiftRegister(znan,  nStage)
  io.zother.zinf  := ShiftRegister(zinf,  nStage)
  io.zother.zzero := ShiftRegister(zzero, nStage)

  val zIsNonTable = xtwo || xhalf || isTaylorSmallPos || isTaylorSmallNeg
  io.zother.zIsNonTable := ShiftRegister(zIsNonTable, nStage)

  val zsgn = io.x.ex < exBias.U
  io.zother.zsgn := ShiftRegister(zsgn,  nStage)

  // --------------------------------------------------------------------------
  // merge Taylor results and special values

  val zex0  = Mux(isTaylorSmallPos, zexTaylorPos,
              Mux(isTaylorSmallNeg, zexTaylorNeg,
              Mux(xtwo,  exBias.U(exW.W), /*xhalf = */ (exBias-1).U(exW.W))))
  val zman0 = Mux(isTaylorSmallPos, zmanTaylorPos, zmanTaylorNeg)
  // here, adding 0es at the LSB of nan/inf/zero does not affect to the result
  // because postprocess only does rounding. zero bits does not change rounding.

//   printf("cir: zman0 = %b\n", zman0)
//   printf("cir: zex0  = %d\n", zex0)

  io.zother.zman := ShiftRegister(zman0, nStage)
  io.zother.zex  := ShiftRegister(zex0,  nStage)
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
  def getStage = nStage

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order
  val extraBits = polySpec.extraBits

  val log2 = (a:Double) => {log(a) / log(2.0)}

  val io = IO(new Bundle {
    val en = Input(UInt(1.W))
    val isln   = Input(Bool()) // If log_e, true. else if log2, false.
    val x      = Flipped(new DecomposedRealOutput(spec))
    val zother = Flipped(new Log2NonTableOutput(spec, polySpec))
    val exadr  = Input(UInt(2.W))
    val xmanbp = Input(UInt(log2Up(manW).W))
    val zres   = Input(UInt(fracW.W)) // polynomial
    val z      = Output(UInt(spec.W.W))
  })

  // --------------------------------------------------------------------------
  // postprocess polynomial result; x is in [1, 2) and [1/2, 1)

  val zSmallExBase = (exBias-manW).U + io.xmanbp
  val zSmallPosEx = Mux(io.zres(fracW-1) === 1.U, zSmallExBase      ,
                    Mux(io.zres(fracW-2) === 1.U, zSmallExBase - 1.U,
                                                  zSmallExBase - 2.U))
  val zSmallNegEx = Mux(io.zres(fracW-1) === 1.U, zSmallExBase - 1.U,
                    Mux(io.zres(fracW-2) === 1.U, zSmallExBase - 2.U,
                                                  zSmallExBase - 3.U))

  val zSmallManBase = Cat(io.zres, 0.U(3.W))
  val zSmallMan0 = Mux(io.zres(fracW-1) === 1.U, zSmallManBase(fracW+2-1, 2),
                   Mux(io.zres(fracW-2) === 1.U, zSmallManBase(fracW+1-1, 1),
                                                 zSmallManBase(fracW+0-1, 0)))
//   printf("cir: zSmallManBase = %b\n", zSmallManBase)
  assert(zSmallMan0.getWidth == fracW)
  // exadr = 1: xexNobias == -1, x in [1/2, 1)
  // exadr = 2: xexNobias ==  0, x in [1, 2)
  // exadr = 0: others
  val zSmallEx0 = Mux(io.exadr === 2.U, zSmallPosEx, zSmallNegEx)

  // --------------------------------------------------------------------------
  // postprocess polynomial result; x is in [2, inf) or (0, 1/2]

  val zLargeInt     = io.zother.zint
  val zLargeFracPos = io.zres
  val zLargeFracNeg = Mux(io.zres === 0.U, maskL(fracW).U(fracW.W), ~io.zres + 1.U)
  val zLargeFrac    = Mux(io.x.ex >= exBias.U, zLargeFracPos, zLargeFracNeg)

//   printf("cir: zLargeInt  = %b\n", zLargeInt)
//   printf("cir: zLargeFrac = %b\n", zLargeFrac)

  // z = log2(2^xex * 1.xman)
  //   = xex + log2(1.xman)
  // if xex < 0, we need to subtract log2(1.xman) from |xex| and set sgn to 1

  val zLargeFull     = Cat(zLargeInt, zLargeFrac)
  val zLargeFullPrec = PriorityEncoder(Reverse(zLargeFull)) // 0es at the MSBs
  val zLargeShiftW   = Mux(zLargeFull(exW+fracW-1) === 1.U, 0.U, zLargeFullPrec)
  val zLargeShifted  = zLargeFull << zLargeShiftW // W = fracW+exW
  // here we multiplied zInt + zFrac by 8 - zLargeShiftW.

//   printf("cir: zLargeFull    = %b\n", zLargeFull   )
//   printf("cir: zLargeShiftW  = %b\n", zLargeShiftW )
//   printf("cir: zLargeShifted = %b\n", zLargeShifted)

  // if log is active && exadr == 0 && it uses table, zLargeShifted should have
  // a value.
  assert(zLargeShifted(fracW+exW-1) === 1.U || io.exadr =/= 0.U ||
         io.zother.zIsNonTable || !io.en)

  assert(extraBits < exW) // normally this holds. normally.
  val zLargeMan0 = zLargeShifted(fracW+exW-2, exW-1)
  val zLargeEx0  = (exBias + exW - 1).U - zLargeShiftW

  // --------------------------------------------------------------------------
  // rounding

  val log2xEx0  = Mux(io.zother.zIsNonTable, io.zother.zex ,
                  Mux(io.exadr === 0.U, zLargeEx0,  zSmallEx0))
  val log2xMan0 = Mux(io.zother.zIsNonTable, io.zother.zman,
                  Mux(io.exadr === 0.U, zLargeMan0, zSmallMan0))
//   printf("log2xEx0  = %b\n", log2xEx0 )
//   printf("log2xMan0 = %b\n", log2xMan0)

  // -------------------------------------------------------------------------
  // ln(x)

  // ln2 = 0.6931... < 1
  val ln2 = math.round(log(2.0) * (1<<(fracW+1))).toLong.U((fracW+1).W)
  assert(ln2(fracW) === 1.U)

  val lnxProd      = Cat(1.U(1.W), log2xMan0) * ln2
  val lnxMoreThan2 = lnxProd((1+fracW)*2-1).asBool
  val lnxShifted   = Mux(lnxMoreThan2, lnxProd(fracW*2,   fracW+extraBits+1),
                                       lnxProd(fracW*2-1, fracW+extraBits  ))
  val lnxInc       = Mux(lnxMoreThan2, lnxProd(fracW+extraBits),
                                       lnxProd(fracW+extraBits-1))
  val lnxRounded   = lnxShifted +& lnxInc
  val lnxMoreThan2AfterRound = lnxRounded(manW)
  val lnxEx0       = log2xEx0 - 1.U + lnxMoreThan2 + lnxMoreThan2AfterRound

  // -------------------------------------------------------------------------
  // log2(x)

  val log2xManRound = Wire(UInt(manW.W))
  val log2xEx       = Wire(UInt(exW.W))

  if(extraBits == 0) {
    log2xManRound := log2xMan0
    log2xEx := log2xEx0
  } else {
    val log2xManRound0 = log2xMan0(fracW-1, extraBits) +& log2xMan0(extraBits-1)
    log2xManRound := log2xManRound0(manW-1, 0)
    log2xEx := log2xEx0 + log2xManRound0(manW)
  }

  val zMan = Mux(io.zother.znan || io.zother.zinf || io.zother.zzero,
                 Cat(io.zother.znan.asUInt, Fill(manW-1, 0.U)),
                 Mux(io.isln, lnxRounded(manW-1, 0), log2xManRound))
  val zEx  = Mux(io.zother.znan || io.zother.zinf, Fill(exW, 1.U(1.W)),
             Mux(io.zother.zzero, 0.U(exW.W),
             Mux(io.isln, lnxEx0, log2xEx)))

  // --------------------------------------------------------------------------
  // select the correct result

  val zSgn = io.zother.zsgn
  val z0   = Cat(zSgn, zEx, zMan)

  assert(zEx .getWidth == exW)
  assert(zMan.getWidth == manW)
  assert(z0  .getWidth == spec.W)

  val z = enable(io.en, z0)

  io.z   := ShiftRegister(z, nStage)
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

  val polySpec = new PolynomialSpec(spec, nOrder, adrW, extraBits,
    enableRangeCheck, enablePolynomialRounding)
  val order = polySpec.order

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val cbits = Log2TableCoeff.getCBits(spec, polySpec)
  val calcW = Log2TableCoeff.getCalcW(spec, polySpec)

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

  val logPre   = Module(new Log2PreProcess (spec, polySpec, stage.preStage))
  val logTab   = Module(new Log2TableCoeff (spec, polySpec, cbits))
  val logOther = Module(new Log2OtherPath  (spec, polySpec, stage.calcStage))
  val logPost  = Module(new Log2PostProcess(spec, polySpec, stage.postStage))

  logPre.io.en  := io.en
  logPre.io.x   := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  val logPreExAdr = logPre.io.adr(logPre.io.adr.getWidth-1, logPre.io.adr.getWidth-2)
  logTab.io.en  := enPCGapReg
  logTab.io.adr := ShiftRegister(logPre.io.adr, pcGap)
  logOther.io.x := xdecPCGapReg
  logOther.io.exadr := ShiftRegister(logPreExAdr, pcGap)

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
  logPost.io.isln   := isln.B
  logPost.io.x      := xdecCPGapReg
  logPost.io.exadr  := ShiftRegister(logPreExAdr, pcGap + nCalcStage + cpGap)
  logPost.io.xmanbp := ShiftRegister(logOther.io.xmanbp, cpGap)

  io.z := logPost.io.z
}


