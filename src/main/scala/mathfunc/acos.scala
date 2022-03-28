//% @file acos.scala
//
// ACos function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

import java.lang.Math.scalb
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
import rial.arith._

import rial.mathfunc._

// Calculate acos(x) if x is in range [-1, 1].
// otherwise, 0 or pi for positive and negative x, respectively.
//
// for x < 0.5, pi/2 - acos(x) is approximated by polynomial.
// for x > 0.5, (pi/2 - acos(x)) / 2 is approximated by polynomial.
// for x close to 1, use puiseux series:
//   let y = 1-x
//   acos(y) = sqrt(2y) * (1 + y/12 + 3y^2/160 + 5y^3/896 + 35y^4/18432 + O(y^5))
//
// In case of Puiseux series, the condition is:
//   35y^4/18432 < 2^-23
//     y^4/526.628... < 2^-23
//     y^4 < 2^-14
//     y   < 2^-4
// So if 1-2^-4 < x, that means that x.man(22,20).andR === 1, use Puiseux series.
//
// Normal table calculates f(x) = pi/2 - acos(x). So
//   acos( x) = pi/2 - pi/2 + acos(x) = pi/2 - [pi/2 - acos(x)]
//   acos(-x) = pi - acos(x)          = pi/2 + [pi/2 - acos(x)]

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------
//
// To use puiseux series, we need to calculate sqrt(2y).
// Also, for the cases where xex < 2, it uses a table pi/2 - acos([0, 0.5]).
// So we need to align xman before calculating acos dx and adr.
//
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
  val exAdrW    = 1

  val io = IO(new Bundle {
    val en  = Input (UInt(1.W))
    val x   = Flipped(new DecomposedRealOutput(spec))
    val adr = Output(UInt((exAdrW+adrW).W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }

    // if x is close to 1, we switch to puiseux series expansion. To do that,
    // we need to calculate sqrt(2*(1-x)). This output indicates we need to
    // calc sqrt.
    val useSqrt = if(order != 0) {Some(Output(Bool()))           } else {None}
    val yex     = if(order != 0) {Some(Output(UInt(spec.exW.W))) } else {None}
    val yman    = if(order != 0) {Some(Output(UInt(spec.manW.W)))} else {None}
  })

  val xsgn = enable(io.en, io.x.sgn)
  val xex  = enable(io.en, io.x.ex )
  val xman = enable(io.en, io.x.man)

  if(order == 0) {

    val adrACosShift0   = (exBias - 2 + 1).U - xex
    val adrACosShift    = adrACosShift0(log2Up(manW+1), 0)
    val adrACosShiftOut = adrACosShift0 > (manW+1).U
    val adrACosAligned0 = Mux(adrACosShiftOut, 0.U, Cat(1.U(1.W), xman) >> adrACosShift)
    val adrACosAligned  = adrACosAligned0(manW-1, 0)

    val exAdrACos = ((exBias - 1).U(exW.W) === xex)
    val adr0ACos  = Mux(exAdrACos, xman(manW-1, dxW), adrACosAligned(manW-1, dxW))
    val adrACos   = Cat(exAdrACos.asUInt, adr0ACos)

    val adr = enable(io.en, adrACos)
    io.adr := ShiftRegister(adr, nStage)

  } else { // order != 0
    assert(io.dx     .isDefined)
    assert(io.useSqrt.isDefined)
    assert(io.yex    .isDefined)
    assert(io.yman   .isDefined)

    // ------------------------------------------------------------------------
    // useSqrt condition is:
    //   0          < 1-x < 2^-4
    //   1.0 - 2^-4 < x   < 1
    // in this range, x.ex == -1.
    //
    // |x|11111110|111xxx...x|
    //            ^  ^
    //            |  +- 2^-4 bit
    //            +- hidden bit = 2^-1

    val puiseuxThreshold = MathFuncACosSim.calcPuiseuxThreshold(manW, fracW-manW, 0)
    // 1 for hidden bit, 1 for table exponent
    val puiseuxMSBs = puiseuxThreshold.abs - 1 - 1
    val useSqrt = enable(io.en, (xex === (exBias-1).U) && xman(manW-1, manW-puiseuxMSBs).andR)
    io.useSqrt.get := ShiftRegister(useSqrt, nStage)

    // ---------------------------------------------------------------------------
    // To use sqrt table, we need to determine the mantissa and exponent of
    // y = 1 - x.
    // Since the value of y will be used in the non-table path later, it should
    // also be passed to ACosOtherPath.

    val oneMinusX     = ~xman + 1.U // 1<<manW - xman
    val oneMinusXLen  = manW.U - PriorityEncoder(Reverse(oneMinusX))
    val oneMinusXSftW = log2UpL(manW)
    val oneMinusXSft  = ((manW+1).U - oneMinusXLen)(oneMinusXSftW-1, 0)
    val ymanW1        = (oneMinusX << oneMinusXSft)(manW, 0)
    val yex           = oneMinusXLen + (exBias - manW - 2).U(exW.W)

    io.yex.get  := ShiftRegister(yex(exW-1, 0),     nStage)
    io.yman.get := ShiftRegister(ymanW1(manW-1, 0), nStage)

//     printf("cir: yex    = %d\n", yex.zext - exBias.S)
//     printf("cir: ymanW1 = %b\n", ymanW1)

    assert(!io.en || (ymanW1(manW) === 1.U))

    // ---------------------------------------------------------------------------
    // calc adr and dx

    // exAdr shows:
    // - 0: 0.0 < x < 0.5 (xex <  exBias-1)
    // - 1: 0.5 < x < 1.0 (xex == exBias-1)

    val adrACosShift0   = (exBias - 2 + 1).U - xex
    val adrACosShift    = adrACosShift0(log2Up(manW+1), 0)
    val adrACosShiftOut = adrACosShift0 > (manW+1).U
    val adrACosAligned0 = Mux(adrACosShiftOut, 0.U, Cat(1.U(1.W), xman) >> adrACosShift)
    val adrACosAligned  = adrACosAligned0(manW-1, 0)

    val exAdrACos = ((exBias - 1).U(exW.W) === xex)
    val adr0ACos  = Mux(exAdrACos, xman(manW-1, dxW), adrACosAligned(manW-1, dxW))
    val adrACos   = Cat(exAdrACos.asUInt, adr0ACos)

    // we need sqrt(2y). Sqrt requires only 1 LSB to distinguish x in 1~2 and 2~4.
    // To pass yex+1, the last bit is inverted.
    val exAdrSqrt = ~yex(0)
    val adrSqrt   = Cat(exAdrSqrt, ymanW1(manW-1, dxW))

    val adr = enable(io.en, Mux(useSqrt, adrSqrt, adrACos))
    io.adr := ShiftRegister(adr, nStage)

    val dxACos = Mux(exAdrACos, Cat(          ~xman(dxW-1),           xman(dxW-2, 0)),
                                Cat(~adrACosAligned(dxW-1), adrACosAligned(dxW-2, 0)))
    val dxSqrt = Cat(~ymanW1(dxW-1), ymanW1(dxW-2, 0))
    val dx   = enable(io.en, Mux(useSqrt, dxSqrt, dxACos))
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
  val maxCbit  : Seq[Int], // max coeff width among all math funcs
) extends Module {

  val manW   = spec.manW
  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order
  val exAdrW = 1

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input  (UInt((exAdrW+adrW).W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  val exAdr = io.adr(adrW)
  val adr   = io.adr(adrW - 1, 0)

  if(order == 0) {

    val cbit = MathFuncACosSim.acosTableGeneration( spec, order, adrW, manW, fracW )
      .map( t => {t.getCBitWidth(/*sign mode = */1)} )
      .reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

    val tableIs = VecInit(
      MathFuncACosSim.acosTableGeneration( spec, order, adrW, manW, fracW ).map(t => {
        t.getVectorWithWidth(cbit, /*sign mode = */ 1)
      })
    )
    val tableI = tableIs(exAdr)
    val coeff = getSlices(tableI(adr), cbit)

    assert(maxCbit(0) == fracW)
    assert(coeff(0).getWidth == fracW, f"coeff = ${coeff(0).getWidth}, fracW = ${fracW}")

    io.cs.cs(0) := enable(io.en, coeff(0))

  } else {

    val cbit = MathFuncACosSim.acosTableGeneration( spec, order, adrW, manW, fracW )
      .map( t => {t.getCBitWidth(/*sign mode = */0)} )
      .reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

    val tableIs = VecInit(
      MathFuncACosSim.acosTableGeneration( spec, order, adrW, manW, fracW ).map(t => {
        t.getVectorWithWidth(cbit, /*sign mode = */ 0)
      })
    )
    val tableI = tableIs(exAdr)
    val coeff = getSlices(tableI(adr), cbit)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(0 <= diffWidth)

      val ci  = if(diffWidth != 0) {
        Cat(0.U(diffWidth.W), coeff(i))
      } else {
        coeff(i) // no need to extend; this is the largest value in all the tables
      }
      coeffs.cs(i) := ci
    }
    io.cs := enable(io.en, coeffs)
  }
}

object ACosTableCoeff {
  def getCBits(
    spec:     RealSpec,
    polySpec: PolynomialSpec,
  ): Seq[Int] = {

    val order     = polySpec.order
    val adrW      = polySpec.adrW
    val extraBits = polySpec.extraBits
    val fracW     = polySpec.fracW

    if(order == 0) {
      return Seq(fracW)
    } else {
      return MathFuncACosSim.acosTableGeneration( spec, order, adrW, spec.manW, fracW ).
        map( t => {t.getCBitWidth(/*sign mode = */0)} ).
        reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )
    }
  }
  def getCalcW(
    spec:     RealSpec,
    polySpec: PolynomialSpec,
  ): Seq[Int] = {

    val order     = polySpec.order
    val adrW      = polySpec.adrW
    val extraBits = polySpec.extraBits
    val fracW     = polySpec.fracW

    if(order == 0) {
      return Seq(fracW)
    } else {
      return MathFuncACosSim.acosTableGeneration( spec, order, adrW, spec.manW, fracW ).
        map( t => {t.calcWidth} ).
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

// for x close to 1, use Puiseux series:
//   acos(1-x) = sqrt(2x) * (1 + x/12 + 3x^2/160 + 5x^3/896 + O(x^4))
//                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//                           this part is used
//
// for |x| > 1, set IsPuiseux = true and returns zero. In the postprocess,
// the result of Puiseux will be subtracted from Pi if x.sgn == 1.
//
// If the float spec has narrow manW and polynomial table only has 0th order
// polynomial, then it calculates exponent using a table.
//

class ACosNonTableOutput(val spec: RealSpec, val polySpec: PolynomialSpec) extends Bundle {
  val xsgn          = Output(UInt(1.W))
  val xLessThanHalf = Output(Bool())
  val znan          = Output(Bool())
  val zzero         = Output(Bool())
  val zex           = Output(UInt(spec.exW.W))
  val zman          = if(polySpec.order != 0) {Some(Output(UInt(spec.manW.W)))} else {None}
  // needs pi - z if xsgn
  val zIsPuiseux    = if(polySpec.order != 0) {Some(Output(Bool()))} else {None}
}

class ACosOtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order

  val io = IO(new Bundle {
    val x       = Flipped(new DecomposedRealOutput(spec))
    val yex     = if(order != 0) {Some(Input(UInt(exW.W))) } else {None}
    val yman    = if(order != 0) {Some(Input(UInt(manW.W)))} else {None}
    val useSqrt = if(order != 0) {Some(Input(Bool()))      } else {None}
    val zother  = new ACosNonTableOutput(spec, polySpec)
  })

  // --------------------------------------------------------------------------
  // Special value (x is nan or |x| is larger than 1)

  val znan  = io.x.nan
  val xMoreThan1 = exBias.U <= io.x.ex

  io.zother.znan          := ShiftRegister(znan,       nStage)
  io.zother.zzero         := ShiftRegister(xMoreThan1, nStage)
  io.zother.xsgn          := ShiftRegister(io.x.sgn,   nStage)
  io.zother.xLessThanHalf := ShiftRegister(io.x.ex < (exBias-1).U, nStage)

  if(polySpec.order == 0) { // for BF16

    val cbit = MathFuncACosSim.acosExTableGeneration( spec, order, adrW, manW, fracW )
      .map( t => {t.getCBitWidth(/*sign mode = */1)} )
      .reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

    val tableIs = VecInit(
      MathFuncACosSim.acosExTableGeneration( spec, order, adrW, manW, fracW ).map(t => {
        t.getVectorWithWidth(cbit, /*sign mode = */ 1)
      })
    )

    val exAdr = ((exBias - 1).U(exW.W) === io.x.ex).asUInt
    val tableI = tableIs(exAdr)
    val coeff = getSlices(tableI(io.x.man), cbit)

    val diffWidth = exW - cbit(0)
    assert(0 <= diffWidth)

    val ci = if(diffWidth == 0) { coeff(0) } else {
      Cat(0.U(diffWidth.W), coeff(0))
    }

    io.zother.zex  := ShiftRegister(ci, nStage)

  } else { // for FP32

    // if |x| > 1, return 0 as puiseux result.

    io.zother.zIsPuiseux.get := ShiftRegister(xMoreThan1 || io.useSqrt.get, nStage)

    // --------------------------------------------------------------------------
    // helper function to multiply values

    val multiply = (xmanW1: UInt, ymanW1: UInt) => {
      assert(xmanW1.getWidth == manW + 1)
      assert(ymanW1.getWidth == manW + 1)

      val zProd      = xmanW1 * ymanW1
      val zMoreThan2 = zProd((manW+1)*2 - 1)
      val zShifted   = Mux(zMoreThan2, zProd((manW+1)*2-2, manW+1), zProd((manW+1)*2-3, manW))
      val zRounded   = zShifted +& Mux(zMoreThan2, zProd(manW), zProd(manW-1))
      val zMoreThan2AfterRound = zRounded(manW)
      val zmanW1     = Mux(zMoreThan2AfterRound, Cat(1.U(1.W), 0.U(manW.W)),
                                                 Cat(1.U(1.W), zRounded(manW-1, 0)))
      val zexInc     = zMoreThan2 + zMoreThan2AfterRound
      (zexInc, zmanW1)
    }

    val shiftOut = log2UpL(manW) // for FP32, log2Up(23) = 5

    // --------------------------------------------------------------------------
    // Puiseux
    //
    // let y = 1 - x
    // acos(1-y) = sqrt(2y) * (1 + y/12 + 3y^2/160 + 5y^3/896 + O(y^4))
    //           = sqrt(2y) * ((1 + 1/3 * 2^-2 * y) + 3y^2/160 * (1 + 25/21 * 2^-2 * y))
    //                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    //                        calc this part

    val yex    = io.yex.get
    val ymanW1 = Cat(1.U(1.W), io.yman.get)

    // ----------------------------------------------------------------------
    // 2^-5 * 3/5 y^2
    //                                          +1 for normalize
    val c3over5 = ((SafeLong(3)<<(manW+1)) / SafeLong(5)).toBigInt.U((manW+1).W)
//     val c3over5 = math.round(3.0/5.0 * (1<<(manW+1))).toLong.U((manW+1).W)

    val (ySqExInc, ySqManW1) = multiply(ymanW1, ymanW1)
    val (ySq3over5ExInc, ySq3over5ManW1) = multiply(ySqManW1, c3over5)
    val ySq3over5Ex = yex +& yex - (exBias+1).U + ySqExInc + ySq3over5ExInc

    // ----------------------------------------------------------------------
    // 1 + 25/21 * 2^-2 * y
    //
    // here, y < 2^-4. so 1 + 25/21 * 2^-2 * y never exceeds 2.
    // We don't need to check if it exceeds 2 after addition.

    val c25over21 = ((SafeLong(25)<<manW) / SafeLong(21)).toBigInt.U((manW+1).W)
//     val c25over21 = math.round(25.0/21.0 * (1<<manW)).toLong.U((manW+1).W)

    val (y25over21ExInc, y25over21ManW1) = multiply(ymanW1, c25over21)
    val y25over21Ex       = yex + y25over21ExInc
    val y25over21ShiftVal = ((exBias + 2).U - y25over21Ex)
    val y25over21ShiftOut = y25over21ShiftVal(y25over21ShiftVal.getWidth-1, shiftOut).orR
    val y25over21Shift    = Mux(y25over21ShiftOut, Fill(shiftOut, 1.U(1.W)),
                                                   y25over21ShiftVal(shiftOut-1,0))
    val y25over21Aligned  = y25over21ManW1 >> y25over21Shift
    // assuming y25over21Aligned < 1
    val onePlus25over21yManW1 = Cat(1.U(1.W), y25over21Aligned(manW-1, 0))

    // ----------------------------------------------------------------------
    // (3/5 * y^2) * (1 + 25/21 * 2^-2 * y)
    //                ^^^^^^^^^^^^^^^^^^^^ ex == 0

    val (secondTermExInc, secondTermManW1) =
      multiply(ySq3over5ManW1, onePlus25over21yManW1)
    val secondTermEx = ySq3over5Ex + secondTermExInc

    // ----------------------------------------------------------------------
    // 1/3 * y

    val c1over3 = ((SafeLong(1)<<(manW+2)) / SafeLong(3)).toBigInt.U((manW+1).W)
//     val c1over3 = math.round(1.0/3.0 * (1<<(manW+2))).toLong

    val (yOver3ExInc, yOver3ManW1) = multiply(ymanW1, c1over3)
    val yOver3Ex = yex - 2.U + yOver3ExInc

    // ----------------------------------------------------------------------
    // 1 + 2^-2 * (1/3 * y) + 2^-5 * secondTerm < 2
    //
    // Here, y < 2^-4. so 2^-2 * (1/3 * y) < 2^-7 and secondTerm =
    // (3/5 * y^2) * (1 + 25/21 * 2^-2 * y) < 2^-7. Thus
    // 2^-2 * (1/3 * y) + 2^-5 * secondTerm < 2^-6. It never affect to the hidden bit, 1.

    val firstTermAligned  = yOver3ManW1     >> ((exBias+2-1).U - yOver3Ex)
    val secondTermAligned = secondTermManW1 >> ((exBias+5-1).U - secondTermEx)

    // here we don't add 1<<manW because it will be omitted.
    val puiseuxMan = firstTermAligned( firstTermAligned.getWidth-1, 1) +  firstTermAligned(0) +
                    secondTermAligned(secondTermAligned.getWidth-1, 1) + secondTermAligned(0)

    // XXX: Note that, we later multiply sqrt(2y) to this. But the information of
    //      will not be passed to the postprocess. We add the ex to this before
    //      postprocess.
    //      y < 2^-4, so yex < exBias-1.
    val sqrt2yExNobiasNeg = ((yex.zext - (exBias - 1).S) >> 1)(exW-1, 0)
    val puiseuxEx         = sqrt2yExNobiasNeg + exBias.U(exW.W)

    assert(puiseuxMan.getWidth == manW)

    // sqrt(2y) is multiplied in postprocess.

    val zex  = Mux(xMoreThan1, 0.U, puiseuxEx)
    val zman = Mux(xMoreThan1, 0.U, puiseuxMan)

    io.zother.zman.get := ShiftRegister(zman, nStage)
    io.zother.zex      := ShiftRegister(zex,  nStage)
  }
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
  val shiftOut = log2UpL(manW) // for FP32, log2Up(23) = 5

  val nStage = stage.total
  def getStage = nStage

  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order
  val extraBits = polySpec.extraBits

  val io = IO(new Bundle {
    val en = Input(UInt(1.W))
    // ex and some flags
    val zother = Flipped(new ACosNonTableOutput(spec, polySpec))
    // table interpolation results
    val zres   = Input(UInt(fracW.W))
    // output
    val z      = Output(UInt(spec.W.W))
  })

  val xsgn          = io.zother.xsgn.asBool
  val xLessThanHalf = io.zother.xLessThanHalf
  val znan         = io.zother.znan
  val zzero        = io.zother.zzero
  val zexNonTable  = io.zother.zex
//   printf("cir: zother.xsgn       = %d\n", xsgn      )
//   printf("cir: zother.znan       = %d\n", znan      )
//   printf("cir: zother.zIsPuiseux = %d\n", zIsPuiseux)
//   printf("cir: zexNonTable       = %d\n", zexNonTable )
//   printf("cir: zmanNonTable      = %b\n", zmanNonTable)

  if(order == 0) { // order == 0, for BF16

    val zex0  = Mux(zzero, 0.U, zexNonTable)
    val zman0 = Mux(zzero, 0.U, io.zres)

    // -----------------------------------------------------------------------
    // if x > 0, return acos(x).

    val zmanPos = Wire(UInt(manW.W))
    val zexPos  = Wire(UInt(exW.W))

    if(extraBits == 0) {
      zmanPos := zman0
      zexPos  := zex0
    } else {
      val zmanRounded = zman0(fracW-1, extraBits) +& zman0(extraBits-1)
      val zmanMoreThan2AfterRound = zmanRounded(manW)
      val zmanResult = zmanRounded(manW-1, 0)
      zmanPos := zmanResult
      zexPos  := zex0 + zmanMoreThan2AfterRound
    }

    // -----------------------------------------------------------------------
    // if x < 0, return pi - acos(|x|)

    val zmanNeg = Wire(UInt(manW.W))
    val zexNeg  = Wire(UInt(exW.W))

    val piEx    = 1
    val piFixed = Real.pi(fracW-1).toBigInt
//     val piFixed = math.round(math.Pi * (1<<(fracW-1))).toLong

    val zmanShift   = (piEx + exBias).U(exW.W) - zex0
    val zmanShift0  = zmanShift(log2Up(fracW), 0)
    val zmanShifted = Cat(1.U(1.W), zman0) >> zmanShift0
    val zmanAligned = Mux(zmanShift > fracW.U, 0.U, zmanShifted)

    val res0 = piFixed.U - zmanAligned
    val res0MoreThan2 = res0(fracW)

    val res = Mux(res0MoreThan2.asBool, res0(fracW-1, 0), Cat(res0, 0.U(1.W))(fracW-1, 0))

    if(extraBits == 0) {
      zmanNeg := res
      zexNeg  := exBias.U(exW.W) + res0MoreThan2
    } else {
      val zmanRounded = res(fracW-1, extraBits) +& res(extraBits-1)
      val zmanMoreThan2AfterRound = zmanRounded(manW)
      val zmanResult = zmanRounded(manW-1, 0)
      zmanNeg := zmanResult
      zexNeg  := exBias.U(exW.W) + res0MoreThan2 + zmanMoreThan2AfterRound
    }

    // -----------------------------------------------------------------------
    // merge them

    val zsgn = znan & xsgn.asUInt // if znan, keep sign of nan. if non-nan(znan==0), 0.
    val zex  = Mux(znan, Fill(exW, 1.U(1.W)),
               Mux(xsgn, zexNeg, zexPos))
    val zman = Mux(znan, Cat(1.U(1.W), Fill(manW-1, 0.U(1.W))),
               Mux(xsgn, zmanNeg, zmanPos))

    val z  = enable(io.en, Cat(zsgn, zex, zman))

    io.z := ShiftRegister(z, nStage)

  } else { // order != 0, for FP32
    assert(io.zother.zIsPuiseux.isDefined)
    assert(io.zother.zman.isDefined      )
    val zIsPuiseux   = io.zother.zIsPuiseux.get
    val zmanNonTable = io.zother.zman.get

    val halfPiFixed = (Real.pi / Real.two)(fracW).toBigInt.U((fracW+1).W)
//     val halfPiFixed = math.round(Pi * 0.5 * (1 << fracW)).U((fracW+1).W)
    val piman = (new RealGeneric(spec, Pi)).man.toLong.U(manW.W)
    val piex  = (new RealGeneric(spec, Pi)).ex.U(exW.W)

    val zOutOfRangeMan = Mux(xsgn, piman, 0.U(manW.W))
    val zOutOfRangeEx  = Mux(xsgn, piex,  0.U(exW.W))

    val res0 = Cat(0.U(1.W), Mux(xLessThanHalf, Cat(0.U(1.W), io.zres), Cat(io.zres, 0.U(1.W))))
    val res  = halfPiFixed + Mux(xsgn, res0, ~res0 + 1.U)

    val shift = (fracW+2).U - (res.getWidth.U - PriorityEncoder(Reverse(res)))
    val resShifted = (res << shift)(fracW+1, 1) - (BigInt(1)<<fracW).U

    val zexTable  = Mux(zzero, zOutOfRangeEx,  (exBias+1).U(exW.W) - shift)
    val zmanTable = Mux(zzero, zOutOfRangeMan, (resShifted >> extraBits) + resShifted(extraBits-1))
  //   printf("cir: zexTable       = %d\n", zexTable )
  //   printf("cir: zmanTable      = %b\n", zmanTable)

    // ---------------------------------------------------------------------------
    // Puiseux: multiply table result `sqrt(2y)` and OtherPath result. then calc
    //          (pi - z) or z. Note that if |x| > 1, the zNonTable is 0.

    val sqrt2y = Cat(1.U(1.W), io.zres) // fracW
    val puiseuxTermManW1 = Cat(1.U(1.W), zmanNonTable)
  //   printf("cir: sqrt2y           = %b\n", sqrt2y   )
  //   printf("cir: puiseuxTermManW1 = %b\n", puiseuxTermManW1)

    val puiseuxProd = sqrt2y * puiseuxTermManW1
    val puiseuxMoreThan2 = puiseuxProd(manW+1+fracW+1-1)
    val puiseuxShifted   = Mux(puiseuxMoreThan2, puiseuxProd(manW+1+fracW+1-2, fracW+1),
                                                 puiseuxProd(manW+1+fracW+1-3, fracW))
    val puiseuxRounded   = puiseuxShifted +&
                           Mux(puiseuxMoreThan2, puiseuxProd(fracW), puiseuxProd(fracW-1))
    val puiseuxMoreThan2AfterRound = puiseuxRounded(manW)

    val puiseuxManW1     = Mux(puiseuxMoreThan2AfterRound,
                               Cat(1.U(1.W), 0.U(manW.W)),
                               Cat(1.U(1.W), puiseuxRounded(manW-1, 0)))

    val puiseuxEx        = zexNonTable + puiseuxMoreThan2 + puiseuxMoreThan2AfterRound
  //   printf("cir: puiseuxEx      = %d\n", puiseuxEx   )
  //   printf("cir: puiseuxManW1   = %b\n", puiseuxManW1)

    // pi - z

    val piExNobias = 1
    val piManW1    = Real.pi(manW-piExNobias).toBigInt.U((1+manW).W)
//     val piManW1    = math.round(Pi * (1<<(manW-piExNobias))).toLong.U((1+manW).W)

    val puiseuxShiftVal = (exBias+1).U(exW.W) - puiseuxEx
    val puiseuxShiftOut = puiseuxShiftVal(puiseuxShiftVal.getWidth-1, shiftOut).orR
    val puiseuxShift    = Mux(puiseuxShiftOut, Fill(shiftOut, 1.U(1.W)), puiseuxShiftVal(shiftOut-1, 0))

    val puiseuxAligned = puiseuxManW1 >> puiseuxShift
    val puiseuxSub     = piManW1 - puiseuxAligned

    val zexPuiseux  = Mux(xsgn, (piExNobias+exBias).U(exW.W), puiseuxEx)
    val zmanPuiseux = Mux(xsgn, puiseuxSub(manW-1, 0),        puiseuxManW1(manW-1, 0))
  //   printf("cir: zexPuiseux   = %d\n", zexPuiseux )
  //   printf("cir: zmanPuiseux  = %b\n", zmanPuiseux)

    // ---------------------------------------------------------------------------
    // select

    val zsgn = znan & xsgn.asUInt // if znan, keep sign of nan. if non-nan(znan==0), 0.
    val zex  = Mux(znan, Fill(exW, 1.U(1.W)),
               Mux(zIsPuiseux, zexPuiseux,  zexTable(exW-1, 0)))
    val zman = Mux(znan, Cat(1.U(1.W), Fill(manW-1, 0.U(1.W))),
               Mux(zIsPuiseux, zmanPuiseux(manW-1, 0), zmanTable(manW-1, 0)))

    val z  = enable(io.en, Cat(zsgn, zex, zman))

    io.z := ShiftRegister(z, nStage)
  }
}

// -------------------------------------------------------------------------
//                      _     _                _
//   ___ ___  _ __ ___ | |__ (_)_ __   ___  __| |
//  / __/ _ \| '_ ` _ \| '_ \| | '_ \ / _ \/ _` |
// | (_| (_) | | | | | | |_) | | | | |  __/ (_| |
//  \___\___/|_| |_| |_|_.__/|_|_| |_|\___|\__,_|
// -------------------------------------------------------------------------

class ACosGeneric(
  val spec     : RealSpec,
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage    : MathFuncPipelineConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
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

  val sqrtCbits = SqrtTableCoeff.getCBits(spec, polySpec)
  val sqrtCalcW = SqrtTableCoeff.getCalcW(spec, polySpec)
  val acosCbits = ACosTableCoeff.getCBits(spec, polySpec)
  val acosCalcW = ACosTableCoeff.getCalcW(spec, polySpec)

  def getCbit  = Seq(sqrtCbits, acosCbits).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )
  def getCalcW = Seq(sqrtCalcW, acosCalcW).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  val cbits = getCbit
  val calcW = getCalcW

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

  val acosPre   = Module(new ACosPreProcess (spec, polySpec, stage.preStage))
  val acosTab   = Module(new ACosTableCoeff (spec, polySpec, cbits))
  val acosOther = Module(new ACosOtherPath  (spec, polySpec, stage.calcStage))
  val acosPost  = Module(new ACosPostProcess(spec, polySpec, stage.postStage))
  val sqrtTab   = if(order != 0) {Some(Module(new SqrtTableCoeff (spec, polySpec, cbits)))} else {None}

  val acosPreUseSqrtPCGapReg = if(order != 0) {
    Some(ShiftRegister(acosPre.io.useSqrt.get, pcGap))
  } else {None}

  val acosPreAdrPCGapReg     = ShiftRegister(acosPre.io.adr,     pcGap)

  acosPre.io.en        := io.en
  acosPre.io.x         := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  acosTab.io.en        := enPCGapReg && (!acosPreUseSqrtPCGapReg.getOrElse(false.B))
  acosTab.io.adr       := acosPreAdrPCGapReg
  if(order != 0) {
    sqrtTab.get.io.en  := enPCGapReg && acosPreUseSqrtPCGapReg.get
    sqrtTab.get.io.adr := acosPreAdrPCGapReg
  }
  acosOther.io.x       := xdecPCGapReg
  if(order != 0) {
    acosOther.io.useSqrt.get := acosPreUseSqrtPCGapReg.get
    acosOther.io.yex.get     := ShiftRegister(acosPre.io.yex.get,  pcGap)
    acosOther.io.yman.get    := ShiftRegister(acosPre.io.yman.get, pcGap)
  }

  // after preprocess
  // XXX Since `PCReg`s are delayed by nPreStage, the timing is the same as acosPre output.
  assert(acosPre.io.adr === 0.U               || enPCReg)
  assert(acosPre.io.dx.getOrElse(0.U) === 0.U || enPCReg)

  // table takes input from preprocess, so table output is delayed compared to acos input.
  assert(acosTab.io.cs.asUInt === 0.U || (enPCGapReg && !acosPreUseSqrtPCGapReg.getOrElse(false.B)))
  if(order != 0) {
    assert(sqrtTab.get.io.cs.asUInt === 0.U || (enPCGapReg &&  acosPreUseSqrtPCGapReg.getOrElse(false.B)))
  }

  // --------------------------------------------------------------------------

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, cbits, stage.calcStage))

  if(order != 0) {
    polynomialEval.io.dx.get := ShiftRegister(acosPre.io.dx.get, pcGap)
    polynomialEval.io.coeffs.cs := (acosTab.io.cs.cs.asUInt | sqrtTab.get.io.cs.cs.asUInt).
      asTypeOf(new MixedVec(cbits.map{w => UInt(w.W)}))
  } else {
    polynomialEval.io.coeffs.cs := acosTab.io.cs.cs
  }

  val polynomialResultCPGapReg = ShiftRegister(polynomialEval.io.result, cpGap)

  acosPost.io.en     := enCPGapReg
  acosPost.io.zother := ShiftRegister(acosOther.io.zother, cpGap)
  acosPost.io.zres   := polynomialResultCPGapReg

  io.z := acosPost.io.z
}


