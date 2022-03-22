//% @file acos.scala
//
// ACos function
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

import rial.mathfunc._

// Calculate acos(x) if x is in range [-1, 1].
// otherwise, 0 or pi for positive and negative x, respectively.
//
// for small x, pi/2 - acos(x) = x + x^3/6 + 3x^5/40 + O(x^7)
// for x < 0.5, pi/2 - acos(x) is approximated by polynomial.
// for x > 0.5, acos(x) is approximated by polynomial.
// for x close to 1, use puiseux series:
//   let y = 1-x
//   acos(y) = sqrt(2y) * (1 + y/12 + 3y^2/160 + 5y^3/896 + 35y^4/18432 + O(y^5))
//
// In case of Taylor series, the condition where pi/2 - acos(x) has enough
// precision is:
//   3x^5/40 < 2^-23
//       x^5 < 2^-20
//       x   < 2^-4.
// So if x < 2^-4, that means that x.ex < 2^-5, use Taylor series.
//
// In case of Puiseux series, the condition is:
//   35y^4/18432 < 2^-23
//     y^4/526.628... < 2^-23
//     y^4 < 2^-14
//     y   < 2^-4
// So if 1-2^-4 < x, that means that x.man(22,20).andR === 1, use Puiseux series.
//
// Normal table calculates f(x) = pi/2 - acos(x).
//   acos( x) = pi/2 - pi/2 + acos(x) = pi/2 - [pi/2 - acos(x)]
//   acos(-x) = pi - acos(x)          = pi/2 + [pi/2 - acos(x)]

// -------------------------------------------------------------------------
//  _ __  _ __ ___ _ __  _ __ ___   ___ ___  ___ ___
// | '_ \| '__/ _ \ '_ \| '__/ _ \ / __/ _ \/ __/ __|
// | |_) | | |  __/ |_) | | | (_) | (_|  __/\__ \__ \
// | .__/|_|  \___| .__/|_|  \___/ \___\___||___/___/
// |_|            |_|
// -------------------------------------------------------------------------

// To use puiseux series, we need to calculate sqrt(2y)
//
//
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
  val exAdrW    = MathFuncACosSim.calcExAdrW(spec)

  val io = IO(new Bundle {
    val en  = Input (UInt(1.W))
    val x   = Flipped(new DecomposedRealOutput(spec))
    val adr = Output(UInt((exAdrW+adrW).W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }

    // if x is close to 1, we switch to puiseux series expansion. To do that,
    // we need to calculate sqrt(2*(1-x)). This output indicates we need to
    // calc sqrt.
    val useSqrt = Output(Bool())
    val yex     = Output(UInt(spec.exW.W))
    val yman    = Output(UInt(spec.manW.W))
  })

  val xsgn = enable(io.en, io.x.sgn)
  val xex  = enable(io.en, io.x.ex )
  val xman = enable(io.en, io.x.man)

  // useSqrt condition is:
  //   0          < 1-x < 2^-4
  //   1.0 - 2^-4 < x   < 1
  // in this range, x.ex == -1.
  //
  // |x|11111110|111xxx...x|
  //            ^  ^
  //            |  +- 2^-4 bit
  //            +- hidden bit = 2^-1

  val useSqrt = enable(io.en, (xex === (exBias-1).U) && xman(manW-1, manW-3).andR)
  io.useSqrt := ShiftRegister(useSqrt, nStage)

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

  io.yex  := ShiftRegister(yex(exW-1, 0),     nStage)
  io.yman := ShiftRegister(ymanW1(manW-1, 0), nStage)

//   printf("cir: yex    = %d\n", yex.zext - exBias.S)
//   printf("cir: ymanW1 = %b\n", ymanW1)

  assert(!io.en || (ymanW1(manW) === 1.U))

  val exAdrACos = ((exBias - 1).U(exW.W) - xex)(exAdrW-1, 0)
  val adrACos   = Cat(exAdrACos, xman(manW-1, dxW))

  // we need sqrt(2y). Sqrt requires only 1 LSB to distinguish x in 1~2 and 2~4.
  // To pass yex+1, the last bit is inverted.
  val exAdrSqrt = Cat(Fill(exAdrW-1, 0.U(1.W)), ~yex(0))
  val adrSqrt   = Cat(exAdrSqrt, ymanW1(manW-1, dxW))

  val adr = enable(io.en, Mux(useSqrt, adrSqrt, adrACos))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dxACos = Cat(  ~xman(dxW-1),   xman(dxW-2, 0))
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
  val taylorOrder: Int,
) extends Module {

  val manW   = spec.manW
  val adrW   = polySpec.adrW
  val fracW  = polySpec.fracW
  val order  = polySpec.order
  val exAdrW = MathFuncACosSim.calcExAdrW(spec)

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val adr = Input  (UInt((exAdrW+adrW).W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  val exAdr = io.adr(exAdrW + adrW - 1, adrW)
  val adr   = io.adr(adrW - 1, 0)

  val taylorThreshold = MathFuncACosSim.calcTaylorThreshold(manW, taylorOrder)

  if(order == 0) {
    val tbl = VecInit( (-1 to taylorThreshold.toInt by -1).map( exponent => {
      VecInit( (0L to (1L<<adrW)-1L).map(
        n => {
          val x = scalb(1.0 + n.toDouble/(1L<<adrW), exponent.toInt)
          val y = round(acos(x)*(1L<<fracW))
          y.U((fracW+1).W)
        } ) )
      } ) )
    assert(maxCbit(0) == fracW)

    io.cs.cs(0) := enable(io.en, tbl(exAdr)(adr))

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
    taylorOrder: Int,
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
    taylorOrder: Int,
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

// for x close to 0, use Taylor series:
//   pi/2 - acos(x) = x + x^3/6 + 3x^5/40 + O(x^7)
//                    ^^^^^^^^^
//                    this part is used
//
// for x close to 1, use Puiseux series:
//   acos(1-x) = sqrt(2x) * (1 + x/12 + 3x^2/160 + 5x^3/896 + O(x^4))
//                           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
//                           this part is used
//
// for |x| > 1, set IsPuiseux = true and returns zero. In the postprocess,
// the result of Puiseux will be subtracted from Pi if x.sgn == 1.
//

class ACosNonTableOutput(val spec: RealSpec) extends Bundle {
  val xsgn       = Output(UInt(1.W))
  val znan       = Output(Bool())
  val zIsTaylor  = Output(Bool()) // needs pi/2 +/- z if xsgn
  val zIsPuiseux = Output(Bool()) // needs pi - z     if xsgn
  val zex        = Output(UInt(spec.exW.W))
  val zman       = Output(UInt(spec.manW.W))
}

class ACosOtherPath(
  val spec     : RealSpec, // Input / Output floating spec
  val polySpec : PolynomialSpec,
  val stage    : PipelineStageConfig,
  val taylorOrder: Int
) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val x       = Flipped(new DecomposedRealOutput(spec))
    val yex     = Input(UInt(exW.W))
    val yman    = Input(UInt(manW.W))
    val useSqrt = Input(Bool()) // == use Puiseux series
    val zother  = new ACosNonTableOutput(spec)
  })

  // --------------------------------------------------------------------------
  // Special value (nan, |x| larger than 1)

  val znan  = io.x.nan
  io.zother.znan := ShiftRegister(znan, nStage)
  io.zother.xsgn := ShiftRegister(io.x.sgn, nStage)

  // if |x| > 1, return 0 as puiseux result.
  val xMoreThan1 = exBias.U <= io.x.ex

  io.zother.zIsPuiseux := ShiftRegister(xMoreThan1 || io.useSqrt, nStage)

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
  // Taylor
  //
  // pi/2 - acos(x) = x + x^3/6 + O(x^5)
  //                = x(1 + x^2/6) if x << 1

  val taylorThreshold = MathFuncACosSim.calcTaylorThreshold(manW, 3)
  val isTaylor        = io.x.ex < (exBias + taylorThreshold).U
  val isConstant      = io.x.ex < (exBias - manW).U

  io.zother.zIsTaylor := ShiftRegister(isTaylor, nStage)

  val zTaylorEx  = Wire(UInt(exW.W))
  val zTaylorMan = Wire(UInt(manW.W))

  if(taylorOrder < 3) {
    zTaylorEx  := Mux(isConstant, 0.U(exW.W),  io.x.ex)
    zTaylorMan := Mux(isConstant, 0.U(manW.W), io.x.man)
  } else {
    assert(taylorOrder < 5)
    val xex    = io.x.ex
    val xmanW1 = Cat(1.U(1.W), io.x.man)
    val (xSqExInc, xSqManW1) = multiply(xmanW1, xmanW1)

    val c1over6Ex    = (exBias-3).U(exW.W)
    val c1over6ManW1 = math.round(1.0/6.0 * (1L<<(manW+3))).U((manW+1).W)
    val (xSq6thExInc, xSq6thManW1) = multiply(xSqManW1, c1over6ManW1)

    val xSq6thEx    = Cat(xex, 0.U(1.W)) - (3 + exBias).U((exW+1).W) + xSqExInc + xSq6thExInc
    val xSq6thShiftVal = (exBias.U - xSq6thEx)
    val xSq6thShiftOut = xSq6thShiftVal(xSq6thShiftVal.getWidth-1, shiftOut).orR
    val xSq6thShift    = Mux(xSq6thShiftOut, Fill(shiftOut, 1.U(1.W)),
                                             (exBias.U - xSq6thEx)(shiftOut-1, 0))

    val xSq6thAligned      = xSq6thManW1 >> xSq6thShift
    val xSq6thPlusOneManW1 = Cat(1.U(1.W), xSq6thAligned(manW-1, 0)) // assuming xSq6thAligned < 1
    assert(xSq6thPlusOneManW1.getWidth == manW+1)

    // in the postprocess, the result will be added/subtracted to pi/2.
    val (taylorExInc, taylorManW1) = multiply(xmanW1, xSq6thPlusOneManW1)
    val taylorEx = xex + taylorExInc

    zTaylorEx  := Mux(isConstant, 0.U(exW.W),  taylorEx)
    zTaylorMan := Mux(isConstant, 0.U(manW.W), taylorManW1(manW-1, 0))
  }

  // --------------------------------------------------------------------------
  // Puiseux
  //
  // let y = 1 - x
  // acos(1-y) = sqrt(2y) * (1 + y/12 + 3y^2/160 + 5y^3/896 + O(y^4))
  //           = sqrt(2y) * ((1 + 1/3 * 2^-2 * y) + 3y^2/160 * (1 + 25/21 * 2^-2 * y))
  //                        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  //                        calc this part

  val yex    = io.yex
  val ymanW1 = Cat(1.U(1.W), io.yman)

  // ----------------------------------------------------------------------
  // 2^-5 * 3/5 y^2
  //                                          +1 for normalize
  val c3over5 = math.round(3.0/5.0 * (1<<(manW+1))).toLong.U((manW+1).W)

  val (ySqExInc, ySqManW1) = multiply(ymanW1, ymanW1)
  val (ySq3over5ExInc, ySq3over5ManW1) = multiply(ySqManW1, c3over5)
  val ySq3over5Ex = yex +& yex - (exBias+1).U + ySqExInc + ySq3over5ExInc

  // ----------------------------------------------------------------------
  // 1 + 25/21 * 2^-2 * y
  //
  // here, y < 2^-4. so 1 + 25/21 * 2^-2 * y never exceeds 2.
  // We don't need to check if it exceeds 2 after addition.

  val c25over21 = math.round(25.0/21.0 * (1<<manW)).toLong

  val (y25over21ExInc, y25over21ManW1) = multiply(ymanW1, c25over21.U((manW+1).W))
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

  val c1over3 = math.round(1.0/3.0 * (1<<(manW+2))).toLong

  val (yOver3ExInc, yOver3ManW1) = multiply(ymanW1, c1over3.U((manW+1).W))
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

  val zPuiseuxEx  = Mux(xMoreThan1, 0.U, puiseuxEx)
  val zPuiseuxMan = Mux(xMoreThan1, 0.U, puiseuxMan)

  // --------------------------------------------------------------------------
  // select Taylor/Puiseux

  val zex  = Mux(isTaylor, zTaylorEx,  zPuiseuxEx)
  val zman = Mux(isTaylor, zTaylorMan, zPuiseuxMan)

  io.zother.zman := ShiftRegister(zman, nStage)
  io.zother.zex  := ShiftRegister(zex,  nStage)
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
    val zother = Flipped(new ACosNonTableOutput(spec))
    // table interpolation results
    val zres   = Input(UInt(fracW.W))
    // output
    val z      = Output(UInt(spec.W.W))
  })

  val xsgn         = io.zother.xsgn
  val zIsTaylor    = io.zother.zIsTaylor
  val zIsPuiseux   = io.zother.zIsPuiseux
  val znan         = io.zother.znan
  val zexNonTable  = io.zother.zex
  val zmanNonTable = io.zother.zman
//   printf("cir: zother.xsgn       = %d\n", xsgn      )
//   printf("cir: zother.znan       = %d\n", znan      )
//   printf("cir: zother.zIsTaylor  = %d\n", zIsTaylor )
//   printf("cir: zother.zIsPuiseux = %d\n", zIsPuiseux)
//   printf("cir: zexNonTable       = %d\n", zexNonTable )
//   printf("cir: zmanNonTable      = %b\n", zmanNonTable)

  // ---------------------------------------------------------------------------
  // taylor and table (pi/2 +/- z). Note that x < 2^-23, the zNonTable is 0

  val halfPiFixed = math.round(Pi * 0.5 * (1 << fracW)).U((fracW+1).W)

  val taylorShiftVal = exBias.U(exW.W) - zexNonTable
  val taylorShiftOut = taylorShiftVal(taylorShiftVal.getWidth-1, shiftOut).orR
  val taylorShift    = Mux(taylorShiftOut, Fill(shiftOut, 1.U(1.W)), taylorShiftVal(shiftOut-1, 0))
  val taylorAligned  = Cat(Cat(1.U(1.W), zmanNonTable), 0.U((fracW-manW).W)) >> taylorShift
  assert(taylorAligned.getWidth == fracW+1)

//   printf("cir: halfPiFixed       = %b\n", halfPiFixed  )
//   printf("cir: taylorAligned     = %b\n", taylorAligned)

  val res0 = Mux(zIsTaylor, taylorAligned, Cat(io.zres, 0.U(1.W)))
  val res  = halfPiFixed + Mux(xsgn.asBool, Cat(0.U(1.W), res0), Cat(1.U(1.W), ~res0 + 1.U))

  val shift = (fracW+2).U - (res.getWidth.U - PriorityEncoder(Reverse(res)))
  val resShifted = (res << shift)(fracW+1, 1) - (1<<fracW).U

  val zexTable  = (exBias+1).U(exW.W) - shift
  val zmanTable = (resShifted >> extraBits) + resShifted(extraBits-1)
//   printf("cir: zexTable       = %d\n", zexTable )
//   printf("cir: zmanTable      = %b\n", zmanTable)

  // ---------------------------------------------------------------------------
  // Puiseux: (pi - z) or z. Note that if |x| > 1, the zNonTable is 0.

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
  val piManW1    = math.round(Pi * (1<<(manW-piExNobias))).toLong.U((1+manW).W)

  val puiseuxShiftVal = (exBias+1).U(exW.W) - puiseuxEx
  val puiseuxShiftOut = puiseuxShiftVal(puiseuxShiftVal.getWidth-1, shiftOut).orR
  val puiseuxShift    = Mux(puiseuxShiftOut, Fill(shiftOut, 1.U(1.W)), puiseuxShiftVal(shiftOut-1, 0))

  val puiseuxAligned = puiseuxManW1 >> puiseuxShift
  val puiseuxSub     = piManW1 - puiseuxAligned

  val zexPuiseux  = Mux(xsgn.asBool, (piExNobias+exBias).U(exW.W), puiseuxEx)
  val zmanPuiseux = Mux(xsgn.asBool, puiseuxSub(manW-1, 0),        puiseuxManW1(manW-1, 0))
//   printf("cir: zexPuiseux   = %d\n", zexPuiseux )
//   printf("cir: zmanPuiseux  = %b\n", zmanPuiseux)

  // ---------------------------------------------------------------------------
  // select

  val zsgn = znan & xsgn // if znan, keep sign of nan. if non-nan(znan==0), 0.
  val zex  = Mux(znan, Fill(exW, 1.U(1.W)),
             Mux(zIsPuiseux, zexPuiseux,  zexTable(exW-1, 0)))
  val zman = Mux(znan, Cat(1.U(1.W), Fill(manW-1, 0.U(1.W))),
             Mux(zIsPuiseux, zmanPuiseux(manW-1, 0), zmanTable(manW-1, 0)))

  val z  = enable(io.en, Cat(zsgn, zex, zman))

  io.z := ShiftRegister(z, nStage)
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
  val taylorOrder: Int,
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
  val acosCbits = ACosTableCoeff.getCBits(spec, polySpec, taylorOrder)
  val acosCalcW = ACosTableCoeff.getCalcW(spec, polySpec, taylorOrder)

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
  val acosTab   = Module(new ACosTableCoeff (spec, polySpec, cbits, taylorOrder))
  val sqrtTab   = Module(new SqrtTableCoeff (spec, polySpec, cbits))
  val acosOther = Module(new ACosOtherPath  (spec, polySpec, stage.calcStage, taylorOrder))
  val acosPost  = Module(new ACosPostProcess(spec, polySpec, stage.postStage))

  val acosPreUseSqrtPCGapReg = ShiftRegister(acosPre.io.useSqrt, pcGap)
  val acosPreAdrPCGapReg     = ShiftRegister(acosPre.io.adr,     pcGap)

  acosPre.io.en        := io.en
  acosPre.io.x         := xdecomp.io.decomp
  // ------ Preprocess-Calculate ------
  acosTab.io.en        := io.en && (!acosPreUseSqrtPCGapReg)
  acosTab.io.adr       := acosPreAdrPCGapReg
  sqrtTab.io.en        := io.en || acosPreUseSqrtPCGapReg
  sqrtTab.io.adr       := acosPreAdrPCGapReg
  acosOther.io.x       := xdecPCGapReg
  acosOther.io.useSqrt := acosPreUseSqrtPCGapReg
  acosOther.io.yex     := ShiftRegister(acosPre.io.yex,  pcGap)
  acosOther.io.yman    := ShiftRegister(acosPre.io.yman, pcGap)

  // after preprocess
  // XXX Since `PCReg`s are delayed by nPreStage, the timing is the same as acosPre output.
  assert(acosPre.io.adr === 0.U               || enPCReg)
  assert(acosPre.io.dx.getOrElse(0.U) === 0.U || enPCReg)

  // table takes input from preprocess, so table output is delayed compared to acos input.
  assert(acosTab.io.cs.asUInt === 0.U || (enPCGapReg && !acosPreUseSqrtPCGapReg))
  assert(sqrtTab.io.cs.asUInt === 0.U || (enPCGapReg && !acosPreUseSqrtPCGapReg))

  // --------------------------------------------------------------------------

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, cbits, stage.calcStage))

  if(order != 0) {
    polynomialEval.io.dx.get := ShiftRegister(acosPre.io.dx.get, pcGap)
  }
  polynomialEval.io.coeffs.cs := (acosTab.io.cs.cs.asUInt | sqrtTab.io.cs.cs.asUInt).
    asTypeOf(new MixedVec(cbits.map{w => UInt(w.W)}))

  val polynomialResultCPGapReg = ShiftRegister(polynomialEval.io.result, cpGap)

  acosPost.io.en     := enCPGapReg
  acosPost.io.zother := ShiftRegister(acosOther.io.zother, cpGap)
  acosPost.io.zres   := polynomialResultCPGapReg

  io.z := acosPost.io.z
}


