//% @file atan2.scala
//
// ATan2 function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

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

// atan2(y, x), (floating, floating) input, floating output
//
// atan2(y, x) = atan(y/x) +/- (pi/0). we always need to calculate y/x. so we
// require reciprocal table to calculate 1.m/1.m and also requires atan table
// to calculate atan(y/x).
//
// To reduce the size of atan table, we use atan(1/x) = pi/2 - atan(x).
//
class ATan2Generic(
  val spec : RealSpec, // Input / Output floating spec
  val nOrderRec  : Int, val adrWRec  : Int, val extraBitsRec  : Int, // polynomial spec
  val nOrderATan : Int, val adrWATan : Int, val extraBitsATan : Int,
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val nStage = stage.total

  def getParam() = { (spec.exW, spec.manW,
    nOrderRec, adrWRec, extraBitsRec, nOrderATan, adrWATan, extraBitsATan,
    nStage, enableRangeCheck, enablePolynomialRounding) }

  def getStage() = nStage

  val manW  = spec.manW
  val exW   = spec.exW
  val exBias = spec.exBias

  val io = IO(iodef = new Bundle {
    val x   = Input (UInt(spec.W.W))
    val y   = Input (UInt(spec.W.W))
    val z   = Output(UInt(spec.W.W))
  })

  if(adrWRec  >= manW && nOrderRec  != 0) {
    println("WARNING: reciprocal table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
  }
  if(adrWATan >= manW && nOrderATan != 0) {
    println("WARNING: atan table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
  }

  val orderRec  = if (adrWRec  >= manW) { 0 } else { nOrderRec  }
  val orderATan = if (adrWATan >= manW) { 0 } else { nOrderATan }

  val bpRec  = manW+extraBitsRec
  val bpATan = manW+extraBitsATan

  val (xsgn,  xex,  xman) = FloatChiselUtil.decompose(spec, io.x)
  val (ysgn,  yex,  yman) = FloatChiselUtil.decompose(spec, io.y)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(spec, io.x)
  val (yzero, yinf, ynan) = FloatChiselUtil.checkValue(spec, io.y)
  val xExNobias = xex.zext() - exBias.S((exW+1).W)
  val yExNobias = yex.zext() - exBias.S((exW+1).W)

//   printf("x = %d|%d(%d)|%b\n", xsgn, xex, xExNobias, xman)
//   printf("y = %d|%d(%d)|%b\n", ysgn, yex, yExNobias, yman)

  val znan = (xnan || ynan) || (xzero && yzero)

  val pi_1over4 = new RealGeneric(spec, Pi * 0.25)
  val pi_over_2 = new RealGeneric(spec, Pi * 0.50)
  val pi_3over4 = new RealGeneric(spec, Pi * 0.75)
  val pi        = new RealGeneric(spec, Pi)

  val zzero     = ((xinf && !yinf) || (!xzero && yzero)) && (xsgn === 0.U(1.W))
  val zpi       = ((xinf && !yinf) || (!xzero && yzero)) && (xsgn === 1.U(1.W))
  val zhalfpi   = (!xinf && yinf) || (xzero && !yzero)
  val z1piover4 = (xinf && yinf && (xsgn === 0.U(1.W)))
  val z3piover4 = (xinf && yinf && (xsgn === 1.U(1.W)))
  val isSpecialValue = znan || zzero || zpi || zhalfpi || z1piover4 || z3piover4
  assert(znan ^ zzero ^ zpi ^ zhalfpi ^ z1piover4 ^ z3piover4 ^ !isSpecialValue)

  val zSgn = ysgn

  // ==========================================================================
  // calc y / x and x / y
//   printf("calculating y/x and x/y ...\n")

  val yOverXSgn = xsgn ^ ysgn
  val yOverXEx0 = yExNobias -& xExNobias
  val xOverYEx0 = xExNobias -& yExNobias

  // calc 1.m/1.m \in (0.5, 2.0)

  val one_over_x_manW1 = Wire(UInt((manW+1+extraBitsRec).W))
  val one_over_y_manW1 = Wire(UInt((manW+1+extraBitsRec).W))

  if(orderRec == 0) {
    val adrX = xman
    val adrY = yman

    // 1.0 / 1.m -> (0.5, 1]
    val tbl = VecInit( (0L to (1L<<adrWRec)-1L).map(
      n => {
        val x = 1.0 + n.toDouble / (1L<<adrWRec)
        val y = math.round((1.0 / x) * (1L<<manW))
        assert(y <= (1L << manW))
        y.U((manW+1).W)
      }
    ) )

    one_over_x_manW1 := tbl(adrX)
    one_over_y_manW1 := tbl(adrY)

  } else {

    val adrX = ~xman(manW-1, manW-adrWRec)
    val adrY = ~yman(manW-1, manW-adrWRec)
    val dX  = Cat(xman(manW-adrWRec-1),~xman(manW-adrWRec-2,0)).asSInt
    val dY  = Cat(yman(manW-adrWRec-1),~yman(manW-adrWRec-2,0)).asSInt

    // Create 1/x table for inverted x
    val eps = pow(2.0, -manW)
    val tableD = new FuncTableDouble( x => 1.0/(2.0-(x+eps)), orderRec )
    tableD.addRange(0.0, 1.0, 1<<adrWRec)
    val tableI = new FuncTableInt( tableD, bpRec )
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(1)

    val coeffXAll = coeffTable(adrX)
    val coeffYAll = coeffTable(adrY)
    val coeffX    = getSlices(coeffXAll, coeffWidth)
    val coeffY    = getSlices(coeffYAll, coeffWidth)
    val coeffXS   = coeffX.map( x => x.zext ) // all positive (x is inverted)
    val coeffYS   = coeffY.map( x => x.zext )

    def hornerC( c: SInt, z: SInt, dx: SInt, enableRounding: Boolean = false ) : SInt = {
      val zw   = z.getWidth
      val dxBp = dx.getWidth-1
      val dxw  = min(zw, dxBp)
      val dxl  = dx(dxBp, dxBp-dxw).asSInt
      val prod = dxl * z
      val prod_sft = dropLSB(dxw, prod)
      val sum = c +& prod_sft // Extend for safe
      sum
    }
    def polynomialEvaluationC( c : Seq[SInt], dx : SInt, enableRounding: Boolean = false ) = {
      val res = c.init.foldRight(c.last)( (c,z) => hornerC( c, z, dx, enableRounding ) )
      res
    }
    val sX0 = polynomialEvaluationC( coeffXS, dX, enablePolynomialRounding ).asUInt
    val sY0 = polynomialEvaluationC( coeffYS, dY, enablePolynomialRounding ).asUInt

    val oneOverXOvf = sX0 > (1<<bpRec).U
    val oneOverYOvf = sY0 > (1<<bpRec).U
    val oneOverXUdf = sX0 < (1<<(bpRec-1)).U
    val oneOverYUdf = sY0 < (1<<(bpRec-1)).U

    one_over_x_manW1 := Mux(xman === 0.U, (1 << bpRec).U((bpRec + 1).W),
      MuxCase(sX0, Seq(
        oneOverXOvf -> (1 <<  bpRec   ).U((bpRec + 1).W),
        oneOverXUdf -> (1 << (bpRec-1)).U((bpRec + 1).W)
      )))
    one_over_y_manW1 := Mux(yman === 0.U, (1 << bpRec).U((bpRec + 1).W),
      MuxCase(sY0, Seq(
        oneOverYOvf -> (1 <<  bpRec   ).U((bpRec + 1).W),
        oneOverYUdf -> (1 << (bpRec-1)).U((bpRec + 1).W)
      )))
  }
  val xmanW1 = 1.U(1.W) ## xman
  val ymanW1 = 1.U(1.W) ## yman

  val yOverXMan0 = ymanW1 * one_over_x_manW1
  val xOverYMan0 = xmanW1 * one_over_y_manW1

  assert(yOverXMan0.getWidth == bpRec+1+manW+1)
  assert(xOverYMan0.getWidth == bpRec+1+manW+1)
  assert(yOverXMan0(bpRec+1+manW+1-1) === 0.U(1.W)) // < 2
  assert(xOverYMan0(bpRec+1+manW+1-1) === 0.U(1.W))

  // max 1.m / 1.m = 1.111/1.000 < 2
  // min 1.m / 1.m = 1.000/1.111 > 0.5
  //
  //   1.xxx    : W = bp+1
  // * 1.yyy    : W = manW+1
  //  zz.zzzzzz : W = bp+1+manW+1
  //  ^
  //  bp+1+manW+1-1
  val bpDiv           = bpRec + manW
  assert(yOverXMan0.getWidth == bpRec+1 + manW+1)
  assert(xOverYMan0.getWidth == bpRec+1 + manW+1)
  assert(yOverXMan0(bpRec+1 + manW+1 - 1) === 0.U(1.W))
  assert(xOverYMan0(bpRec+1 + manW+1 - 1) === 0.U(1.W))
  val yOverXLessThan1 = yOverXMan0(bpRec+1 + manW+1 - 2) === 0.U(1.W)
  val xOverYLessThan1 = xOverYMan0(bpRec+1 + manW+1 - 2) === 0.U(1.W)
  val roundBits   = bpDiv - manW
  val yOverXShift = Mux(yOverXLessThan1, yOverXMan0(bpDiv-1, roundBits-1), yOverXMan0(bpDiv, roundBits))
  val xOverYShift = Mux(xOverYLessThan1, xOverYMan0(bpDiv-1, roundBits-1), xOverYMan0(bpDiv, roundBits))
  val yOverXRound = Mux(yOverXLessThan1, yOverXMan0(roundBits-2), yOverXMan0(roundBits-1))
  val xOverYRound = Mux(xOverYLessThan1, xOverYMan0(roundBits-2), xOverYMan0(roundBits-1))

  val yOverXManW1 = yOverXShift + yOverXRound
  val xOverYManW1 = xOverYShift + xOverYRound
  val yOverXMan   = yOverXManW1(manW-1, 0);
  val xOverYMan   = xOverYManW1(manW-1, 0);
  assert(yOverXManW1.getWidth == manW+1)
  assert(xOverYManW1.getWidth == manW+1)
  assert(yOverXManW1(manW) === 1.U(1.W))
  assert(xOverYManW1(manW) === 1.U(1.W))

  val yOverXEx = Mux(yOverXLessThan1, yOverXEx0 - 1.S, yOverXEx0)
  val xOverYEx = Mux(xOverYLessThan1, xOverYEx0 - 1.S, xOverYEx0)

//   printf("y/x = %d|%d|%b\n", yOverXSgn, yOverXEx, yOverXMan)
//   printf("x/y = %d|%d|%b\n", yOverXSgn, xOverYEx, xOverYMan)

  // ==========================================================================
  // constant (2^manW < x -> z = pi/2)
//   printf("calculating constant case ...\n")

  val isConstant = manW.S < yOverXEx

  val zConstantEx  = pi_over_2.ex .toLong.U(exW.W)
  val zConstantMan = pi_over_2.man.toLong.U(manW.W)

  // ==========================================================================
  // linear (x^3 / 6 < x * 2^-manW)
//   printf("calculating linear case ...\n")

  val linearThreshold = -math.round(math.ceil(manW / 2.0 + 1.0)) // -13, if FP32
  val isYOverXLinear  = !isConstant && yOverXEx < linearThreshold.S
  val isXOverYLinear  = !isConstant && xOverYEx < linearThreshold.S

  // if y/x < linearThreshold, then z = y/x. if y/x is too small, z = 0.

  val yOverXUdf = yOverXEx < -exBias.S
  val yOverXExBiased = (yOverXEx + exBias.S)(exW-1, 0)
  val zYOverXLinearEx  = Mux(yOverXUdf, 0.U(exW.W),  yOverXExBiased)
  val zYOverXLinearMan = Mux(yOverXUdf, 0.U(manW.W), yOverXMan)

  // if x/y > -linearThreshold, then z = pi/2 - x/y

  // always anchor pi/2, 2^0*1.57... = 2^0*1.1001...
  // normally linearThreshold is smaller than -2, so it does not change the exponent

  val log2ManW  = log2Up(manW) // = 5 (FP32)
  val halfPiMan = 1.U(1.W) ## pi_over_2.man.toLong.U(manW.W) ## 0.U(3.W)

  val xOverYLinearShift   = ~xOverYEx(log2ManW-1, 0) - 2.U
  val xOverYLinearAligned = xOverYManW1 >> xOverYLinearShift
  val zXOverYLinearMan0   = halfPiMan - xOverYLinearAligned

  val zXOverYLinearEx  = exBias.U(exW.W)
  val zXOverYLinearMan = zXOverYLinearMan0(manW+3-1, 3) + zXOverYLinearMan0(2)

  // ==========================================================================
  // use table (2^linearThreshold <= x < 2^-linearThreshold)
//   printf("calculating table case ...\n")

  val isYOverXTable = (linearThreshold.S <= yOverXEx) && (yOverXEx < 0.S) && !isSpecialValue
  val isXOverYTable = (linearThreshold.S <= xOverYEx) && (xOverYEx < 0.S) && !isSpecialValue

  val exAdrW   = log2Up(abs(linearThreshold).toInt)
  val exAdrMax = -linearThreshold - 1
  assert(exAdrMax > 0)
  val exAdrYOverX0 = (-yOverXEx - 1.S)(exAdrW-1, 0)
  val exAdrXOverY0 = (-xOverYEx - 1.S)(exAdrW-1, 0)
  val exAdrYOverX  = Mux(exAdrYOverX0 > exAdrMax.U, 0.U, exAdrYOverX0)
  val exAdrXOverY  = Mux(exAdrXOverY0 > exAdrMax.U, 0.U, exAdrXOverY0)

  val calcW = manW + extraBitsATan
  val atanYOverXEx    = Wire(SInt((exW+1).W))
  val atanXOverYEx    = Wire(SInt((exW+1).W))
  val atanYOverXManW1 = Wire(UInt((manW+1).W))
  val atanXOverYManW1 = Wire(UInt((manW+1).W))

  if(orderATan == 0) {
    val adrYOverX = yOverXMan
    val adrXOverY = xOverYMan

    val tbl = VecInit( (-1 to linearThreshold.toInt ).map( exponent => {
      VecInit((0L to (1L<<adrWATan)-1L).map( n => {
        // atan(x) < x.
        val x = scalb(1.0 + n.toDouble / (1L<<adrWATan), exponent)
        val y = math.round(scalb(atan(x), -exponent-1) * (1L<<calcW))
        assert(y < (1L << calcW))
        y.U((calcW+1).W)
      } ) )
    } ) )

    val atanYOverXScaled = tbl(exAdrYOverX)(adrYOverX)
    val atanXOverYScaled = tbl(exAdrXOverY)(adrXOverY)
    val atanYOverXScaledMoreThan1 = atanYOverXScaled(calcW) === 1.U
    val atanXOverYScaledMoreThan1 = atanXOverYScaled(calcW) === 1.U

    atanYOverXEx    := Mux(atanYOverXScaledMoreThan1, yOverXEx, yOverXEx - 1.S)
    atanXOverYEx    := Mux(atanXOverYScaledMoreThan1, xOverYEx, xOverYEx - 1.S)

    if(extraBitsATan > 0) {
      atanYOverXManW1 := Mux(atanYOverXScaledMoreThan1, atanYOverXScaled(calcW, extraBitsATan), atanYOverXScaled(calcW-1, extraBitsATan-1))
      atanXOverYManW1 := Mux(atanXOverYScaledMoreThan1, atanXOverYScaled(calcW, extraBitsATan), atanXOverYScaled(calcW-1, extraBitsATan-1))
    } else {
      atanYOverXManW1 := Mux(atanYOverXScaledMoreThan1, atanYOverXScaled(calcW, 0), Cat(atanYOverXScaled(calcW-1, 0), 0.U(1.W)))
      atanXOverYManW1 := Mux(atanXOverYScaledMoreThan1, atanXOverYScaled(calcW, 0), Cat(atanXOverYScaled(calcW-1, 0), 0.U(1.W)))
    }
  } else {

    val dxbp = manW - adrWATan - 1
    val adrYOverX = yOverXMan(manW-1, manW-adrWATan)
    val adrXOverY = xOverYMan(manW-1, manW-adrWATan)
    val dYOverX   = Cat(~yOverXMan(manW-adrWATan-1), yOverXMan(manW-adrWATan-2,0)).asSInt
    val dXOverY   = Cat(~xOverYMan(manW-adrWATan-1), xOverYMan(manW-adrWATan-2,0)).asSInt

    val coeffWidth = (-1 to linearThreshold.toInt by -1).map(exponent => {
      val tableD = new FuncTableDouble( x => scalb(atan(scalb(1.0 + x, exponent)), -exponent-1), orderATan )
      tableD.addRange(0.0, 1.0, 1<<adrWATan)
      val tableI = new FuncTableInt( tableD, calcW ) // convert float table into int
      val w = tableI.getCBitWidth(/*sign mode = */0)
      w
    }).reduce( (lhs, rhs) => {
      lhs.zip(rhs).map( x => max(x._1, x._2))
    })
    val coeffTables = VecInit((-1 to linearThreshold.toInt by -1).map(exponent => {
      val tableD = new FuncTableDouble( x => scalb(atan(scalb(1.0 + x, exponent)), -exponent-1), orderATan )
      tableD.addRange(0.0, 1.0, 1<<adrWATan)
      val tableI = new FuncTableInt( tableD, calcW ) // convert float table into int
      tableI.getVectorWithWidth(coeffWidth, /*sign mode = */0)
    }))

    val coeffYOverXTable = coeffTables(exAdrYOverX)
    val coeffXOverYTable = coeffTables(exAdrXOverY)
    val coeffYOverX      = getSlices(coeffYOverXTable(adrYOverX), coeffWidth)
    val coeffXOverY      = getSlices(coeffXOverYTable(adrXOverY), coeffWidth)
    val coeffSYOverX     = coeffYOverX.map( x => x.asSInt )
    val coeffSXOverY     = coeffXOverY.map( x => x.asSInt )

    def hornerC( c: SInt, z: SInt, dx: SInt, enableRounding: Boolean = false ) : SInt = {
      val zw   = z.getWidth
      val dxBp = dx.getWidth-1
      val dxw  = min(zw, dxBp)
      val dxl  = dx(dxBp, dxBp-dxw).asSInt
      val prod = dxl * z
      val prod_sft = dropLSB(dxw, prod)
      val sum = c +& prod_sft // Extend for safe
      sum
    }
    def polynomialEvaluationC( c : Seq[SInt], dx : SInt, enableRounding: Boolean = false ) = {
      val res = c.init.foldRight(c.last)( (c,z) => hornerC( c, z, dx, enableRounding ) )
      res
    }
    val s0YOverX = polynomialEvaluationC( coeffSYOverX, dYOverX, enablePolynomialRounding ).asUInt
    val s0XOverY = polynomialEvaluationC( coeffSXOverY, dXOverY, enablePolynomialRounding ).asUInt

    val log2CalcW    = log2Up(calcW)
    val s0YOverXLenS = s0YOverX.getWidth.S - PriorityEncoder(Reverse(s0YOverX)).zext()
    val s0XOverYLenS = s0XOverY.getWidth.S - PriorityEncoder(Reverse(s0XOverY)).zext()
    val s0YOverXLen  = s0YOverXLenS(log2CalcW-1, 0)
    val s0XOverYLen  = s0XOverYLenS(log2CalcW-1, 0)
    assert(0.S <= s0YOverXLenS)
    assert(0.S <= s0XOverYLenS)

    val s0YOverXShift = (calcW+1).U(log2CalcW.W) - s0YOverXLen
    val s0XOverYShift = (calcW+1).U(log2CalcW.W) - s0XOverYLen
    assert(0.U <= s0YOverXShift && s0YOverXShift <= calcW.U)
    assert(0.U <= s0XOverYShift && s0XOverYShift <= calcW.U)

    val resYOverXEx  = yOverXEx + 1.S - s0YOverXShift.zext
    val resXOverYEx  = xOverYEx + 1.S - s0XOverYShift.zext
    val resYOverXMan = (s0YOverX(calcW, 0) << s0YOverXShift)(calcW, 0)
    val resXOverYMan = (s0XOverY(calcW, 0) << s0XOverYShift)(calcW, 0)

    // ------------------------------------------------------------------------
    // shift/round atan(y/x) and return it

    atanYOverXEx := resYOverXEx
    if(extraBitsATan > 0) {
      atanYOverXManW1 := resYOverXMan(manW+extraBitsATan, extraBitsATan) + resYOverXMan(extraBitsATan-1)
    } else {
      atanYOverXManW1 := resYOverXMan
    }

    // ------------------------------------------------------------------------
    // calc pi/2 - atan(x/y) = atan(y/x)

    val halfPiMan = ((pi_over_2.man.toLong + (1L<<manW)) << 3).U((manW+1+3).W)

    val resXOverYAlign0 = Wire(UInt((manW+1+3).W))
    if(extraBitsATan >= 3) {
        resXOverYAlign0 := (resXOverYMan >> (extraBitsATan - 3))(manW+1+3-1, 0)
    } else {
        resXOverYAlign0 := resXOverYMan << (3 - extraBitsATan)
    }
    val log2ManW3 = log2Up(manW+1+3)
    val resXOverYShift = -resXOverYEx(log2ManW3-1, 0)
    val resXOverYAligned = resXOverYAlign0 >> resXOverYShift

    val halfPiMinusAtanXOverYMan = halfPiMan - resXOverYAligned
    val halfPiMinusAtanXOverYLessThan1 = halfPiMinusAtanXOverYMan(manW+3) === 0.U(1.W)

    atanXOverYEx    := Mux(halfPiMinusAtanXOverYLessThan1, -1.S((exW+1).W), 0.S((exW+1).W))
    atanXOverYManW1 := Mux(halfPiMinusAtanXOverYLessThan1,
      halfPiMinusAtanXOverYMan(manW+2, 2) + halfPiMinusAtanXOverYMan(1),
      halfPiMinusAtanXOverYMan(manW+3, 3) + halfPiMinusAtanXOverYMan(2))
  }

  val zYOverXTableEx0 = atanYOverXEx + exBias.S((exW+1).W)
  val zXOverYTableEx0 = atanXOverYEx + exBias.S((exW+1).W)
  when(isYOverXTable) {assert(zYOverXTableEx0 > 0.S)}
  when(isXOverYTable) {assert(zXOverYTableEx0 > 0.S)}
  val zYOverXTableEx  = zYOverXTableEx0(exW-1, 0)
  val zXOverYTableEx  = zXOverYTableEx0(exW-1, 0)

  when(isYOverXTable) {assert(atanYOverXManW1(manW) === 1.U(1.W))}
  when(isXOverYTable) {assert(atanXOverYManW1(manW) === 1.U(1.W))}
  val zYOverXTableMan = atanYOverXManW1(manW-1, 0)
  val zXOverYTableMan = atanXOverYManW1(manW-1, 0)

  // ==========================================================================
  // select correct one
//   printf("selecting ...\n")

  assert(isYOverXLinear ^ isYOverXTable ^ isXOverYTable ^ isXOverYLinear ^ isConstant ^ znan)
//   printf("isYOverXLinear = %b, zYOverXLinearMan = %b\n", isYOverXLinear, zYOverXLinearMan)
//   printf("isXOverYLinear = %b, zXOverYLinearMan = %b\n", isXOverYLinear, zXOverYLinearMan)
//   printf("isYOverXTable  = %b, zYOverXTableMan  = %b\n", isYOverXTable , zYOverXTableMan )
//   printf("isXOverYTable  = %b, zXOverYTableMan  = %b\n", isXOverYTable , zXOverYTableMan )
//   printf("isConstant     = %b, zConstantMan     = %b\n", isConstant    , zConstantMan    )

  val zMan0 = MuxCase(zYOverXLinearMan, Seq(
    isYOverXTable  -> zYOverXTableMan,
    isXOverYTable  -> zXOverYTableMan,
    isXOverYLinear -> zXOverYLinearMan,
    isConstant     -> zConstantMan
    ))
  val zEx0 = MuxCase(zYOverXLinearEx, Seq(
    isYOverXTable  -> zYOverXTableEx,
    isXOverYTable  -> zXOverYTableEx,
    isXOverYLinear -> zXOverYLinearEx,
    isConstant     -> zConstantEx
    ))
//   printf("zEx0             = %d\n", zEx0   )
//   printf("zMan0            = %b\n", zMan0  )

  // ---------------------------------------------------------------------
  // calc pi - atan(y/x) if x < 0

  // FIXME(perf) do not add exBias before this
  val zShift0  = -(zEx0.zext() - exBias.S)
  val zShift   = zShift0(log2ManW-1, 0)
  val zAligned = (1.U(1.W) ## zMan0 ## 0.U(2.W)) >> zShift
  val piManW1  = (((1L<<manW) + pi.man.toLong) << 3).U((manW+1+3).W)

  val piMinusZMan0 = piManW1 - zAligned

  val piMinusZMan0LessThan2 = piMinusZMan0(manW+3) === 0.U(1.W)

  val piMinusZMan =
    Mux(piMinusZMan0LessThan2, piMinusZMan0(manW+1+3-3, 2), piMinusZMan0(manW+1+3-2, 3)) +
    Mux(piMinusZMan0LessThan2, piMinusZMan0(1), piMinusZMan0(2))
  val piMinusZEx  = (1 + exBias).U(exW.W) - Mux(piMinusZMan0LessThan2, 1.U(1.W), 0.U(1.W))

//   printf("piMinusZMan      = %b(w = %d)\n", piMinusZMan , piMinusZMan.getWidth.U)
//   printf("piMinusZEx       = %b(w = %d)\n", piMinusZEx  , piMinusZEx .getWidth.U)

  // ---------------------------------------------------------------------
  // construct z

  val zMan = Mux(xsgn === 0.U(1.W), zMan0, piMinusZMan)
  val zEx  = Mux(xsgn === 0.U(1.W), zEx0,  piMinusZEx)
//   printf("xsgn             = %b(w = %d)\n", xsgn , xsgn.getWidth.U )
//   printf("zEx              = %d(w = %d)\n", zEx  , zEx .getWidth.U )
//   printf("zMan             = %b(w = %d)\n", zMan , zMan.getWidth.U )

  val z0 = MuxCase(Cat(zSgn, zEx, zMan), Seq(
    zzero     -> 0.U(spec.W.W),
    zpi       -> Cat(zSgn, pi       .ex.toLong.U(exW.W), pi       .man.toLong.U(manW.W)),
    zhalfpi   -> Cat(zSgn, pi_over_2.ex.toLong.U(exW.W), pi_over_2.man.toLong.U(manW.W)),
    z1piover4 -> Cat(zSgn, pi_1over4.ex.toLong.U(exW.W), pi_1over4.man.toLong.U(manW.W)),
    z3piover4 -> Cat(zSgn, pi_3over4.ex.toLong.U(exW.W), pi_3over4.man.toLong.U(manW.W)),
    znan      -> Cat(0.U(1.W), maskL(exW).U(exW.W), 1.U(1.W), 0.U((manW-1).W))
    ))

  io.z := ShiftRegister(z0, nStage)
  //printf("x=%x z=%x\n", io.x, io.z)
}

class ATan2F32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends ATan2Generic(RealSpec.Float32Spec,2,8,2,2,8,2,nStage,enableRangeCheck, enablePolynomialRounding) {
}
