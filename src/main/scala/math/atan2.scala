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

// atan2(y, x), (floating, floating) input, floating output.
//
// = ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|
//   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y|
//   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y|
//   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y|
//           ^^^^^^^^^^^^^^^^^^^
//           this part is always positive in all the cases
//
class ATan2Generic(
  val spec : RealSpec, // Input / Output floating spec
  val nOrderRec  : Int, val adrWRec  : Int, val extraBitsRec  : Int, // polynomial spec
  val nOrderATan : Int, val adrWATan : Int, val extraBitsATan : Int,
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {
//   printf("=================================\n")

  val nStage = stage.total

  def getParam = { (spec.exW, spec.manW,
    nOrderRec, adrWRec, extraBitsRec, nOrderATan, adrWATan, extraBitsATan,
    nStage, enableRangeCheck, enablePolynomialRounding) }

  def getStage = nStage

  val manW  = spec.manW
  val exW   = spec.exW
  val exBias = spec.exBias

  val io = IO(new Bundle{
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

  val znan      =  (xnan ||  ynan) || ( xzero &&  yzero)
  val zzero     = ((xinf && !yinf) || (!xzero &&  yzero)) && (xsgn === 0.U(1.W))
  val zpi       = ((xinf && !yinf) || (!xzero &&  yzero)) && (xsgn === 1.U(1.W))
  val zhalfpi   = (!xinf &&  yinf) || ( xzero && !yzero)
  val z1piover4 = ( xinf &&  yinf  && ( xsgn === 0.U(1.W)))
  val z3piover4 = ( xinf &&  yinf  && ( xsgn === 1.U(1.W)))
  val isSpecialValue = znan || zzero || zpi || zhalfpi || z1piover4 || z3piover4
  assert(znan ^ zzero ^ zpi ^ zhalfpi ^ z1piover4 ^ z3piover4 ^ !isSpecialValue)

  val zSgn = ysgn

  // ==========================================================================
  // calc y / x or x / y
//   printf("calculating y/x and x/y ...\n")

  val swapped = io.x(spec.W-2, 0) < io.y(spec.W-2, 0) // without sign bit
  val xexS  = Mux(swapped, yExNobias, xExNobias)
  val yexS  = Mux(swapped, xExNobias, yExNobias)
  val xmanS = Mux(swapped, yman,      xman)
  val ymanS = Mux(swapped, xman,      yman)

  val yOverXEx0 = yexS -& xexS

  // calc 1.m/1.m \in (0.5, 2.0)

  val calcWRec = manW + extraBitsRec
  val one_over_x_manW1 = Wire(UInt((calcWRec+1).W))

  if(orderRec == 0) {
    val adrX = xmanS

    // 1.0 / 1.m -> (0.5, 1]
    val tbl = VecInit( (0L to (1L<<adrWRec)-1L).map(
      n => {
        val x = 1.0 + n.toDouble / (1L<<adrWRec)
        val y = math.round((1.0 / x) * (1L<<calcWRec))
        assert(y <= (1L << calcWRec))
        y.U((calcWRec+1).W)
      }
    ) )

    one_over_x_manW1 := tbl(adrX)

  } else {

    val adrX = ~xmanS(manW-1, manW-adrWRec)
    val dX  = Cat(xmanS(manW-adrWRec-1),~xmanS(manW-adrWRec-2,0)).asSInt

    // Create 1/x table for inverted x
    val eps = pow(2.0, -manW)
    val tableD = new FuncTableDouble( x => 1.0/(2.0-(x+eps)), orderRec )
    tableD.addRange(0.0, 1.0, 1<<adrWRec)
    val tableI = new FuncTableInt( tableD, bpRec )
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(1)

    val coeffXAll = coeffTable(adrX)
    val coeffX    = getSlices(coeffXAll, coeffWidth)
    val coeffXS   = coeffX.map( x => x.zext ) // all positive (x is inverted)

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

    val oneOverXOvf = sX0 > (1<<bpRec).U
    val oneOverXUdf = sX0 < (1<<(bpRec-1)).U

    one_over_x_manW1 := Mux(xmanS === 0.U, (1 << bpRec).U((bpRec + 1).W),
      MuxCase(sX0, Seq(
        oneOverXOvf -> (1 <<  bpRec   ).U((bpRec + 1).W),
        oneOverXUdf -> (1 << (bpRec-1)).U((bpRec + 1).W)
      )))
  }
  val ymanW1     = 1.U(1.W) ## ymanS
  val yOverXMan0 = ymanW1 * one_over_x_manW1

  assert(yOverXMan0.getWidth == bpRec+1+manW+1)
  assert(yOverXMan0(bpRec+1+manW+1-1) === 0.U(1.W)) // < 2

//   printf("1/x_man = %b\n", one_over_x_manW1)

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
  assert(yOverXMan0(bpRec+1 + manW+1 - 1) === 0.U(1.W))
  val yOverXLessThan1 = yOverXMan0(bpRec+1 + manW+1 - 2) === 0.U(1.W)
  val roundBits   = bpDiv - manW
  val yOverXShift = Mux(yOverXLessThan1, yOverXMan0(bpDiv-1, roundBits-1), yOverXMan0(bpDiv, roundBits))
  val yOverXRound = Mux(yOverXLessThan1, yOverXMan0(roundBits-2), yOverXMan0(roundBits-1))

  val yOverXManW1Round = yOverXShift +& yOverXRound
  val yOverXMoreThan2 = yOverXManW1Round(manW+1) === 1.U

  val sameMan = xman === yman

  val yOverXManW1 = Mux(sameMan,         (1<<manW).U((manW+1).W),
                    Mux(yOverXMoreThan2, yOverXManW1Round(manW+1, 1) + yOverXManW1Round(0),
                                         yOverXManW1Round(manW, 0)))
  val yOverXMan   = yOverXManW1(manW-1, 0);

  assert(yOverXManW1.getWidth == manW+1)
  assert(yOverXManW1(manW) === 1.U(1.W))

  val yOverXEx = Mux(
    !sameMan && yOverXLessThan1 && !yOverXMoreThan2, yOverXEx0 - 1.S, yOverXEx0)

//   printf("y/x ManW1Round = %b\n", yOverXManW1Round)
//   printf("y/x ManW1 = %b\n", yOverXManW1)
//   printf("y/x Ex    = %d\n", yOverXEx)

//   printf("y/x = %d|%d|%b\n", ysgn ^ xsgn, yOverXEx, yOverXMan)

  // ==========================================================================
  // calculate atan

  // --------------------------------------------------------------------------
  // linear (x^3 / 6 < x * 2^-manW)

//   printf("calculating linear case ...\n")

  val linearThreshold = -math.round(math.ceil(manW / 2.0 + 1.0)) // -13, if FP32
  val isYOverXLinear  = yOverXEx < linearThreshold.S

  val atanYOverXLinearEx    = yOverXEx    // Wire(SInt((exW+1).W))
  val atanYOverXLinearManW1 = yOverXManW1 // Wire(UInt((manW+1).W))
//   atanYOverXLinearEx    := 
//   atanYOverXLinearManW1 := 

  // --------------------------------------------------------------------------
  // use table (2^linearThreshold <= x < 2^-linearThreshold)
//   printf("calculating table case ...\n")

  val isYOverXTable = (linearThreshold.S <= yOverXEx) && (yOverXEx < 0.S) && !isSpecialValue

  val exAdrW   = log2Up(abs(linearThreshold).toInt)
  val exAdrMax = -linearThreshold - 1
  assert(exAdrMax > 0)
  val exAdrYOverX0 = (-yOverXEx - 1.S)(exAdrW-1, 0)
  val exAdrYOverX  = Mux(exAdrYOverX0 > exAdrMax.U, 0.U(exAdrW.W), exAdrYOverX0)

  val calcW = manW + extraBitsATan
  val atanYOverXTableEx    = Wire(SInt((exW+1).W))
  val atanYOverXTableManW1 = Wire(UInt((manW+1).W))

  if(orderATan == 0) {
    val adrYOverX = yOverXMan

    val tbl = VecInit( (-1 to linearThreshold.toInt by -1 ).map( exponent => {
      VecInit((0L to (1L<<adrWATan)-1L).map( n => {
        // atan(x) < x.
        val x = scalb(1.0 + n.toDouble / (1L<<adrWATan), exponent)
        val y = math.round(scalb(atan(x), -exponent-1) * (1L<<calcW))
        assert(y < (1L << calcW))
        y.U((calcW+1).W)
      } ) )
    } ) )

    val atanYOverXScaled = tbl(exAdrYOverX)(adrYOverX)
    assert(calcW+1 == atanYOverXScaled.getWidth)
    val atanYOverXShift = PriorityEncoder(Reverse(atanYOverXScaled))

    val resYOverXEx  = yOverXEx + 1.S - atanYOverXShift.zext
    val resYOverXMan = (atanYOverXScaled << atanYOverXShift)(calcW, calcW-manW)
    assert(resYOverXMan.getWidth == manW+1)
    assert(resYOverXMan(manW) === 1.U)

    atanYOverXTableEx    := resYOverXEx
    atanYOverXTableManW1 := resYOverXMan

  } else {

    val dxbp = manW - adrWATan - 1
    val adrYOverX = yOverXMan(manW-1, manW-adrWATan)
    val dYOverX   = Cat(~yOverXMan(manW-adrWATan-1), yOverXMan(manW-adrWATan-2,0)).asSInt

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
    val coeffYOverX      = getSlices(coeffYOverXTable(adrYOverX), coeffWidth)
    val coeffSYOverX     = coeffYOverX.map( x => x.asSInt )

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

    val log2CalcW    = log2Up(calcW)
    val s0YOverXLenS = s0YOverX.getWidth.S - PriorityEncoder(Reverse(s0YOverX)).zext()
    val s0YOverXLen  = s0YOverXLenS(log2CalcW-1, 0)
    assert(0.S <= s0YOverXLenS)

    val s0YOverXShift = (calcW+1).U(log2CalcW.W) - s0YOverXLen
    assert(0.U <= s0YOverXShift && s0YOverXShift <= calcW.U)

    val resYOverXEx  = yOverXEx + 1.S - s0YOverXShift.zext
    val resYOverXMan = (s0YOverX(calcW, 0) << s0YOverXShift)(calcW, 0)

    // ------------------------------------------------------------------------
    // shift/round atan(y/x) and return it

    atanYOverXTableEx := resYOverXEx
    if(extraBitsATan > 0) {
      atanYOverXTableManW1 := resYOverXMan(manW+extraBitsATan, extraBitsATan) + resYOverXMan(extraBitsATan-1)
    } else {
      atanYOverXTableManW1 := resYOverXMan
    }
  }

  // --------------------------------------------------------------------------
  // merge atan results

  val pi_1over4 = new RealGeneric(spec, Pi * 0.25) // |x| === |y|
  val atanEx    = Wire(SInt((exW+1).W))
  val atanManW1 = Wire(UInt((manW+1).W))
  atanEx := MuxCase((pi_1over4.ex.toLong).S((exW+1).W), Seq(
    isYOverXLinear -> atanYOverXLinearEx,
    isYOverXTable  -> atanYOverXTableEx
    ))
  atanManW1 := MuxCase((pi_1over4.man.toLong + (1L<<manW)).U((spec.manW+1).W), Seq(
    isYOverXLinear -> atanYOverXLinearManW1,
    isYOverXTable  -> atanYOverXTableManW1
    ))
  val atanMan = atanManW1(manW-1, 0)
  val atanExBiased0 = exBias.S((exW+1).W) + atanEx
  val atanExBiased  = Mux(atanExBiased0 < 0.S, 0.U(exW), atanExBiased0(exW-1,0))

//   printf("atanYOverXLinearEx = %d\n", atanYOverXLinearEx)
//   printf("atanYOverXTableEx  = %d\n", atanYOverXTableEx )
//   printf("pi/4.ex            = %d\n", pi_1over4.ex.toLong.S)

  // ==========================================================================
  // select correct one
  // = ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|
  //   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y|
  //   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y|
  //   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y|
//   printf("selecting ...\n")

  val pi_over_2 = new RealGeneric(spec, Pi * 0.50)
  val pi_3over4 = new RealGeneric(spec, Pi * 0.75)
  val pi        = new RealGeneric(spec, Pi)

  assert(isYOverXLinear ^ isYOverXTable ^ znan)
//   printf("isYOverXLinear = %b, zYOverXLinearEx = %b, zYOverXLinearMan = %b\n", isYOverXLinear, zYOverXLinearEx, zYOverXLinearMan)
//   printf("isYOverXTable  = %b, zYOverXTableEx  = %b, zYOverXTableMan  = %b\n", isYOverXTable , zYOverXTableEx , zYOverXTableMan )

  // --------------------------------------------------------------------------
  // calc pi - atan

  val log2ManW3     = log2Up(1+manW+3)
  val atanShift0Pi  = -atanEx
  val atanShiftPi   = atanShift0Pi(log2ManW3-1, 0)
  val atanAlignedPi = Mux(atanShift0Pi > (manW+1+3).S, 0.U((manW+1+3).W),
                          (0.U(1.W) ## atanManW1 ## 0.U(2.W)) >> atanShiftPi)
  val piManW1       = (((1L<<manW) + pi.man.toLong) << 3).U((manW+1+3).W)

  val piMinusATanMan0 = piManW1 - atanAlignedPi

//   printf("piManW1      = %b(w = %d)\n", piManW1      , piManW1     .getWidth.U)
//   printf("zAligned     = %b(w = %d)\n", zAligned     , zAligned    .getWidth.U)
//   printf("piMinusZMan0 = %b(w = %d)\n", piMinusZMan0 , piMinusZMan0.getWidth.U)

  val piMinusATanMan0LessThan2 = piMinusATanMan0(manW+3) === 0.U(1.W)

  val piMinusATanMan = Mux(piMinusATanMan0LessThan2,
    piMinusATanMan0(manW+1+3-3, 2) + piMinusATanMan0(1),
    piMinusATanMan0(manW+1+3-2, 3) + piMinusATanMan0(2))
  // roundmorethan2?
  val piMinusATanEx  = exBias.U(exW.W) + (!piMinusATanMan0LessThan2).asUInt

  // --------------------------------------------------------------------------
  // calc pi/2 +/- atan
  val atanShift0HalfPi  = -atanEx
  val atanShiftHalfPi   = atanShift0HalfPi(log2ManW3-1, 0)
  val atanAlignedHalfPi = Mux(atanShift0HalfPi > (manW+1+3).S, 0.U((manW+1+3).W),
                              (atanManW1 ## 0.U(3.W)) >> atanShiftHalfPi)
  val halfPiManW1       = (((1L<<manW) + pi_over_2.man.toLong) << 3).U((manW+1+3).W)

//   printf("halfPiManW1 = %b\n", halfPiManW1)
//   printf("atanAligned = %b\n", atanAlignedHalfPi)
//   printf("atanEx      = %b\n", atanEx)
//   printf("atanManW1   = %b\n", atanManW1)

  val halfPiMinusATanMan0 = halfPiManW1 - atanAlignedHalfPi
  val halfPiMinusATanMan0LessThan1 = halfPiMinusATanMan0(manW+3) === 0.U
  val halfPiMinusATanManRound = Mux(halfPiMinusATanMan0LessThan1,
    halfPiMinusATanMan0(manW+1+3-2, 2) +& halfPiMinusATanMan0(1),
    halfPiMinusATanMan0(manW+1+3-1, 3) +& halfPiMinusATanMan0(2))
  val halfPiMinusATanManRoundMoreThan2 = halfPiMinusATanManRound(manW+1)
  val halfPiMinusATanEx  = exBias.U(exW.W) - halfPiMinusATanMan0LessThan1.asUInt + halfPiMinusATanManRoundMoreThan2.asUInt
  val halfPiMinusATanMan = halfPiMinusATanManRound(manW-1, 0)
//   printf("halfPiMinusATanMan = %b\n", halfPiMinusATanMan)

  val halfPiPlusATanMan0  = halfPiManW1 +& atanAlignedHalfPi
  val halfPiPlusATanMan0MoreThan2 = halfPiPlusATanMan0(manW+4) === 1.U
  val halfPiPlusATanManRound = Mux(halfPiPlusATanMan0MoreThan2,
    halfPiPlusATanMan0(manW+1+3,   4) +& halfPiPlusATanMan0(3),
    halfPiPlusATanMan0(manW+1+3-1, 3) +& halfPiPlusATanMan0(2))
  val halfPiPlusATanManRoundMoreThan2 = halfPiPlusATanManRound(manW+1)
  val halfPiPlusATanEx   = exBias.U(exW.W) + halfPiPlusATanMan0MoreThan2.asUInt + halfPiPlusATanManRoundMoreThan2.asUInt
  val halfPiPlusATanMan  = halfPiPlusATanManRound(manW-1, 0)

  // --------------------------------------------------------------------------

//   printf("     atan = %d|%b\n", atanExBiased.zext  - exBias.S, atanMan        )
//   printf("pi  -atan = %d|%b\n", piMinusATanEx.zext - exBias.S, piMinusATanMan )
//   printf("pi/2-atan = %d|%b\n", halfPiMinusATanEx.zext - exBias.S, halfPiMinusATanMan)
//   printf("pi/2+atan = %d|%b\n", halfPiPlusATanEx.zext  - exBias.S, halfPiPlusATanMan)
  val zMan = Mux(swapped,
    Mux(xsgn === 1.U, halfPiPlusATanMan, halfPiMinusATanMan),
    Mux(xsgn === 1.U, piMinusATanMan,    atanMan)
    )
  val zEx  = Mux(swapped,
    Mux(xsgn === 1.U, halfPiPlusATanEx, halfPiMinusATanEx),
    Mux(xsgn === 1.U, piMinusATanEx,    atanExBiased)
    )

  // ---------------------------------------------------------------------
  // construct z

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
