//% @file acos.scala
//
// ACos function
// Copyright (C) Makoto Taiji RIKEN BDR 2020
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

// acos(x), floating input, floating output

class ACosGeneric(
  val spec : RealSpec, // Input / Output floating spec
  val nOrder: Int, val adrW : Int, val extraBits : Int, // Polynomial spec
  val stage : PipelineStageConfig,
  val enableRangeCheck : Boolean = true,
  val enablePolynomialRounding : Boolean = false,
) extends Module {

  val nStage = stage.total

  def getParam() = { (spec.exW, spec.manW, nOrder, adrW, extraBits, nStage,
    enableRangeCheck, enablePolynomialRounding) }

  def getStage() = nStage

  val manW  = spec.manW
  val exW   = spec.exW
  val exBias = spec.exBias

  val io = IO(iodef = new Bundle {
    val x   = Input (UInt(spec.W.W))
    val z   = Output(UInt(spec.W.W))
  })

  val order =
    if (adrW>=manW) {
      if (nOrder != 0)
        println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
      0
    } else {
      nOrder
    }

  val bp = manW+extraBits // for rounding

  val (xsgn,  xex,  xman) = FloatChiselUtil.decompose(spec, io.x)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(spec, io.x)
  val xExNobias = xex.zext() - exBias.S((exW+1).W)
  printf("x = %d|%d(%d)|%b\n", xsgn, xex, xExNobias, xman)

  val znan = xnan || xinf || (xExNobias === 0.S && xman =/= 0.U) || (0.S < xExNobias)
  val zzero = (xExNobias === 0.S && xman === 0.U) && xsgn === 0.U
  val zpi   = (xExNobias === 0.S && xman === 0.U) && xsgn === 1.U
  assert(znan ^ zzero ^ zpi ^ (!(znan || zzero || zpi)))

  val zSgn = 0.U(1.W)

  // hereafter, xex < 0
  //
  // acos(x) = pi/2 - x - x^3/6 - 3x^5/40
  //         < pi/2 - x - x^3/4 -  x^5/8
  val constThreshold  = -manW                                    // -23, if FP32
  val linearThreshold = math.round(math.ceil((-manW - 1) / 3.0)) //  -8

  val isConstant = xExNobias < constThreshold .S((exW+1).W)
  val isLinear   = xExNobias < linearThreshold.S((exW+1).W)

  // --------------------------------------------------------------------------
  // constant

  val halfpi = new RealGeneric(spec, Pi * 0.5)
  val zExConstant  = halfpi.ex .U(exW.W)
  val zManConstant = halfpi.man.toLong.U(manW.W)

  // --------------------------------------------------------------------------
  // linear (pi/2 - x), pi/2 = 1.57... = 2^0 * 1.57...

  val constThresholdDigit = log2Up(abs(constThreshold))
//   printf(f"constThreshold      = ${constThreshold     }\n")
//   printf(f"constThresholdDigit = ${constThresholdDigit}\n")

  val xmanW1          = xman + (1<<manW).U((4+manW).W)
  val linearExDiff    = (~xExNobias(constThresholdDigit-1, 0)) - 2.U(constThresholdDigit.W)
  val linearXAligned  = Mux(xsgn === 1.U, (xmanW1 >> linearExDiff), ~(xmanW1 >> linearExDiff) + 1.U)
  val linearManDiff   = (((1<<manW) + halfpi.man.toLong) << 3).U((manW+1+3).W) + linearXAligned
  // since pi/2 = 1.1001... and xEx < linearThreshold, mandiff always larger than 1 and never be larger than 2.

  val linearLSB    = linearManDiff(3)
  val linearRound  = linearManDiff(2)
  val linearSticky = linearManDiff(1) | linearManDiff(0)
  val linearInc    = FloatChiselUtil.roundIncBySpec(
    RoundSpec.roundToEven, linearLSB, linearRound, linearSticky)

  assert(linearManDiff(manW+3-1, 3) =/= maskL(manW).U(manW.W))

  val zExLinear  = exBias.U(exW.W)
  val zManLinear = linearManDiff(manW+3-1, 3) + linearInc

  // --------------------------------------------------------------------------
  // otherwise, use table

  val zman0 = Wire(UInt(manW.W))
  val zex0  = Wire(UInt(exW.W))

  val exAdr0 = -xExNobias-1.S
  val exAdr  = Mux(exAdr0 < 0.S, 0.U, exAdr0.asUInt) // if x == 1.0, exAdr0 == -1.

  val calcW = manW + extraBits
  if (order<=0) { // Fixed Table
    val adr = xman(adrW-1,0)
    val tbl = VecInit( (-1 to linearThreshold.toInt by -1).map( exponent => {
      VecInit( (0L to (1L<<adrW)-1L).map(
        n => {
          val x = scalb(1.0 + n.toDouble/(1L<<adrW), exponent.toInt)
          val y = round(acos(x)*(1L<<calcW))
          //println(f"$n $x $y")
          y.U((calcW+1).W)
        } ) )
      } ) )

    val piFixed = (math.round(Pi * (1<<calcW))).toLong.U((calcW+2).W)
    val res00 = tbl(exAdr)(adr)
    val res0  = Mux(xsgn === 1.U, piFixed - res00, 0.U(1.W) ## res00)

    val shift0 = PriorityEncoder(Reverse(res0))
    assert(shift0 < (calcW+1).U)
    val log2CalcW = log2Up(calcW+1)
    val shift = shift0(log2CalcW-1, 0)

    val res   = res0 << shift
    assert(res(calcW+1) === 1.U(1.W))

    printf("adr        = %b\n", adr       )
    printf("res0       = %b\n", res0      )
    printf("shift      = %b\n", shift     )
    printf("resShifted = %b\n", res(calcW, calcW-manW+1))

    val resRound = res(calcW, calcW-manW+1) +& res(calcW-manW)

    zex0  := (exBias+1).U(exW.W) - shift
    zman0 := Mux(resRound(manW) === 1.U, maskL(manW).U(manW.W), resRound(manW-1, 0))
  } else {
    val adr = xman(manW-1, manW-adrW)
    val d   = Cat(~xman(manW-adrW-1), xman(manW-adrW-2,0)).asSInt

    printf("man = %b\n", xman)
    printf("adr = %b\n", adr)
    printf(" d  = %b\n", d)

    val coeffWidth = (-1 to linearThreshold.toInt by -1).map(exponent => {
      val tableD = new FuncTableDouble( x => (Pi * 0.5) - acos(scalb(1.0 + x, exponent)), order )
      tableD.addRange(0.0, 1.0, 1<<adrW)
      val tableI = new FuncTableInt( tableD, bp ) // convert float table into int
      val w = tableI.getCBitWidth(/*sign mode = */1)
      println(f"//  width : "+w.mkString("", ", ", ""))
      w
    }).reduce( (lhs, rhs) => {
      lhs.zip(rhs).map( x => max(x._1, x._2))
    })
    val coeffTables = VecInit((-1 to linearThreshold.toInt by -1).map(exponent => {
      val tableD = new FuncTableDouble( x => (Pi * 0.5) - acos(scalb(1.0 + x, exponent)), order )
      tableD.addRange(0.0, 1.0, 1<<adrW)
      val tableI = new FuncTableInt( tableD, bp ) // convert float table into int
      tableI.getVectorWithWidth(coeffWidth, /*sign mode = */1)
    }))
    val coeffTable = coeffTables(exAdr)
    val coeff  = getSlices(coeffTable(adr), coeffWidth)
    val coeffS = coeff.map( x => x.zext )

    def hornerC( c: SInt, z: SInt, dx: SInt, enableRounding: Boolean = false ) : SInt = {
      val zw   = z.getWidth
      val dxBp = dx.getWidth-1
      val dxw  = min(zw, dxBp)
      val dxl  = dx(dxBp, dxBp-dxw).asSInt
      val prod = dxl * z
      val prod_sft = dropLSB(dxw, prod)
      val sum = c +& prod_sft // Extend for safe
      //printf("sum=%x c=%x z=%x dx=%x dxl=%x prod=%x prod_sft=%x\n",sum, c, z, dx, dxl, prod, prod_sft)
      sum
    }

    def polynomialEvaluationC( c : Seq[SInt], dx : SInt, enableRounding: Boolean = false ) = {
      val res = c.init.foldRight(c.last)( (c,z) => hornerC( c, z, dx, enableRounding ) )
      res
    }
    val halfPiFixed = math.round(Pi * 0.5 * (1 << calcW)).U((calcW+1).W)

    val res0S = polynomialEvaluationC( coeffS, d, enablePolynomialRounding ).asUInt
    val res0  = Mux(res0S(res0S.getWidth-1), 0.U((res0S.getWidth-1).W), res0S(res0S.getWidth-2, 0))
    val res   = halfPiFixed + Mux(xsgn === 1.U(1.W), res0, ~res0 + 1.U)
    val shift = (calcW+2).U - (res.getWidth.U - PriorityEncoder(Reverse(res)))
    val resShifted = (res << shift)(calcW+1, 1) - (1<<calcW).U

    printf("cir:res0       = %b(w=%d)\n", res0, res0.getWidth.asUInt)
    printf("cir:halfpi     = %b\n", halfPiFixed)
    printf("cir:res        = %b\n", res)
    printf("cir:shift      = %d\n", shift)
    printf("cir:resShifted = %b\n", resShifted)

    zex0  := (exBias+1).U(exW.W) - shift
    zman0 := (resShifted >> extraBits) + resShifted(extraBits-1)
  }

  val zEx  = Mux(isConstant, zExConstant,  Mux(isLinear, zExLinear,  zex0))
  val zMan = Mux(isConstant, zManConstant, Mux(isLinear, zManLinear, zman0))

  val pi = new RealGeneric(spec, Pi)

  val z0 = MuxCase(Cat(zSgn, zEx, zMan), Seq(
    zzero -> Cat(0.U(1.W), 0.U(exW.W), 0.U(manW.W)),
    zpi   -> Cat(0.U(1.W), pi.ex.toLong.U(exW.W), pi.man.toLong.U(manW.W)),
    znan  -> Cat(0.U(1.W), maskL(exW).U(exW.W), 1.U(1.W) ## 0.U((manW-1).W))
    ))
  io.z   := ShiftRegister(z0, nStage)
  //printf("x=%x z=%x\n", io.x, io.z)
}

class ACosF32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends ACosGeneric(RealSpec.Float32Spec,2,8,2,nStage,enableRangeCheck, enablePolynomialRounding) {
}
