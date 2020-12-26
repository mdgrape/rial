//% @file exponential.scala
//
// Exponential functions
// Copyright (C) Makoto Taiji RIKEN BDR 2020
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

//
// Power of 2, floating input, floating output
// x = 2^e * 1.m = xi + xf (xi = Int(x), xf = Frac(x))
// 2^x = 2^xi * 2^xf

class Pow2Generic(
  spec : RealSpec, // Input / Output floating spec
  nOrder: Int, adrW : Int, extraBits : Int, // Polynomial spec
  stage : PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends Module {

  val nStage = stage.total

  def getParam() = { (spec.exW, spec.manW, nOrder, adrW, extraBits, nStage,
    enableRangeCheck, enablePolynomialRounding) }

  def getStage() = nStage

  val expW = spec.exW
  val manW = spec.manW
  val W    = expW+manW+1

  val io = IO(iodef = new Bundle {
    val x   = Input(UInt(W.W))
    val z   = Output(UInt(W.W))
  })

  val calcW = manW+extraBits // for rounding
  val exValidW = log2Up(calcW+expW-1) // valid exponent bits

  val sgn = io.x(W-1)
  val ex  = io.x(W-2,manW)
  val man = io.x(manW-1,0)
  val exBias = spec.exBias

  val exOvfUdf = ex >= (expW-1+exBias).U
  val exOvf = exOvfUdf && !sgn
  val exUdf = exOvfUdf && sgn
  val nan   = if (spec.disableNaN) false.B else (ex === Fill(expW,1.U(1.W))) && ( man =/= 0.U )

  // Float to Integer, keep fractional width calcW
  // max exponent is expW-2
  val lsbPadding = expW-2+extraBits
  val manWith1 = 1.U(1.W) ## man ## 0.U(lsbPadding.W)
  val man1W = calcW+expW-1
  // manWith1 width: manW+extraBits+expW-1
  //   binary point: manW+extraBits

  // right shift amount :
  //   0  if ex==exBias+expW-2
  //   shift = exBias+expW-2-ex 
  //   0 if shift >= calcW+expW-1
  //     => ex <= exBias-calcW-1
  val shift_base = (exBias+expW-2) & maskI(exValidW)
  val shiftToZero = ex < (exBias-calcW).U(expW.W)
  val shift = shift_base.U - ex(exValidW-1,0)
  val xshift =manWith1 >> shift
  // binary point at bit calcW; extraBits bit for rounding / extension
  val xu = 0.U(1.W) ## Mux(shiftToZero, 0.U(man1W.W), xshift)
  val xi = Mux(sgn, -xu, xu) // calcW+expW
  val xiInt = xi(calcW+expW-1,calcW) // Integral part - for exponent
  //printf("ex=%d shiftToZero=%b shift=%d xi=%x\n", ex, shiftToZero, shift, xi)

  // Result exponent
  // zero when xi = 0x80??????/0x81??????/
  val zExNegMax = (xiInt(expW-1)===1.U) && (xiInt(expW-2,1)===0.U)
  val zEx0 = exBias.U + xiInt
  val zEx = Mux( exOvf || nan , 0xFF.U,
            Mux( exUdf || zExNegMax, 0.U, zEx0))

  val z0 = Wire(UInt((manW+expW).W))

  if (nOrder<=0) { // Fixed Table
    val adr = xi(calcW-1,0)
    var ovfFlag : Boolean = false
    val tbl = VecInit( (0L to (1L<<calcW)-1L).map(
      n => n.toDouble/(1L<<calcW)).map(
        x => {
          val y=round((pow(2.0,x)-1.0)*(1L<<manW))
          if (y>=(1L<<manW)) {
            ovfFlag = true
            println("WARNING: mantissa reaches to 2")
            maskL(manW).U(manW.W)
          } else {
            y.U(manW.W)
          }
        }
      )
    )
    val zman = tbl(adr)
    val zeroFlush = exOvfUdf || zExNegMax || nan
    z0 := zEx ## Mux(zeroFlush, nan ## 0.U((manW-1).W), zman)
  } else {
    // Fractional part by Polynomial
    val adr = xi(calcW-1,calcW-1-adrW+1)
    val d   = Cat(~xi(calcW-adrW-1),xi(calcW-1-adrW-1,0)).asSInt
    val dW  = calcW - adrW

    val tableD = new FuncTableDouble( x => pow(2.0,x)-1.0, nOrder )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    val tableI = new FuncTableInt( tableD, calcW )
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(1)
    val coeffAll = coeffTable(adr)
    val coeff    = getSlices(coeffAll, coeffWidth)
    // Coefficients should be all positive for 2^x
    val coeffS   = coeff.map( x => x.zext )

    def hornerC( c: SInt, z: SInt, dx: SInt, enableRounding: Boolean = false ) : SInt = {
      val zw   = z.getWidth
      val dxBp = dx.getWidth-1
      val dxw  = min(zw, dxBp)
      val dxl  = dx(dxBp, dxBp-dxw).asSInt
      val prod = dxl * z
      val prod_sft = dropLSB(dxw, prod)
      val sum = c +& prod_sft // Extend for safe
      //printf("sum=%x c=%x z=%x dx=%x dxl=%x prod=%x prod_sft=%x\n",sum, c, z, dx, dxl, prod, prod_sft)
      val dxlw=dxl.getWidth
      sum
    }

    def polynomialEvaluationC( c : Seq[SInt], dx : SInt, enableRounding: Boolean = false ) = {
      val res = c.init.foldRight(c.last)( (c,z) => hornerC( c, z, dx, enableRounding ) )
      res
    }

    val s0 = polynomialEvaluationC( coeffS, d, enablePolynomialRounding ).asUInt
    val zman  = dropLSB(extraBits, s0) + s0(extraBits-1)
    val polynomialOvf = zman.head(2) === 1.U // >=1
    val polynomialUdf = zman.head(1) === 1.U // Negative
    val zeroFlush = exOvfUdf || zExNegMax || polynomialUdf || nan
    val zman_checkovf = Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman(manW-1,0))
    z0 := zEx ## Mux(zeroFlush, nan ## 0.U((manW-1).W), zman_checkovf)
  }


  io.z   := 0.U(1.W) ## ShiftRegister(z0, nStage)
  //printf("x=%x z=%x\n", io.x, io.z)
}

class Pow2F32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends Pow2Generic(RealSpec.Float32Spec,2,8,2,nStage,enableRangeCheck, enablePolynomialRounding) {
}

object Pow2F32_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new Pow2F32(2)) ) )
}
