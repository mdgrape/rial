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

class Pow2Generic(
  expW : Int, manW : Int, // Input / Output floating spec
  nOrder: Int, adrW : Int, extraBits : Int, // Polynomial spec
  stage : PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends Module {

  val nStage = stage.total

  def getParam() = { (expW, manW, nOrder, adrW, extraBits, nStage,
    enableRangeCheck, enablePolynomialRounding) }

  def getStage() = nStage

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
  val exBias = (1L<<(expW-1))-1

  val exOvfUdf = ex >= (expW-1+exBias).U
  val exOvf = exOvfUdf && !sgn
  val exUdf = exOvfUdf && sgn
  val nan   = (ex === Fill(expW,1.U(1.W))) && ( man =/= 0.U )

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
  val shift_base = (exBias+expW-2) & mask(exValidW)
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
    val prod_sft = dropSLSB(dxw, prod)
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
  val zman  = dropULSB(extraBits, s0) + s0(extraBits-1) 
  val polynomialOvf = zman.head(2) === 1.U // >=1
  val polynomialUdf = zman.head(1) === 1.U // Negative 
  val zeroFlush = exOvfUdf || zExNegMax || polynomialUdf || nan
  val zman_checkovf = Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman(manW-1,0))

  val z0 = zEx ## Mux(zeroFlush, nan ## 0.U((manW-1).W), zman_checkovf) 

  io.z   := 0.U(1.W) ## ShiftRegister(z0, nStage)
  //printf("x=%x z=%x\n", io.x, io.z)
}

class Pow2F32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends Pow2Generic(8,23,2,8,2,nStage,enableRangeCheck, enablePolynomialRounding) {
}

//
// First we write fixed precision version, then generalize it
//
// Pow2_f32 : Power of 2 function, 32-bit single input / output
//            x = (-1)^s * 2^e * 1.m
//            -> xi = fixed point representation of x
//               overflow     if xi >=  0x80
//               underflow(0) if xi <= -0x7f
//               e >= 7, s = 0 -> overflow
//               e >= 7, s = 1 -> 0
//            z = 2^xi = 2^(xi.int) * 2^(xi.frac)
//               1 <= 2^(xi.frac) < 2

// Stage : Insert pipeline stages
//         Pipeline stages are inserted after the logics
//         as shift registers. If multiple registers are
//         inserted, they should be balanced in synthesis.

class Pow2F32_obsolete( nStage: Int = 0 ) extends Module {
  
  def getParam() = { (nStage) }

  val io = IO(iodef = new Bundle {
    val x   = Input(UInt(32.W))
    val z   = Output(UInt(32.W))
    val ovf = Output(UInt(1.W))
  })

  val expW = 8
  val manW = 23
  val W    = expW+manW+1
  val calcW = manW+1 // for rounding
  val sgn = io.x(W-1)
  val ex  = io.x(W-2,manW)
  val man = io.x(manW-1,0)
  val exBias = (1<<(expW-1))-1

  val exOvfUdf = ex >= (expW-1+exBias).U
  val exOvf = exOvfUdf && !sgn
  val exUdf = exOvfUdf && sgn
  val nan   = (ex === Fill(expW, 1.U(1.W))) && ( man =/= 0.U )

  // Float to Integer, range -128 ~ 127
  val manWith1 = Cat(1.U(1.W), man) 

  // right shift amount :
  //   0  if ex==6
  //   32 if ex==-26
  val shiftToZero = ex < (-26+exBias).U(8.W)
  val shift = ~(ex(4,0) - (0x1f+7).U)
  // 6-(ex-bias) = -(ex-(bias+7)) - 1 = ~(ex-(bias+7))+1-1
  val xshift =(1.U(1.W) ## man ## 0.U(8.W)) >> shift(4,0)
  // binary point at bit 24; 1 bit for rounding / extension
  val xu = Mux(shiftToZero, 0.U(31.W), xshift).pad(32)
  val xi = Mux(sgn, -xu, xu) // 32bit

  // Result exponent
  // zero when xi = 0x80??????/0x81??????/
  val zExNegMax = (xi(31,24)===0x80.U)||(xi(31,24)===0x81.U)
  val zEx0 = 0x7f.U + xi(31, 24) 
  val zEx = Mux( exOvf || nan , 0xFF.U,
            Mux( exUdf || zExNegMax, 0.U, zEx0))

  // Fractional part by 2nd Order Polynomial
  val nOrder = 2
  val adrW   = 8

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

  val prod1 = d * coeffS(2)
  val prod1_sft = dropSLSB(dW-1, prod1)
  val s1    = coeffS(1) +& prod1_sft // Extend for safe
  //println(prod1, prod1_sft, s1)

  val prod0 = d * s1
  val prod0_sft = dropSLSB(dW-1, prod0)
  val s0    = (coeffS(0) + prod0_sft).asUInt
  val zman  = dropULSB(1, s0) + s0(0) 
  val polynomialOvf = zman.head(1).asBool
  val zeroFlush = exOvfUdf || zExNegMax
  val zman_checkovf = Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman(manW-1,0))

  val z0 = zEx ## Mux(zeroFlush, 0.U(manW.W), zman_checkovf) 

  io.z   := 0.U(1.W) ## ShiftRegister(z0, nStage)
  io.ovf := ShiftRegister(exOvf, nStage)

}

object Pow2F32_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new Pow2F32(2)) ) )
}
