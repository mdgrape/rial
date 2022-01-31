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

  val expW = spec.exW
  val manW = spec.manW
  val W    = expW+manW+1

  val io = IO(new Bundle{
    val x   = Input(UInt(W.W))
    val z   = Output(UInt(W.W))
  })

  val padding = if (adrW>=manW) {
    if (nOrder != 0) {
      println("ERROR: table address width >= mantissa width, but polynomial order is not zero.")
      sys.exit(1)
      0
    } else {
      adrW-manW
    }
  } else {
    extraBits
  }

  val bp    = manW+padding // for rounding
  val exValidW = log2Up(bp+expW-1) // valid exponent bits

  val sgn = io.x(W-1)
  val ex  = io.x(W-2,manW)
  val man = io.x(manW-1,0)
  val exBias = spec.exBias

  val exOvfUdf = ex >= (expW-1+exBias).U
  val exOvf = exOvfUdf && !sgn
  val exUdf = exOvfUdf && sgn
  val nan   = if (spec.disableNaN) false.B else (ex === Fill(expW,1.U(1.W))) && ( man =/= 0.U )

  // Float to Integer, keep fractional width bp
  // max exponent is expW-2
  val lsbPadding = expW-2+padding
  val manWith1 = 1.U(1.W) ## man ## 0.U(lsbPadding.W)
  val man1W = manW+padding+expW-1
  // manWith1 width: manW+extraBits+expW-1
  //   binary point: manW+extraBits

  // right shift amount :
  //   0  if ex==exBias+expW-2
  //   shift = exBias+expW-2-ex 
  //   0 if shift >= bp+expW-1
  //     => ex <= exBias-bp-1
  val shift_base = (exBias+expW-2) & maskI(exValidW)
  val shiftToZero = ex < (exBias-bp).U(expW.W)
  val shift = shift_base.U - ex(exValidW-1,0)
  val xshift =manWith1 >> shift
  // bp includes extraBits bit for rounding / extension
  val xu = 0.U(1.W) ## Mux(shiftToZero, 0.U(man1W.W), xshift)
  val xi = Mux(sgn, -xu, xu) // bp+expW
  val xiInt = xi(bp+expW-1,bp) // Integral part - for exponent
  //printf("ex=%d shiftToZero=%b shift=%d xi=%x\n", ex, shiftToZero, shift, xi)

  // Result exponent
  // zero when xi = 0x80??????/0x81??????/
  val zExNegMax = (xiInt(expW-1)===1.U) && (xiInt(expW-2,1)===0.U)
  val zEx0 = exBias.U + xiInt
  val zEx = Mux( exOvf || nan , Fill(expW, 1.U(1.W)),
            Mux( exUdf || zExNegMax, 0.U, zEx0))

  val z0 = Wire(UInt((manW+expW).W))

  if (nOrder<=0) { // Fixed Table
    val adr = xi(adrW-1,0)
    var ovfFlag : Boolean = false
    val tbl = VecInit( (0L to (1L<<adrW)-1L).map(
      n => {
        val x = n.toDouble/(1L<<adrW)
        val y=round((pow(2.0,x)-1.0)*(1L<<manW))
        //println(f"$n $x $y")
        if (y>=(1L<<manW)) {
          ovfFlag = true
          println("WARNING: mantissa reaches to 2")
          maskL(manW).U(manW.W)
        } else {
          y.U(manW.W)
        }
      }
    ) )
    val zman = tbl(adr)
    val zeroFlush = exOvfUdf || zExNegMax || nan
    z0 := zEx ## Mux(zeroFlush, nan ## 0.U((manW-1).W), zman)
  } else {
    // Fractional part by Polynomial
    val adr = xi(bp-1,bp-1-adrW+1)
    val d   = Cat(~xi(bp-adrW-1),xi(bp-1-adrW-1,0)).asSInt
    val dW  = bp - adrW

    val tableD = new FuncTableDouble( x => pow(2.0,x)-1.0, nOrder )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    val tableI = new FuncTableInt( tableD, bp )
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

class ExponentialGeneric(
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

  val expW = spec.exW
  val manW = spec.manW
  val W    = expW+manW+1

  val io = IO(new Bundle{
    val x   = Input(UInt(W.W))
    val z   = Output(UInt(W.W))
  })

  val xm = 1.U(1.W) ## io.x(manW-1,0)
  val xe = io.x(W-2,manW)
  val xs = io.x(W-1)
  val zeroOrInfOrNaN = xe.andR || !xe.orR
  val xePlus1 = Mux(zeroOrInfOrNaN, xe, xe+1.U)

  // Multiply log_2 e
  val log2e = 1.0d/java.lang.Math.log(2.0d)
  val log2eMan = BigInt(round(java.lang.Math.scalb(log2e, manW))).U((manW+1).W)

  val prod = xm * log2eMan
  val le2  = prod(manW*2+1)
  val ym = Mux(le2, prod(manW*2,manW+1), prod(manW*2-1, manW))
  val r  = Mux(le2, prod(manW), prod(manW-1))
  val stickey = prod(manW-2,0).orR || (le2 && prod(manW-1))
  val inc = r && (stickey || ym(0))
  val ymRound = ym +& inc
  val ye = Mux(ymRound(manW) | le2, xePlus1, xe)

  val yman = Mux(zeroOrInfOrNaN, io.x(manW-1,0), ymRound(manW-1,0))
  val y = xs ## ye ## yman

  val pow2 = Module(new Pow2Generic( spec, nOrder, adrW, extraBits, stage,
    enableRangeCheck, enablePolynomialRounding ) )

  //printf("%x %x %x %x %x %x\n", io.x, prod, ym, y, ye, pow2.io.z)

  pow2.io.x := y
  io.z := pow2.io.z

}

////////////////////////////////////////////////////////////////////////
// The followings are example implementations

class Pow2F32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends Pow2Generic(RealSpec.Float32Spec,2,8,2,nStage,enableRangeCheck, enablePolynomialRounding) {
}

class Pow2BF16(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends Pow2Generic(RealSpec.BFloat16Spec,0,8,0,nStage,enableRangeCheck, enablePolynomialRounding) {
}

class ExpF32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends ExponentialGeneric(RealSpec.Float32Spec,2,8,2,nStage,enableRangeCheck, enablePolynomialRounding) {
}

class ExpBF16(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends ExponentialGeneric(RealSpec.BFloat16Spec,0,8,0,nStage,enableRangeCheck, enablePolynomialRounding) {
}

object Pow2F32_driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new Pow2F32(2)) ) )
}
