//% @file reciprocal.scala
//
// Reciprocal function
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

// 1/x, floating input, floating output
// x = 2^e * 1.m = 
// 1/x = 2^(-e-1) * 2/1.m

class ReciprocalGeneric(
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

  val io = IO(iodef = new Bundle {
    val x   = Input(UInt(W.W))
    val z   = Output(UInt(W.W))
  })

  val order =
    if (adrW>=manW) {
      if (nOrder != 0)
        println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
      0
    } else {
      nOrder
    }

  val bp    = manW+extraBits // for rounding

  val sgn = io.x(W-1)
  val ex  = io.x(W-2,manW)
  val man = io.x(manW-1,0)
  val exBias = spec.exBias

  val manNotZero = man.orR.asUInt
  val exZero = !ex.orR
  val exInf  = ex.andR
  val nan   = if (spec.disableNaN) false.B else (ex === Fill(expW,1.U(1.W))) && ( man =/= 0.U )

  // 2/(1+e) - 2(1-0.5*e) = 2/(1+e) [ 1 - (1-0.5e)(1+e) ]
  //   = 2/(1+e) [ -0.5e + 0.5 e^2 ] <0
  // there is no possibility results reaches to 2
  //   if mantissa is not zero.

  // In usual case, there is possibility of underflow
  // but no overflow is expected.
  // Need to care either one of overflow or underflow
  val invExMax = exBias-(1-exBias)
  val invExMin = exBias-(maskI(expW)-1-exBias)-1
  val zex0 = (exBias<<1).U((expW+1).W)-ex-manNotZero
  val ovfudf = zex0(expW)

  val (inf,zero) = 
    if (invExMax>=maskI(expW)) {
      // Possible overflow, no underflow
      ( ovfudf || nan || exZero, exInf)
    } else if (invExMin<=0) {
      // Possible underflow, no overflow
      ( nan || exZero, ovfudf || exInf)
    } else {
      // Never happens
      ( nan || exZero, exInf)
    }

  val zEx = Mux(inf, maskU(expW),
    Mux(zero, 0.U(expW.W), zex0(expW-1,0)))

  val zSgn = sgn && !(inf||zero)

  val z0 = Wire(UInt(W.W))

  if (order<=0) { // Fixed Table
    val adr = man(adrW-1,0)
    val tbl = VecInit( (0L to (1L<<adrW)-1L).map(
      n => {
        val x = 1.0+n.toDouble/(1L<<adrW)
        val y = round((2.0/x-1.0)*(1L<<manW))
        //println(f"$n $x $y")
        if ((n!=0)&&(y>=(1L<<manW))) {
          println("WARNING: mantissa reaches to 2")
          maskL(manW).U(manW.W)
        } else {
          y.U(manW.W)
        }
      }
    ) )
    val zman = tbl(adr)
    val zeroFlush = inf || zero || nan
    z0 := zSgn ## zEx ## Mux(zeroFlush, nan ## 0.U((manW-1).W), zman)
  } else {
    // Fractional part by Polynomial
    val adr = ~man(manW-1, manW-adrW)
    val d   = Cat(man(manW-adrW-1),~man(manW-adrW-2,0)).asSInt

    // Create table for 1/(2-x)
    val eps = pow(2.0, -manW)
    val tableD = new FuncTableDouble( x => 2.0/(2.0-(x+eps))-1.0, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    val tableI = new FuncTableInt( tableD, bp )
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(1)
    val coeffAll = coeffTable(adr)
    val coeff    = getSlices(coeffAll, coeffWidth)
    // Coefficients should be all positive for 2(2-x)
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
    val zeroFlush = inf || zero || polynomialUdf || nan
    val zman_checkovf = Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman(manW-1,0))
    z0 := zSgn ## zEx ## Mux(zeroFlush, nan ## 0.U((manW-1).W), zman_checkovf)
  }

  io.z   := ShiftRegister(z0, nStage)
  //printf("x=%x z=%x\n", io.x, io.z)
}

class ReciprocalF32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false,
) extends ReciprocalGeneric(RealSpec.Float32Spec,2,8,2,nStage,enableRangeCheck, enablePolynomialRounding) {
}
