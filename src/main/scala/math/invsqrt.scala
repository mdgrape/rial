//% @file sqrt.scala
//
// square root function
// Copyright (C) Toru Niina RIKEN BDR 2021
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
import rial.arith.FloatChiselUtil

// 1/sqrt(x): floating => floating
//
// 1.m \in [1.0, 2.0)
// =>   sqrt(1.m) \in [1.0, 1.414),     sqrt(2*1.m) \in [1.414, 2.0)
// => 1/sqrt(1.m) \in [1/1.414, 1.0), 1/sqrt(2*1.m) \in [1/2.0, 1/1.414)
// => 2/sqrt(1.m) \in [2/1.414, 2.0), 2/sqrt(2*1.m) \in [1.0, 2/1.414)
//
// 1/sqrt(2^ 5 * 1.m) = 2^-2 * 1/sqrt(2*1.m) = 2^-3 * 2/sqrt(2*1.m)
// 1/sqrt(2^ 4 * 1.m) = 2^-2 * 1/sqrt(1.m)   = 2^-3 * 2/sqrt(1.m)
// 1/sqrt(2^-4 * 1.m) = 2^ 2 * 1/sqrt(1.m)   = 2^ 1 * 2/sqrt(1.m)
// 1/sqrt(2^-5 * 1.m) = 2^ 3 * 1/sqrt(2*1.m) = 2^ 2 * 2/sqrt(2*1.m)
class InvSqrtGeneric(
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

  val io = IO(iodef = new Bundle {
    val x = Input (UInt(spec.W.W))
    val z = Output(UInt(spec.W.W))
  })

  val manW   = spec.manW
  val exW    = spec.exW
  val exBias = spec.exBias

  val order =
    if (adrW>=manW) { // all the mantissa bits are used as a table address
      if (nOrder != 0)
        println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
      0
    } else {
      nOrder
    }

  val bp  = manW+extraBits // for rounding

  val (xsgn,  xex,  xman) = FloatChiselUtil.decompose(spec, io.x)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(spec, io.x)
  val xneg = if(spec.disableSign) {false.B} else {xsgn === 1.U(1.W)}

//   printf("x = %d|%d|%d\n", xsgn, xex, xman)
//   printf("xzero = %d, xinf = %d, xnan = %d\n", xzero, xinf, xnan)

  // 2^ExMin --invsqrt-> 2^-ExMin/2, does not overflow
  // 2^ExMax --invsqrt-> 2^-ExMax/2, does not underflow

  val znan  = xnan
  val zinf  = xzero || (xneg && !xzero)
  val zzero = xinf

  val xExNobias = xex - exBias.U
  val xExHalf   = xExNobias(exW-1) ## xExNobias(exW-1, 1) // arith right shift
  val zEx0      = ~xExHalf // -(xex>>1)-1 = ~(xex>>1)+1-1 = ~(xex>>1)
  val zEx  = Mux(zinf || znan, maskU(exW), Mux(zzero, 0.U(exW.W), zEx0 + exBias.U))
  val zSgn = 0.U(1.W)

  // --------------------------------------------------------------------------
  // mantissa calculation

  val man = io.x(manW, 0) // include LSB of x.ex

  val z0 = Wire(UInt(spec.W.W))

  if (order<=0) { // generate fixed-sized table. no calculation needed.
    val adr = man

    val tbl = VecInit( (0L to 1L<<(adrW+1)).map(
      n => {
        val x = if (n < (1L<<adrW)) {
          (n.toDouble / (1L<<(adrW+1))) * 4.0 + 2.0 // 0.0~0.499 -> 2.0~3.999
        } else {
          (n.toDouble / (1L<<(adrW+1))) * 2.0       // 0.5~0.999 -> 1.0~1.999
        }
        val y = round((2.0 / math.sqrt(x)-1.0) * (1L<<manW))
        if (y >= (1L<<manW)) {
          println("WARNING: mantissa reaches to 2")
          maskL(manW).U(manW.W)
        } else if(y < 0.0) { // not used, actually
          0.U(manW.W)
        } else {
          y.U(manW.W)
        }
      })
    )
    val zman = tbl(adr)
    val zeroFlush = zinf || zzero || znan
    z0 := zSgn ## zEx ## Mux(zeroFlush, znan ## 0.U((manW-1).W), zman)
  } else {
    // Fractional part by Polynomial
    val adr = man((manW+1)-1, (manW+1)-(adrW+1))
    val d   = Cat(~man((manW+1)-(adrW+1)-1), man((manW+1)-(adrW+1)-2,0)).asSInt

    val tableD = new FuncTableDouble(
      x => if(x<0.5) { (2.0/sqrt(x*4.0+2.0))-1.0 } else { (2.0/sqrt(x*2.0))-1.0 }, order )
    tableD.addRange(0.0, 1.0, 1<<(adrW+1)) // gen table, dividing the range [0,1]
    val tableI = new FuncTableInt( tableD, bp ) // convert float table

    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    // signMode=0 : always include sign bit
    // signMode=1 : 2's complement and no sign bit (if possible)
    // signMode=2 : absolute and no sign bit (if possible)

    val coeff  = getSlices(coeffTable(adr), coeffWidth)
    val coeffS = coeff.map( x => x.asSInt )

    def hornerC( c: SInt, z: SInt, dx: SInt, enableRounding: Boolean = false ) : SInt = {
      val zw   = z.getWidth
      val dxBp = dx.getWidth-1
      val dxw  = min(zw, dxBp)
      val dxl  = dx(dxBp, dxBp-dxw).asSInt // ?
      val prod = dxl * z
      val prod_sft = dropLSB(dxw, prod)
      val sum = c +& prod_sft // Extend for safe
//       printf("----------------------------------------\n")
//       printf("horner  : ci + dx * z; z is the current sum\n")
//       printf("c       : %b(%d), W=%d\n", c              , c              , c              .getWidth.U)
//       printf("zw      : %b(%d), W=%d\n", zw      .asUInt, zw      .asUInt, zw      .asUInt.getWidth.U)
//       printf("dxBp    : %b(%d), W=%d\n", dxBp    .asUInt, dxBp    .asUInt, dxBp    .asUInt.getWidth.U)
//       printf("dxw     : %b(%d), W=%d\n", dxw     .asUInt, dxw     .asUInt, dxw     .asUInt.getWidth.U)
//       printf("dx      : %b(%d), W=%d\n", dx             , dx             , dx             .getWidth.U)
//       printf("dxl     : %b(%d), W=%d\n", dxl            , dxl            , dxl            .getWidth.U)
//       printf("z       : %b(%d), W=%d\n", z              , z              , z              .getWidth.U)
//       printf("prod    : %b(%d), W=%d\n", prod           , prod           , prod           .getWidth.U)
//       printf("prod_sft: %b(%d), W=%d\n", prod_sft       , prod_sft       , prod_sft       .getWidth.U)
//       printf("sum     : %b(%d), W=%d\n", sum            , sum            , sum            .getWidth.U)
      sum
    }

    def polynomialEvaluationC( c : Seq[SInt], dx : SInt, enableRounding: Boolean = false ) = {
      val res = c.init.foldRight(c.last)( (c,z) => hornerC( c, z, dx, enableRounding ) )
      res
    }

    val s0 = polynomialEvaluationC( coeffS,  d, enablePolynomialRounding ).asUInt

//     printf("Odd : %b(%d), width=%d\n", zodd , zodd , zodd .getWidth.U)
//     printf("Even: %b(%d), width=%d\n", zeven, zeven, zeven.getWidth.U)
//     printf("xEx : %b(%d), width=%d\n", xExNobias, xExNobias, xExNobias.getWidth.U)
//     printf("res0: %b(%d), width=%d\n", res0,  res0, res0.getWidth.U)
//     printf("rres: %b(%d), width=%d\n", rres,  rres, rres.getWidth.U)
//     printf("s0  : %b(%d), width=%d\n", s0,    s0, s0.getWidth.U)

    val zman  = dropLSB(extraBits, s0) + s0(extraBits-1)
    val polynomialOvf = zman(manW)
    val zeroFlush = zinf || zzero || znan
    val zman_checkovf = Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman(manW-1,0))

//     printf("zman: %b(%d), width=%d\n", zman , zman , zman .getWidth.U)
//     printf("xExNobias(0) = %b, zman.head(1) = %b, zman.head(2) = %b\n", xExNobias(0), zman.head(1), zman.head(2))
//     printf("z0: %b|%b(%d)|%b(%d)\n", zSgn, zEx, zEx, zman_checkovf, zman_checkovf)

    z0 := zSgn ## zEx ## Mux(zeroFlush, znan ## 0.U((manW-1).W), zman_checkovf)
  }

  io.z   := ShiftRegister(z0, nStage)
  //printf("x=%x z=%x\n", io.x, io.z)
}

class InvSqrtFP32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false)
  extends InvSqrtGeneric(RealSpec.Float32Spec, 2, 8, 2, nStage,
                         enableRangeCheck, enablePolynomialRounding) {
}
