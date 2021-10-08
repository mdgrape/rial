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

// sqrt(x): floating => floating
// x = 2^e * 1.m
// sqrt(x) = 2^( e   /2) * sqrt(1.m)           ... e % 2 == 0
//         = 2^((e-1)/2) * sqrt(1.m) * sqrt(2) ... e % 2 == 1
class SqrtGeneric(
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
//   printf("x = %d|%d(%d)|%d\n", xsgn, xex, xex-spec.exBias.U, xman)
  // special values

  // sqrt never becomes inf nor underflows. 100 -> 10, 1/100 -> 1/10.
  // sqrt(-0) is not NaN, it is zero.
  val znan  = xnan || (xneg && !xzero)
  val zinf  = xinf
  val zzero = xzero

  // 8   = 2^3   => 2sqrt(2)    = 2    * 1.414:  3(0011) ->  1(0001)
  // 1/8 = 1/2^3 => 1/2sqrt(2)  = 2^-2 * 1.414: -3(1101) -> -2(1110)
  // 32  = 2^5   => 2^2sqrt(2)  = 2^2  * 1.414:  5(0101) ->  2(0010)
  // 1/32= 1/2^5 => 1/2^2sqrt(2)= 2^-3 * 1.414: -5(1011) -> -3(1101)
  val xExNobias = xex - exBias.U
  val zex0 = xExNobias(exW-1) ## xExNobias(exW-1, 1) // arith right shift
  val zEx  = Mux(zinf || znan, maskU(exW), Mux(zzero, 0.U(exW.W), zex0 + exBias.U))
  val zSgn = 0.U(1.W) // always positive.

  val z0 = Wire(UInt(spec.W.W))

  if (order<=0) { // generate fixed-sized table. no calculation needed.
    val adr = xman(manW-1, manW-adrW) // take MSBs as an address

    // sqrt(1.m) (1.x * 2^4 -> sqrt(1.x) * 2^2)
    val tbl_even = VecInit( (0L to (1L<<adrW)-1L).map(
      n => {
        val x = 1.0 + n.toDouble / (1L<<adrW)  // value from 1 to 2
        val y = round((math.sqrt(x)-1.0) * (1L<<manW))
        if (y >= (1L<<manW)) {
          println("WARNING: mantissa reaches to 2")
          maskL(manW).U(manW.W)
        } else {
          y.U(manW.W)
        }
      })
    )
    // sqrt(1.m) * 1.414 (1.x * 2^5 -> sqrt(1.x) * sqrt(2) * 2^2)
    val tbl_odd = VecInit( (0L to (1L<<adrW)-1L).map(
      n => {
        val x = 1.0 + n.toDouble / (1L<<adrW)  // value from 1 to 2
        val y = round((math.sqrt(x * 2.0)-1.0) * (1L<<manW))
        if (y >= (1L<<manW)) {
          println("WARNING: mantissa reaches to 2")
          maskL(manW).U(manW.W)
        } else {
          y.U(manW.W)
        }
      })
    )
    val zman = Mux(xExNobias(0) === 1.U(1.W), tbl_odd(adr), tbl_even(adr))
    val zeroFlush = zinf || zzero || znan
    z0 := zSgn ## zEx ## Mux(zeroFlush, znan ## 0.U((manW-1).W), zman)
  } else {
    // Fractional part by Polynomial
    val adr = ~xman(manW-1, manW-adrW)
    val d   = Cat(xman(manW-adrW-1),~xman(manW-adrW-2,0)).asSInt
//     printf("xman = %b(%d)\n", xman, xman)
//     printf("d    = %b(%d)\n", d, d)
//     printf("adr  = %b(%d)\n", adr, adr)

    val eps = pow(2.0, -manW)
    val tableDeven = new FuncTableDouble(x => 2.0 - math.sqrt(2.0 -     (x+eps)), order)
    val tableDodd  = new FuncTableDouble(x => 2.0 - math.sqrt(4.0 - 2.0*(x+eps)), order)
    tableDeven.addRange(0.0, 1.0, 1<<adrW) // gen table, dividing the range [0,1]
    tableDodd .addRange(0.0, 1.0, 1<<adrW) // into 1<<adrW subregions
    val tableIeven = new FuncTableInt( tableDeven, bp ) // convert float table
    val tableIodd  = new FuncTableInt( tableDodd,  bp ) // into mantissa-Int
    val (coeffTableEven, coeffWidthEven) = tableIeven.getVectorUnified(/*sign mode =*/1)
    val (coeffTableOdd,  coeffWidthOdd ) = tableIodd .getVectorUnified(1)
    // signMode=0 : always include sign bit
    // signMode=1 : 2's complement and no sign bit (if possible)
    // signMode=2 : absolute and no sign bit (if possible)

    val coeffEven  = getSlices(coeffTableEven(adr), coeffWidthEven)
    val coeffOdd   = getSlices(coeffTableOdd (adr), coeffWidthOdd)
    val coeffSEven = coeffEven.map( x => x.zext )
    val coeffSOdd  = coeffOdd .map( x => x.zext )

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

    val zodd  = polynomialEvaluationC( coeffSOdd,  d, enablePolynomialRounding ).asUInt
    val zeven = polynomialEvaluationC( coeffSEven, d, enablePolynomialRounding ).asUInt

    val res0 = Mux(xExNobias(0) === 1.U(1.W), zodd, zeven) // 2-sqrt(...)
    val rres = res0 - (2.U << (manW + extraBits))          //  -sqrt(...)
    val s0   = ~rres + 1.U - (1.U << (manW + extraBits))   //   sqrt(...)

//     printf("Odd : %b(%d), width=%d\n", zodd , zodd , zodd .getWidth.U)
//     printf("Even: %b(%d), width=%d\n", zeven, zeven, zeven.getWidth.U)
//     printf("xEx : %b(%d), width=%d\n", xExNobias, xExNobias, xExNobias.getWidth.U)
//     printf("res0: %b(%d), width=%d\n", res0,  res0, res0.getWidth.U)
//     printf("rres: %b(%d), width=%d\n", rres,  rres, rres.getWidth.U)
//     printf("s0  : %b(%d), width=%d\n", s0,    s0, s0.getWidth.U)

    val zman  = dropLSB(extraBits, s0) + s0(extraBits-1)
    val polynomialOvf = zman.head(2) === 1.U
    val polynomialUdf = zman.head(1) === 1.U
    val zeroFlush = zinf || zzero || polynomialUdf || znan
    val zman_checkovf = Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman(manW-1,0))

//     printf("zman: %b(%d), width=%d\n", zman , zman , zman .getWidth.U)
//     printf("xExNobias(0) = %b, zman.head(1) = %b, zman.head(2) = %b\n", xExNobias(0), zman.head(1), zman.head(2))
//     printf("z0: %b|%b(%d)|%b(%d)\n", zSgn, zEx, zEx, zman_checkovf, zman_checkovf)

    z0 := zSgn ## zEx ## Mux(zeroFlush, znan ## 0.U((manW-1).W), zman_checkovf)
  }

  io.z   := ShiftRegister(z0, nStage)
  //printf("x=%x z=%x\n", io.x, io.z)
}

class SqrtFP32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false)
  extends SqrtGeneric(RealSpec.Float32Spec, 2, 8, 2, nStage,
                      enableRangeCheck, enablePolynomialRounding) {
}
