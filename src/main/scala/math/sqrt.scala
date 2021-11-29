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
// - if x < 0, returns 0.
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

  // sqrt never becomes inf nor underflows. 100 -> 10, 1/100 -> 1/10.
  // sqrt(-0) is not NaN, it is zero.
  val znan  = xnan
  val zinf  = xinf
  val zzero = xzero || xneg

  val xExNobias = xex - exBias.U
  val zex0 = xExNobias(exW-1) ## xExNobias(exW-1, 1) // arith right shift
  val zEx  = Mux(zinf || znan, maskU(exW),
             Mux(zzero,        0.U(exW.W),
                               zex0 + exBias.U))

  val zSgn = 0.U(1.W) // always positive.

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
        val y = round((math.sqrt(x)-1.0) * (1L<<manW))
        if (y >= (1L<<manW)) {
          maskL(manW).U(manW.W)
        } else if (y <= 0.0) {
          0.U(manW.W)
        } else {
          y.U(manW.W)
        }
      })
    )
    val zman = tbl(adr)
    val zeroFlush = zinf || zzero || znan
    z0 := Cat(zSgn, zEx, Mux(zeroFlush, znan ## 0.U((manW-1).W), zman))
  } else {
    // Fractional part by Polynomial
    val adr = man((manW+1)-1, (manW+1)-(adrW+1))
    val d   = Cat(~man((manW+1)-(adrW+1)-1), man((manW+1)-(adrW+1)-2,0)).asSInt
    assert(adr.getWidth == adrW+1)
    assert(d  .getWidth == manW-adrW)
    printf("man = %b(%d)\n", man, man)
    printf("adr = %b(%d)\n", adr, adr)
    printf("dx = %b(%d)\n", d, d)

    val tableD = new FuncTableDouble(
      x => if(x<0.5) { sqrt(x*4.0+2.0)-1.0 } else { sqrt(x*2.0)-1.0 }, order )

    tableD.addRange(0.0, 1.0, 1<<(adrW+1))
    val tableI = new FuncTableInt( tableD, bp )

    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)

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
      sum
    }

    def polynomialEvaluationC( c : Seq[SInt], dx : SInt, enableRounding: Boolean = false ) = {
      val res = c.init.foldRight(c.last)( (coeff,z) => hornerC( coeff, z, dx, enableRounding ) )
      res
    }

    val res = polynomialEvaluationC( coeffS,  d, enablePolynomialRounding ).asUInt

    val zman = dropLSB(extraBits, res) + res(extraBits-1)
    val polynomialOvf = zman(manW)
    val zeroFlush     = zinf || zzero || znan
    val zman_checkovf = Mux(zeroFlush,     Cat(znan, 0.U((manW-1).W)),
                        Mux(polynomialOvf, Fill(manW,1.U(1.W)), zman(manW-1,0)))

//     printf("zman: %b(%d), width=%d\n", zman , zman , zman .getWidth.U)
//     printf("xExNobias(0) = %b, zman.head(1) = %b, zman.head(2) = %b\n", xExNobias(0), zman.head(1), zman.head(2))
//     printf("z0: %b|%b(%d)|%b(%d)\n", zSgn, zEx, zEx, zman_checkovf, zman_checkovf)

    z0 := Cat(zSgn, zEx, zman_checkovf)
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
