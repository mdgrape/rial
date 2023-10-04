//% @file tableEval.scala
//
// polynomial approx using coefficient table
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._

import rial.arith.RealSpec
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._

/** An I/O bundle for [[rial.math.PolynomialEval]] Module and table coeff module
 *  of functions.
 *  The direction is Input. To use it in a function, flip it.
 *
 *  @constructor create a new TableCoeffInput.
 *  @param cbit bitwidths of coefficients.
 */
class TableCoeffInput( val cbit: Seq[Int] ) extends Bundle {
  val cs = Input(MixedVec(cbit.map{w => UInt(w.W)}))
}

/** Config class for [[rial.math.PolynomialEval]].
 *
 *  @constructor create a new PolynomialSpec.
 *  @param spec the spec of floating point number.
 *  @param nOrder the order of polynomial.
 *  @param adrW the bit width of table address.
 *  @param dxW0 the input bit width of polynomial.
 */
class PolynomialSpec (
  val manW            : Int,
  val nOrder          : Int,
  val adrW            : Int,
  val extraBits       : Int,
      dxW0            : Option[Int] = None,
  val enableRangeCheck: Boolean = true,
  val enableRounding  : Boolean = false,
) {
  def order : Int = {if(adrW == manW) {0} else {nOrder}}
  def fracW : Int = {manW + extraBits}
  def dxW   : Int = {dxW0.getOrElse(manW - adrW)}
}

/** Evaluates polynomial function.
 *
 *  @constructor create a new PolynomialEval module.
 *  @param spec the spec of floating point number.
 *  @param polySpec the spec of polynomial.
 *  @param cbit bitwidths of polynomial coefficients.
 *  @param stage the number of pipeline stages inside this module.
 */
class PolynomialEval(
  val spec:           RealSpec,
  val polySpec:       PolynomialSpec,
  val cbit:           Seq[Int],
  val stage:          PipelineStageConfig,
) extends Module {

  val nStage = stage.total

  val manW  = spec.manW
  val dxW   = polySpec.dxW
  val fracW = polySpec.fracW
  val order = polySpec.order

  val io = IO(new Bundle {
    val coeffs = new TableCoeffInput(cbit)
    val dx     = if (order != 0) { Some(Input(UInt(dxW.W))) } else { None }
    val result = Output(UInt(fracW.W))
  })

  val res = Wire(UInt(fracW.W))

  if(order == 0) {

    res := io.coeffs.cs(0)

  } else {

    val dxS = io.dx.get.asSInt

    def hornerC( c: SInt, z: SInt, dx: SInt ) : SInt = {
      val zw   = z.getWidth
      val dxBp = dx.getWidth-1
      val dxw  = min(zw, dxBp)
      val dxl  = dx(dxBp, dxBp-dxw).asSInt
      val prod = dxl * z
      val prod_sft = dropLSB(dxw, prod)
//       val sum = c +& prod_sft // Extend for safe
      val sum = c + prod_sft

//       printf("cir: z   = %b\n", z)
//       printf("cir: dx  = %b\n", dx)
//       printf("cir: dxBp= %d\n", dxBp.U)
//       printf("cir: zw  = %d\n", zw.U)
//       printf("cir: dxw = %d\n", dxw.U)
//       printf("cir: dxl = %b\n", dxl)
//       printf("cir: c   = %b\n", c)
//       printf("cir: sum = %b\n", sum)

      sum
    }

    val coeffS = io.coeffs.cs.map( x => x.asSInt )

    val resS = coeffS.init.foldRight(coeffS.last)(
      (c,z) => hornerC( c, z, dxS )
    )

    // Since some of the math funcs use the extraBits as a part of mantissa
    // to normalize the result, we cannot round it here.
    res := resS.asUInt
  }
  io.result := ShiftRegister(res, nStage)
}

