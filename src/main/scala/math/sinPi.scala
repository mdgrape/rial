//% @file sinPi.scala
//
// square root function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import chisel3._
import chisel3.util._

import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.FloatChiselUtil
import rial.arith._

//
// sinPi(x): floating -> floating = sin(pi * x)
//
class SinPiGeneric(
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

  val (xsgn_, xex_, xman_) = FloatChiselUtil.decompose(spec, io.x)
  val (xzero, xinf, xnan)  = FloatChiselUtil.checkValue(spec, io.x)
  val xExNobias = xex_.zext - exBias.S((exW + 1).W)

//   printf("x = %d|%d(%d)|%d\n", xsgn_, xex_, xExNobias, xman_)

  val out_of_bounds = (1.S <= xExNobias) && (xman_ =/= 0.U) // 2 < |x|
  val znan  = xnan || xinf || out_of_bounds

  //         y
  //         ^
  //         | .--.
  // _.______|'____'.______.__ x
  // -1'.__.'|      1'.__.' 2
  //         |
  //           '-'--'------'
  //            | -1  ex=0
  //            -2

  // --------------------------------------------------------------------------
  // convert everything into [0, 0.5) (ex = [-inf to -2])

  val manExceedsHalf = xman_(manW-1) === 1.U
  val halfpos  = xman_               // width = manW
  val halfneg  = (1<<manW).U - xman_ // width = manW+1
  val shiftpos = PriorityEncoder(Reverse(halfpos)) + 1.U
  val shiftneg = PriorityEncoder(Reverse(halfneg)) // already +1-ed

  val halfEx0  = Mux(manExceedsHalf, halfneg, halfpos)
  val shiftEx0 = Mux(manExceedsHalf, shiftneg, shiftpos)

  val xman = MuxCase(xman_, Seq(
    (xExNobias ===  0.S) -> ((halfEx0 << shiftEx0) - (1 << manW).U),
    (xExNobias === -1.S) -> ((halfneg << shiftneg) - (1 << manW).U)
  ))
  val xex  = Mux(xExNobias === -1.S, (-1.S((exW+1).W) - shiftneg.zext),
             Mux(xExNobias ===  0.S, Mux(halfEx0 === 0.U, -exBias.S((exW+1).W), -shiftEx0.zext),
                                     xExNobias))
//   printf("circ xex  = %b(%d)\n", xex,  xex)
//   printf("circ xman = %b(%d)\n", xman, xman)

  val zzero = ((0.S <= xExNobias) && xman_ === 0.U) || // integer input
              (xex === -exBias.S  && xman  === 0.U)    // converted into 0
  val zone  = (xex === -1.S       && xman  === 0.U)    // pi/2

  // RealGeneric does not allow -0
  val zSgn = !zzero && ((xsgn_ === 1.U) || (xExNobias === 0.S))

  // --------------------------------------------------------------------------
  // linear approximation around zero

  val linearThreshold = (-math.ceil((manW + 1) / 2)).toInt // -12, if FP32
  val tableExRange    = -2 - linearThreshold + 1
  val pi = new RealGeneric(spec, Pi)

  val isLinear = (xex < linearThreshold.S((exW+1).W))

  val linearProdEx        = (pi.ex-exBias).S((exW+1).W) + xex
  val linearProdMan       = (pi.man + (1<<manW)).toLong.U((manW+1).W) * (xman + (1<<manW).U((manW+1).W))
  val linearProdbp        = manW + manW
  val linearProdMoreThan2 = linearProdMan(linearProdbp+1)
  val linearProdRoundBits = linearProdbp - manW

  val linearProdSticky    = linearProdMan(linearProdRoundBits-2, 0).orR |
                           (linearProdMoreThan2 & linearProdMan(linearProdRoundBits-1))
  val linearProdRound     = Mux(linearProdMoreThan2, linearProdMan(linearProdRoundBits),
                                                     linearProdMan(linearProdRoundBits-1))
  val linearProdShift     = Mux(linearProdMoreThan2, linearProdMan(linearProdbp, linearProdRoundBits+1),
                                                     linearProdMan(linearProdbp-1, linearProdRoundBits))
  val linearProdLSB       = linearProdShift(0)

  val linearProdInc       = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven,
                                linearProdLSB, linearProdRound, linearProdSticky)
  val linearProdRounded   = linearProdShift +& linearProdInc
  val linearProdMoreThan2AfterRound = linearProdRounded(manW)

  val zManLinear = linearProdRounded(manW-1, 0)
  val zExLinear0 = linearProdEx + (linearProdMoreThan2 | linearProdMoreThan2AfterRound).zext + exBias.S((exW+1).W)
  assert(zExLinear0 >= 0.S)
  val zExLinear  = zExLinear0(exW-1, 0)

  // --------------------------------------------------------------------------
  // interpolation using table
  //

  val zExTable  = Wire(UInt(exW.W))
  val zManTable = Wire(UInt(manW.W))

  val exAdrW = log2UpL(tableExRange)
  val exOfs = -xex - 2.S
  val exAdr = MuxCase(exOfs.asUInt()(exAdrW-1, 0), Seq(
    (exOfs < 0.S           ) -> 0.U(exAdrW.W),
    (exOfs > tableExRange.S) -> 0.U(exAdrW.W)
  ))
  if (order<=0) {
    val adr = xman(adrW-1, 0)

    val tbl = VecInit( (-2 to linearThreshold by -1).map( exponent => {
      VecInit((0L to (1L << adrW)).map(
        n => {
          val x = n.toDouble / (1L<<adrW)
          val y = round(scalb(math.sin(Pi * scalb(1.0 + x, exponent)), -exponent-3) * (1L<<bp))
          assert(y < (1L<<bp))
          y.U(bp.W)
        }))
      })
    )
    val zres0 = tbl(exAdr)(adr)
    val zresLessThanHalf = zres0(manW-1) === 0.U
    val zex0 = Mux(zresLessThanHalf, xex + (1+exBias).S, xex + (2+exBias).S)
    val zex  = Mux(zex0 < 0.S, 0.U(exW), zex0(exW-1, 0))

    val zman = Mux(zresLessThanHalf,
                  Cat(zres0, 0.U(2.W))(manW-1, 0) - (1<<manW).U,
                  Cat(zres0, 0.U(1.W))(manW-1, 0) - (1<<manW).U)

    zExTable  := zex
    zManTable := zman
  } else { // second order table
    val adr = ~xman(manW-1, manW-adrW)
    val d   = Cat(xman(manW-adrW-1),~xman(manW-adrW-2,0)).asSInt

    val eps = pow(2.0, -manW)
    val coeffWidth = (-2 to linearThreshold by -1).map(exponent => {
      val tableD = new FuncTableDouble( x => 1.0 - sin(Pi * scalb(2.0 - (x+eps), exponent)), order)
      tableD.addRange(0.0, 1.0, 1<<adrW)
      val tableI = new FuncTableInt( tableD, bp ) // convert float table into int
      val w = tableI.getCBitWidth(/*sign mode = */1)
      println(f"//  width : "+w.mkString("", ", ", ""))
      w
    }).reduce( (lhs, rhs) => {
      lhs.zip(rhs).map( x => max(x._1, x._2))
    })
    println(f"//  maxwidth : "+coeffWidth.mkString("", ", ", ""))

    val coeffTables = VecInit((-2 to linearThreshold by -1).map(exponent => {
      val tableD = new FuncTableDouble( x => 1.0 - sin(Pi * scalb(2.0 - (x+eps), exponent)), order)
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
      val dxl  = dx(dxBp, dxBp-dxw).asSInt // ?
      val prod = dxl * z
      val prod_sft = dropLSB(dxw, prod)
      val sum = c +& prod_sft // Extend for safe
      sum
    }

    def polynomialEvaluationC( c : Seq[SInt], dx : SInt, enableRounding: Boolean = false ) = {
      val res = c.init.foldRight(c.last)( (c,z) => hornerC( c, z, dx, enableRounding ) )
      res
    }

    val res0  = polynomialEvaluationC( coeffS,  d, enablePolynomialRounding ).asUInt
    val res   = ((1 << bp).U - res0)(bp-1, 0)
    val shift = PriorityEncoder(Reverse(res)) + 1.U
    val s0    = (res << shift) - (1 << bp).U
    val zman  = dropLSB(extraBits, s0) +& s0(extraBits-1)

    zExTable  := exBias.U - shift
    zManTable := Mux(zman(manW) === 1.U, Fill(manW, 1.U(1.W)), zman) // check overflow
  }
//   printf("zExTable  = %b(%d)\n", zExTable , zExTable)
//   printf("zManTable = %b(%d)\n", zManTable, zManTable)

  val zeroFlush = zzero || zone || znan

  val zMan = Mux(zeroFlush, Cat(znan, Fill(manW-1, 0.U(1.W))),
             Mux(isLinear,  zManLinear, zManTable))
  val zEx  = Mux(znan,  Fill(exW, 1.U(1.W)),
             Mux(zzero, Fill(exW, 0.U(1.W)),
             Mux(zone,  exBias.U(exW.W),
             Mux(isLinear, zExLinear, zExTable))))
  val z0 = zSgn ## zEx ## zMan
  io.z   := ShiftRegister(z0, nStage)
}

class SinPiFP32(
  nStage: PipelineStageConfig = new PipelineStageConfig,
  enableRangeCheck : Boolean = true,
  enablePolynomialRounding : Boolean = false)
  extends SinPiGeneric(RealSpec.Float32Spec, 2, 8, 6, nStage,
                      enableRangeCheck, enablePolynomialRounding) {
}
