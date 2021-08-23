package rial.simd

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import rial.arith._
import rial.util._
import rial.util.ScalaUtil._

/**
  * Compute negation of packed floating point numbers.
  *
  */
class NegPackedFPGeneric(
  len : Int, // length of SIMD pack
  xSpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig,
  val enableDebug : Boolean = false
) extends MultiIOModule with DebugControlSlave {

  val nStage = stage.total

  def getParam() = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage() = nStage

  if( xSpec.disableSign || zSpec.disableSign) {
    throw new RuntimeException(f"ERROR (${this.getClass.getName}): sign is disabled")
  }

  val io = IO(new Bundle {
    val x = Input (Vec(len, UInt(xSpec.W.W)))
    val z = Output(Vec(len, UInt(zSpec.W.W)))
  })

  // prepare parameters
  val dBias = zSpec.exBias - xSpec.exBias  // xEx - xBias == zEx - zBias
  val maxEx = maskI(xSpec.exW) - 1 + dBias // maximum non-inf number
  val minEx = 1                    + dBias // minimum normalized number
  val exW0  = max(log2Up(maxEx + 1), log2Up(abs(minEx)+1)) // required width
  val exW   = if (minEx < 0) exW0+1 else exW0 // if it has negative value,
                                              // the range of value is doubled
  if(enableDebug) {
    println(xSpec.W, " -> ", zSpec.W)
    println("dBias = ",  dBias)
    println("maxEx = ",  maxEx)
    println("minEx = ",  minEx)
    println("exW0  = ",  exW0)
    println("exW   = ",  exW)
  }

  for (i <- 0 until len) {
    val (xsgn, xex,  xman) = FloatChiselUtil.decompose (xSpec, io.x(i))
    val (xex0, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, io.x(i))
    val xzero = !(io.x(i)(xSpec.W-2, 0).orR.asBool)

    val zsgn = ~xsgn // negate

    val zex0   = xex + dBias.S(exW)
    val exNeg  = (minEx < 0).B && (zex0(exW-1) =/= 0.U)
    val exZero = !zex0.orR.asBool
    val exZN   = exZero || exNeg || xzero

    // x is already inf ? or conversion result exceeds range of zSpec ?
    val exInf = !xnan && (if (exW > zSpec.manW) {
      xinf || (zex0(zSpec.exW-1, 0).andR | zex0(exW-1, zSpec.exW).orR)
    } else {
      xinf || zex0(zSpec.exW-1, 0).andR
    })

    val zex = Mux(exZN,         0.U(zSpec.exW),
              Mux(exInf | xnan, Fill(zSpec.exW, 1.U(1.W)),
                                zex0(zSpec.exW-1, 0)))

    val zman_rounded = if(xSpec.manW > zSpec.manW) {
      val zman0  = xman(xSpec.manW - 1, xSpec.manW - zSpec.manW)
      val lsb    = zman0(0)
      val round  = xman(xSpec.manW - zSpec.manW - 1)
      val sticky = xman(xSpec.manW - zSpec.manW - 2, 0).orR
      val inc   = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky, zsgn(0))
      zman0 +& inc
    } else if (xSpec.manW < zSpec.manW) {
      xman(xSpec.manW - 1, 0) ## 0.U(zSpec.manW - xSpec.manW)
    } else {
      xman(xSpec.manW - 1, 0)
    }

    val zman = Mux(xzero || exInf, 0.U(zSpec.manW.W), zman_rounded)

    if(enableDebug) {
      printf("xsgn = %d, zsgn = %d\n", xsgn, zsgn)
      printf("xexp = %x, zexp = %x\n", xex , zex )
      printf("xman = %x, zman = %x\n", xman, zman)
    }

    val z0 = zsgn ## zex ## zman
    io.z(i) := ShiftRegister(z0, nStage)
  }
}
