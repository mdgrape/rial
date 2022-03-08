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
) extends Module with DebugControlSlave {

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

  if( ! xSpec.disableSubnormal || ! zSpec.disableSubnormal) {
    throw new RuntimeException(f"TODO (${this.getClass.getName}): subnormal is not supported")
    // subnormal numbers are considered to be zero, currently
  }

  if (xSpec == zSpec) {
    for (i <- 0 until len) {
      val (xsgn, xex,  xman) = FloatChiselUtil.decompose(xSpec, io.x(i))
      val z0 = (~xsgn) ## xex ## xman
      io.z(i) := ShiftRegister(z0, nStage)
    }
  } else {
    for (i <- 0 until len) {
      val (xsgn, xex,  xman) = FloatChiselUtil.decomposeWithLeading1(xSpec, io.x(i))
      val (xex0, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, io.x(i))
      val xzero = !(io.x(i)(xSpec.W-2, 0).orR.asBool)
      // TODO: 2^-126 0.0010 -> 2^-129 1.000

      val zsgn = ~xsgn

      val (zmanRound, zex0) = if (xSpec.manW <= zSpec.manW) {// e.g. f32 -> f64
        // no rounding required
        val zmanShift = if(xSpec.manW == zSpec.manW) {
          xman(xSpec.manW-1, 0)
        } else {
          xman(xSpec.manW-1, 0) << (zSpec.manW - xSpec.manW)
        }
        // no rounding, so moreThan2AfterRound never be 1
        val zex0 = xex.pad(exW) - xSpec.exBias.U(exW.W) + zSpec.exBias.U(exW.W)
        (zmanShift, zex0)

      } else { // e.g. f64 -> f32

        // rounding required
        val zman0  = xman(xSpec.manW, xSpec.manW - zSpec.manW) // +leading bit
        val lsb    = zman0(0)
        val round  = xman(xSpec.manW - zSpec.manW - 1)
        val sticky = xman(xSpec.manW - zSpec.manW - 2, 0).orR
        val inc    = FloatChiselUtil.roundIncBySpec(roundSpec, lsb, round, sticky, zsgn(0))
        val zmanRound = zman0 +& inc
        val moreThan2 = zmanRound(zSpec.manW+1) // 1.1111|1 -> 10.0000

        val zex0 = xex.pad(exW) - xSpec.exBias.U(exW.W) + zSpec.exBias.U(exW.W) + moreThan2
        (zmanRound(zSpec.manW-1, 0), zex0)
      }

      val exNeg  = (minEx < 0).B && (zex0(exW-1) =/= 0.U) // -> 0
      val exZero = !zex0.orR.asBool                       // -> subnormal (if supported)
      val exZN   = exZero || exNeg || xzero

      val exInf = !xnan && (if (exW > zSpec.exW) {
        xinf || (zex0(zSpec.exW-1, 0).andR | zex0(exW-1, zSpec.exW).orR)
      } else {
        xinf || zex0(zSpec.exW-1, 0).andR
      })
      val zex = Mux(exZN,         0.U(zSpec.exW),
                Mux(exInf | xnan, Fill(zSpec.exW, 1.U(1.W)),
                                  zex0(zSpec.exW-1, 0)))

      // since subnormal is disabled, we can just set mantissa zero if exponent is zero
      val zman = Mux(exZN | exInf, 0.U(zSpec.exW), zmanRound)

      if(enableDebug) {
        printf("%d -> %d: xsgn = %d, zsgn = %d\n", xSpec.W.U, zSpec.W.U, xsgn, zsgn)
        printf("%d -> %d: xex  = %x, zex  = %x\n", xSpec.W.U, zSpec.W.U, xex , zex )
        printf("%d -> %d: xman = %x, zman = %x\n", xSpec.W.U, zSpec.W.U, xman, zman)
      }
      val z0 = zsgn ## zex ## zman
      io.z(i) := ShiftRegister(z0, nStage)
    }
  }
}
