package rial.funcs

import rial.arith.RealSpec
import rial.approx.PolynomialCoefficients
import rial.util.PipelineStageConfig

import chisel3._
import chisel3.util._

class RoundingPostProcInfo(val spec: RealSpec) extends Bundle {
  val zsgn        = UInt(1.W)
  val zex         = UInt(spec.exW.W)
  val zinf        = Bool()
  val znan        = Bool()
  val zIsNonTable = Bool()
}

// It just rounds the polynomial result.
// sqrt, invsqrt, reciprocal requires the same postprocess.
//
class RoundingPostProcess(
  val spec:  RealSpec, // Input / Output floating spec
  val coeff: PolynomialCoefficients,
  val stage: PipelineStageConfig,
) extends Module {

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val adrW   = coeff.bitW.addrW
  val fracW  = coeff.bitW.outBitW
  val order  = coeff.bitW.nOrder

  val extraBits = fracW - manW

  val nStage = stage.total

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
      val info = new RoundingPostProcInfo(spec)
      val zres = UInt(fracW.W)
    }))
    val z = Valid(UInt(spec.W.W))
  })

  val zsgn = io.in.bits.info.zsgn
  val zex  = io.in.bits.info.zex
  val znan = io.in.bits.info.znan
  val zIsNonTable  = io.in.bits.info.zIsNonTable // === znan || zinf || zzero
  val zmanNonTable = Cat(znan, 0.U((manW-1).W))

  val zmanRounded = WireDefault(0.U(manW.W))
  if(extraBits == 0) {
    zmanRounded := io.in.bits.zres
  } else {
    val zman0 = (io.in.bits.zres >> extraBits) +& io.in.bits.zres(extraBits-1)
    val polynomialOvf = zman0(manW)
    zmanRounded := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW-1,0))
  }

  val zman = Mux(zIsNonTable, zmanNonTable, zmanRounded)
  val z = Cat(zsgn, zex, zman)

  io.z := Pipe(io.in.valid, z, nStage)
}
