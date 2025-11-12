package rial.rng

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._

import rial.arith._
import rial.util._

/** generate Float in [0, 1] using one uint by Downey's method.
 *
 * {{{
 * class GenRandomFloat01Close(...) extends Module {
 *   val io = IO(new Bundle {
 *     val rnd = Input(UInt(rndW.W))
 *     val z   = Output(UInt(spec.W.W))
 *   })
 *   // ...
 * }
 * }}}
 *
 * @constructor    create a chisel Module.
 * @param rndW     input random number bit width.
 * @param spec     output floating-point spec.
 * @param stage    pipeline stages.
 */
class GenRandomFloat01Close(
  rndW: Int, // random number width
  spec: RealSpec,
  stage: PipelineStageConfig
) extends Module {

  // In case of UInt(32.W) -> BFloat(1, 8, 7):
  // -> man: 7 bits, downey: 1 bit, ex = 24 bits
  // 
  // In case of UInt(32.W) -> Float32(1, 8, 23):
  // -> man: 23 bits, downey: 1 bit, ex = 8 bits
  // 
  // In case of UInt(64.W) -> Float32(1, 8, 23):
  // -> man: 23 bits, downey: 1 bit, ex = 40 bits
  // 
  // In case of UInt(128.W) -> Float32(1, 8, 23):
  // -> man: 23 bits, downey: 1 bit, ex = 104 bits
  // 
  // rndEx = 1xxx..x -> zex = exBias-1 (0.5~1.0) (+ rndDowney (0.5 or 1.0))
  //         01xx..x -> zex = exBias-2 (1/4~1/2)
  //         0000..1 -> zex = exBias-rndExW
  //         0000..0 -> zex = something like a subnormal number.

  val nStage = stage.total

  val exBias = spec.exBias
  val exW    = spec.exW
  val manW   = spec.manW
  val rndExW = rndW - (manW + 1)

  assert(rndW > manW + 1, f"GenRandomFloat01Close requires rndW (${rndW}) > manW+1 (${manW+1})")

  val io = IO(new Bundle {
    val rnd = Input(UInt(rndW.W))
    val z   = Output(UInt(spec.W.W))
  })

  val rndMantissa = io.rnd(manW-1, 0)
  val rndDowney   = io.rnd(manW)
  val rndExponent = if(rndExW <= exBias-1) {
    io.rnd(rndW-1, manW+1)
  } else {
    io.rnd(spec.exBias-1 + manW+1 -1, manW+1)
  }
  assert(rndExponent.getWidth <= exBias-1)

  val zsgn = WireDefault(0.U(1.W))
  val zex  = WireDefault(0.U(exW.W))
  val zman = WireDefault(0.U(manW.W))

  //       LSB     MSB
  // rndEx = 1xxx..x -> zex = exBias-1 (0.5~1.0) (+ rndDowney (0.5 or 1.0))
  //         01xx..x -> zex = exBias-2 (1/4~1/2)
  //         0000..1 -> zex = exBias-rndExW
  //         0000..0 -> zex = exBias-(rndExW+1) ; use man like subnormal

  when(rndExponent === 0.U) { // subnormal!

    if(rndExW < exBias-1) {

      when(rndMantissa === 0.U) {
        zex  := 0.U
        zman := 0.U
      }.otherwise {
        val baseEx     = (exBias - 1 - rndExW).U // > 0
        val manPE      = PriorityEncoder(Reverse(rndMantissa))
        val manAligned = (rndMantissa << (manPE + 1.U))(manW, 0)

        assert(manAligned(manW) === 1.U, "RNG: rndMantissa = %b, manPE = %d, aligned = %b", rndMantissa, manPE, manAligned)

        zex  := Mux(manPE > baseEx, 0.U, baseEx - manPE)
        zman := manAligned(manW-1, 0)
      }
    } else { // rndExW >= exBias-1. we have looooong random bits input.
      zex  := 0.U
      zman := rndMantissa
    }
  }.otherwise { // normal

    val zexPE = PriorityEncoder(rndExponent)
    zex  := (spec.exBias-1).U(exW.W) - zexPE
    zman := rndMantissa
  }

  val exCorrection = Mux(zman === 0.U, rndDowney, 0.U(1.W))
  val zexCorrected = zex + exCorrection

  // if subnormal numbers are disabled, numbers in range [0, 2^-exMin) will be
  // rounded to zero, so distribution around zero is slightly skewed.
  val zmanCorrected = if(spec.disableSubnormal) {
    Mux(zexCorrected === 0.U, 0.U, zman)
  } else {
    zman
  }

  val z = Cat(zsgn, zexCorrected, zmanCorrected)
  assert(z.getWidth == spec.W)

  io.z := ShiftRegister(z, nStage)
}
