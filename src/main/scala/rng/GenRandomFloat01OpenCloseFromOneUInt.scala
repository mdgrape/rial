package rial.rng

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._

import rial.arith._
import rial.util._

// generate Float in (0, 1] using one uint by Downey's method.
//
// In case of UInt(32.W) -> BFloat(1, 8, 7):
// - man: 7 bits, downey: 1 bit, ex = 24 bits
// - min value = 2^-25 * 1.0
// - max value = 1.0
//
// In case of UInt(32.W) -> Float32(1, 8, 23):
// - man: 23 bits, downey: 1 bit, ex = 8 bits
// - min value = 2^-9 * 1.0
// - max value = 1.0
//
// In case of UInt(64.W) -> Float32(1, 8, 23):
// - man: 23 bits, downey: 1 bit, ex = 40 bits
// - min value = 2^-41 * 1.0
// - max value = 1.0
//
class GenRandomFloat01OpenCloseFromOneUInt(
  rndW: Int, // random number width
  spec: RealSpec,
  stage: PipelineStageConfig
) extends Module {

  val nStage = stage.total

  val exW  = spec.exW
  val manW = spec.manW

  assert(rndW > manW + 1, "GenRandomFloat01OpenCloseFromOneUInt requires " +
    f"rndW (${rndW}) > manW+1 (${manW+1})")

  val io = IO(new Bundle {
    val rnd = Input(UInt(rndW.W))
    val z   = Output(UInt(spec.W.W))
  })

  val rndMantissa = io.rnd(manW-1, 0)
  val rndDowney   = io.rnd(manW)
  val rndExponent = io.rnd(rndW-1, manW+1)

  val rndExW = rndW - (manW + 1)

  assert(rndExW < spec.exBias-1, "GenRandomFloat01OpenCloseFromOneUInt requires " +
    f"rndExW (${rndExW}) < exBias-1 (${spec.exBias-1})")

  val exCorrection = Mux(rndMantissa === 0.U, rndDowney, 0.U(1.W))

  //       LSB     MSB
  // rndEx = 1xxx..x -> zex = exBias-1 (0.5~1.0) (+ rndDowney (0.5 or 1.0))
  //         01xx..x -> zex = exBias-2 (1/4~1/2)
  //         0000..1 -> zex = exBias-rndExW
  //         0000..0 -> zex = exBias-(rndExW+1)

  val zexP = PriorityEncoder(rndExponent)
  val zex0 = Mux(rndExponent === 0.U, rndExW.U, zexP)
  val zex  = (spec.exBias-1).U(exW.W) - zex0 + exCorrection

//   printf(f"GenRandomFloat01(rndW=${rndW}): rndEx = %%b, zex0 = %%d\n", rndExponent, zex0)

  val zsgn = 0.U(1.W)
  val zman = rndMantissa

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

//   printf(f"GenRandomFloat01(rndW=${rndW}): z = %%b|%%d|%%b\n",zsgn, zex, zman)

  io.z := ShiftRegister(z, nStage)
}
