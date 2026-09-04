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
 *     val rnd = Flipped(Valid(UInt(rndW.W)))
 *     val z   = Valid(UInt(spec.W.W))
 *   })
 *   // ...
 * }
 * }}}
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`    the module input
 *  - `"norm"`  after the exponent/mantissa priority encoders and the
 *              subnormal alignment shift
 *  - `"round"` after the Downey exponent correction carried by the extra
 *              random bit
 *  - `"out"`   the module output
 *
 * [[getStage]] gives the resulting latency, and `z.valid` is `rnd.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // two registers, placed by the module
 * new GenRandomFloat01Close(rndW, spec, 2)
 * // hand-picked positions
 * new GenRandomFloat01Close(rndW, spec,
 *   PipelineStageConfig.specify("norm" -> 1, "round" -> 1))
 * // all the registers at the module output
 * new GenRandomFloat01Close(rndW, spec, PipelineStageConfig.atOut(2))
 * // combinational
 * new GenRandomFloat01Close(rndW, spec, PipelineStageConfig.none)
 * }}}
 *
 * @constructor    create a chisel Module.
 * @param rndW     input random number bit width.
 * @param spec     output floating-point spec.
 * @param stage    pipeline registers of this module.
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
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

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = GenRandomFloat01Close.stagePoints

  val stages = stage.resolve(stagePoints)

  def getParam = { (rndW, spec, stage) }

  val nStage = stages.total

  def getStage = nStage

  val exBias = spec.exBias
  val exW    = spec.exW
  val manW   = spec.manW
  val rndExW = rndW - (manW + 1)

  assert(rndW > manW + 1, f"GenRandomFloat01Close requires rndW (${rndW}) > manW+1 (${manW+1})")

  val io = IO(new Bundle {
    val rnd = Flipped(Valid(UInt(rndW.W)))
    val z   = Valid(UInt(spec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")

  private val cutIn = new SharedPipe("in", io.rnd.valid, stages.at("in"))

  val rndIn = cutIn(io.rnd.bits)

  //----------------------------------------------------------------------
  // S1 : split the random bits and normalise the mantissa

  val rndMantissa = rndIn(manW-1, 0)
  val rndDowney   = rndIn(manW)
  val rndExponent = if(rndExW <= exBias-1) {
    rndIn(rndW-1, manW+1)
  } else {
    rndIn(spec.exBias-1 + manW+1 -1, manW+1)
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

        assert(!cutIn.valid || manAligned(manW) === 1.U,
          "RNG: rndMantissa = %b, manPE = %d, aligned = %b", rndMantissa, manPE, manAligned)

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

  //--------------------------- stage "norm" -----------------------------
  private val cutNorm = new SharedPipe("norm", cutIn.valid, stages.at("norm"))

  val zsgnNorm      = cutNorm(zsgn)
  val zexNorm       = cutNorm(zex)
  val zmanNorm      = cutNorm(zman)
  val rndDowneyNorm = cutNorm(rndDowney)

  //----------------------------------------------------------------------
  // S2 : Downey's exponent correction

  val exCorrection = Mux(zmanNorm === 0.U, rndDowneyNorm, 0.U(1.W))
  val zexCorrected = zexNorm + exCorrection

  // if subnormal numbers are disabled, numbers in range [0, 2^-exMin) will be
  // rounded to zero, so distribution around zero is slightly skewed.
  val zmanCorrected = if(spec.disableSubnormal) {
    Mux(zexCorrected === 0.U, 0.U, zmanNorm)
  } else {
    zmanNorm
  }

  //--------------------------- stage "round" ----------------------------
  private val cutRound = new SharedPipe("round", cutNorm.valid, stages.at("round"))

  val zsgnRound          = cutRound(zsgnNorm)
  val zexCorrectedRound  = cutRound(zexCorrected)
  val zmanCorrectedRound = cutRound(zmanCorrected)

  //----------------------------------------------------------------------
  // S3 : assemble

  val z = Cat(zsgnRound, zexCorrectedRound, zmanCorrectedRound)
  assert(z.getWidth == spec.W)

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutRound.valid, stages.at("out"))

  io.z.bits  := cutOut(z)
  io.z.valid := cutOut.valid
}

/** Pipeline stage points of [[rial.rng.GenRandomFloat01Close]].
 *
 * @see [[rial.util.StagePoint]]
 */
object GenRandomFloat01Close {
  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.rng.GenRandomFloat01Close]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("norm", 3), StagePoint("round", 1), StagePoint("out", 1))
}
