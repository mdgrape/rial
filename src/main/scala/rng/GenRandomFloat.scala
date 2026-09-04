//% @file GenRandomFloat.scala
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.rng

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._

import rial.arith._
import rial.util._

/** generate Float in [1, 2).
 *
 * It sets random bits into the mantissa part.
 *
 * Pipeline registers can be enabled at the module boundaries `"in"` and
 * `"out"` only (names for `PipelineStageConfig.specify`). [[getStage]] gives
 * the resulting latency, and `z.valid` is `rnd.valid` delayed by that many
 * cycles.
 *
 * {{{
 * // two registers, placed by the module
 * new GenRandomFloat12(spec, 2)
 * // hand-picked positions
 * new GenRandomFloat12(spec, PipelineStageConfig.specify("in" -> 1, "out" -> 1))
 * // combinational
 * new GenRandomFloat12(spec, PipelineStageConfig.none)
 * }}}
 *
 * @constructor    create a chisel Module.
 * @param spec     output floating-point spec.
 * @param stage    pipeline registers of this module.
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class GenRandomFloat12(
  spec: RealSpec,
  stage: PipelineStageConfig
) extends Module {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = GenRandomFloat12.stagePoints

  val stages = stage.resolve(stagePoints)

  def getParam = { (spec, stage) }

  val nStage = stages.total

  def getStage = nStage

  val exW  = spec.exW
  val manW = spec.manW

  val io = IO(new Bundle {
    val rnd = Flipped(Valid(UInt(manW.W)))
    val z   = Valid(UInt(spec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")

  private val cutIn = new SharedPipe("in", io.rnd.valid, stages.at("in"))

  val rndIn   = cutIn(io.rnd.bits)
  val validIn = cutIn.valid

  //----------------------------------------------------------------------
  // S1 : assemble the result

  val zsgn = 0.U(1.W)
  val zex  = spec.exBias.U(exW.W)
  val zman = rndIn

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  //---------------------------- stage "out" -----------------------------
  io.z := Pipe(validIn, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.rng.GenRandomFloat12]].
 *
 * @see [[rial.util.StagePoint]]
 */
object GenRandomFloat12 {
  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.rng.GenRandomFloat12]]. */
  val stagePoints: Seq[StagePoint] = Seq(StagePoint("out", 1))
}

/** generate Float in [0, 1) based on a uniform distribution of integer.
 *
 * It interprets random integer input as a fixedpoint number rnd / (1<<rndW)
 * and convert it to floating point. the minimum number you can get is
 * 1 / (1<<rndW).
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`   the module input
 *  - `"clz"`  after the leading zero count
 *  - `"norm"` after the normalisation
 *  - `"out"`  the module output
 *
 * [[getStage]] gives the resulting latency, and `z.valid` is `rnd.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // two registers, placed by the module
 * new GenRandomFloat01(rndW, spec, 2)
 * // hand-picked positions
 * new GenRandomFloat01(rndW, spec,
 *   PipelineStageConfig.specify("clz" -> 1, "norm" -> 1))
 * // all the registers at the module output
 * new GenRandomFloat01(rndW, spec, PipelineStageConfig.atOut(2))
 * // combinational
 * new GenRandomFloat01(rndW, spec, PipelineStageConfig.none)
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
class GenRandomFloat01(
  rndW: Int, // random number width
  spec: RealSpec,
  stage: PipelineStageConfig
) extends Module {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = GenRandomFloat01.stagePoints

  val stages = stage.resolve(stagePoints)

  def getParam = { (rndW, spec, stage) }

  val nStage = stages.total

  def getStage = nStage

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  require(rndW > manW)

  val io = IO(new Bundle {
    val rnd = Flipped(Valid(UInt(rndW.W)))
    val z   = Valid(UInt(spec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")

  private val cutIn = new SharedPipe("in", io.rnd.valid, stages.at("in"))

  val rndIn   = cutIn(io.rnd.bits)
  val validIn = cutIn.valid

  //----------------------------------------------------------------------
  // S1 : leading zero count

  val zerosAtMSB = PriorityEncoder(Reverse(rndIn))
  val rndZero    = rndIn === 0.U

  //---------------------------- stage "clz" -----------------------------
  private val cutClz = new SharedPipe("clz", validIn, stages.at("clz"))

  val rndClz        = cutClz(rndIn)
  val zerosAtMSBClz = cutClz(zerosAtMSB)
  val rndZeroClz    = cutClz(rndZero)
  val validClz      = cutClz.valid

  //----------------------------------------------------------------------
  // S2 : normalisation shift

  val manW1 = (rndClz << zerosAtMSBClz)(rndW-1, rndW-1-manW)
  require(manW1.getWidth == manW+1)

  assert(!validClz || manW1(manW) === 1.U)

  //---------------------------- stage "norm" ----------------------------
  private val cutNorm = new SharedPipe("norm", validClz, stages.at("norm"))

  val manW1Norm      = cutNorm(manW1)
  val zerosAtMSBNorm = cutNorm(zerosAtMSBClz)
  val rndZeroNorm    = cutNorm(rndZeroClz)
  val validNorm      = cutNorm.valid

  //----------------------------------------------------------------------
  // S3 : exponent and assemble

  val sgn = 0.U(1.W)
  val ex  = Mux(rndZeroNorm, 0.U(exW.W), (exBias-1).U(exW.W) - zerosAtMSBNorm)
  val man = manW1Norm(manW-1, 0)

  val z = Cat(sgn, ex, man)
  require(z.getWidth == spec.W)

  //---------------------------- stage "out" -----------------------------
  io.z := Pipe(validNorm, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.rng.GenRandomFloat01]].
 *
 * @see [[rial.util.StagePoint]]
 */
object GenRandomFloat01 {
  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.rng.GenRandomFloat01]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("clz", 3), StagePoint("norm", 3), StagePoint("out", 1))
}
