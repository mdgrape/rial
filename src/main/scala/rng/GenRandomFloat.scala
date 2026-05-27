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
 * @constructor    create a chisel Module.
 * @param spec     output floating-point spec.
 * @param stage    pipeline stages.
 */
class GenRandomFloat12(
  spec: RealSpec,
  stage: PipelineStageConfig
) extends Module {

  def getParam = { (spec, stage) }

  val nStage = stage.total

  val exW  = spec.exW
  val manW = spec.manW

  val io = IO(new Bundle {
    val rnd = Flipped(Valid(UInt(manW.W)))
    val z   = Valid(UInt(spec.W.W))
  })

  val zsgn = 0.U(1.W)
  val zex  = spec.exBias.U(exW.W)
  val zman = io.rnd.bits

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := Pipe(io.rnd.valid, z, nStage)
}

/** generate Float in [0, 1) based on a uniform distribution of integer.
 *
 * It interprets random integer input as a fixedpoint number rnd / (1<<rndW)
 * and convert it to floating point. the minimum number you can get is
 * 1 / (1<<rndW).
 *
 * @constructor    create a chisel Module.
 * @param rndW     input random number bit width.
 * @param spec     output floating-point spec.
 * @param stage    pipeline stages.
 */
class GenRandomFloat01(
  rndW: Int, // random number width
  spec: RealSpec,
  stage: PipelineStageConfig
) extends Module {

  def getParam = { (rndW, spec, stage) }

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  require(rndW > manW)

  val io = IO(new Bundle {
    val rnd = Flipped(Valid(UInt(rndW.W)))
    val z   = Valid(UInt(spec.W.W))
  })

  val rnd = io.rnd.bits

  val zerosAtMSB = PriorityEncoder(Reverse(rnd))

  val manW1 = (rnd << zerosAtMSB)(rndW-1, rndW-1-manW)
  require(manW1.getWidth == manW+1)

  assert(!io.rnd.valid || manW1(manW) === 1.U)

  val sgn = 0.U(1.W)
  val ex  = Mux(rnd === 0.U, 0.U(exW.W), (exBias-1).U(exW.W) - zerosAtMSB)
  val man = manW1(manW-1, 0)

  val z = Cat(sgn, ex, man)
  require(z.getWidth == spec.W)

  io.z := Pipe(io.rnd.valid, z, nStage)
}
