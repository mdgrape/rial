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

// generate Float in [1, 2).
class GenRandomFloat12(
  rndW: Int, // random number width
  spec: RealSpec,
  stage: PipelineStageConfig
  ) extends Module {

  def getParam = { (rndW, spec, stage) }

  val nStage = stage.total

  val exW  = spec.exW
  val manW = spec.manW

  // required number of input integers to generate float
  val nRndInt = (manW / rndW) + 1

  val io = IO(new Bundle {
    val rnds = Input(Vec(nRndInt, UInt(rndW.W)))
    val z    = Output(UInt(spec.W.W))
  })

  val zsgn = 0.U(1.W)
  val zex  = spec.exBias.U(exW.W)
  val zman = (io.rnds.asUInt)(manW-1, 0)

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := ShiftRegister(z, nStage)
}

// generate Float in [0, 1) based on a uniform distribution of integer.
//
// Since bits in a random number vector is a multiple of the width of a random
// integer, the precision of resulting FP increases a bit compared to the case
// when we simply do `GenRandomFloat12 - 1.0`.
class GenRandomFloat01(
  rndW: Int, // random number width
  spec: RealSpec,
  stage: PipelineStageConfig
  ) extends Module {

  def getParam = { (rndW, spec, stage) }

  val nStage = stage.total

  val exW  = spec.exW
  val manW = spec.manW
  val exBias = spec.exBias

  // minimum required number of input integers to generate float
  val nRndInt = (manW / rndW) + 1

  val io = IO(new Bundle {
    val rnds = Input(Vec(nRndInt, UInt(rndW.W)))
    val z    = Output(UInt(spec.W.W))
  })

  // consider rndBits is a fraction in [0, 1.0). 1.0 corresponds to 1 << rndBits.W.
  val rndBits    = io.rnds.asUInt
  val zerosAtLSB = PriorityEncoder(rndBits)

  // normalize rndBits using zerrosAtLSB. First we normalize LSB as the
  // leading 1 in FP and then reverse it.
  val zsgn = 0.U(1.W)
  val zex        = Mux(rndBits === 0.U, 0.U(exW.W), (exBias-1).U(exW.W) - zerosAtLSB)
  val zmanW1     = Reverse((rndBits >> zerosAtLSB)(manW, 0))

  assert(exBias - (rndBits.getWidth+1) >= 0)
  assert(zmanW1(manW) === 1.U)

  val z = Cat(zsgn, zex, zmanW1.tail(1))
  assert(z.getWidth == spec.W)

  io.z := ShiftRegister(z, nStage)
}

// generate Float in [0, 1).
// This generates all the possible normalized numbers in [0, 1).
class GenRandomFloat01Full(
  rndW: Int, // random number width
  spec: RealSpec,
  stage: PipelineStageConfig
  ) extends Module {

  def getParam = { (rndW, spec, stage) }

  val nStage = stage.total

  val exW  = spec.exW
  val manW = spec.manW
  val exBias = spec.exBias

  // If the exponent differs by 1, the region that has smaller exponent has
  // 2 times less numbers than the larger exponent region. That means that,
  // there are the same number of FP numbers between a region that corresponds
  // to exponent = n and exponent = m where m < n. So the probability of getting
  // exponent = n is 1/2.
  // So, to determine the exponent, we do "coin toss" until we get "head".
  // Each bit of random numbers should have no correlation, so we can use those
  // bits as a result of coin toss. This means that, the exponent can be
  // calculated as the 0s at one end (or 1s, whatever).

  val nBitsForEx = 0 - (0 - exBias) // including subnormal
  val nRndInt    = ((nBitsForEx + manW) / rndW) + 1

  val io = IO(new Bundle {
    val rnds = Input(Vec(nRndInt, UInt(rndW.W)))
    val z    = Output(UInt(spec.W.W))
  })

  val rndBits   = io.rnds.asUInt

  val rndForEx  = rndBits(nBitsForEx-1, 0)
  val rndForMan = rndBits(manW-1 + nBitsForEx, nBitsForEx)

  val zsgn = 0.U(1.W)
  val zex  = Mux(rndForEx === 0.U, 0.U(exW.W),
                 (exBias - 1).U(exW.W) - PriorityEncoder(rndForEx))
  val zman0 = rndForMan
  val zman = if(spec.disableSubnormal) {
    Mux(zex === 0.U, 0.U(manW.W), zman0)
  } else {
    zman0
  }

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := ShiftRegister(z, nStage)
}

// generate Float in [0, 1).
// This generates all the possible normalized numbers in [0, 1).
class GenRandomFloat01FullPipe(
  rndW: Int, // random number width
  spec: RealSpec,
  ) extends Module {

  def getParam = { (rndW, spec) }

  val exW  = spec.exW
  val manW = spec.manW
  val exBias = spec.exBias

  // If the exponent differs by 1, the region that has smaller exponent has
  // 2 times less numbers than the larger exponent region. That means that,
  // there are the same number of FP numbers between a region that corresponds
  // to exponent = n and exponent = m where m < n. So the probability of getting
  // exponent = n is 1/2.
  // So, to determine the exponent, we do "coin toss" until we get "head".
  // Each bit of random numbers should have no correlation, so we can use those
  // bits as a result of coin toss. This means that, the exponent can be
  // calculated as the 0s at one end (or 1s, whatever).

  val nRndsForEx  = (exBias / rndW) + (if(exBias % rndW == 0) {0} else {1})
  val nRndsForMan = (manW   / rndW) + (if(manW   % rndW == 0) {0} else {1})

  val io = IO(new Bundle {
    val rnd   = Input(UInt(rndW.W))
    val valid = Output(Bool())
    val z     = Output(UInt(spec.W.W))
  })

  val doneEx   = RegInit(false.B)
  val countMan = RegInit(nRndsForMan.U)

  val zsgn  = 0.U(1.W)
  val zex   = RegInit((exBias - 1).U(exW.W))
  val zman  = RegInit(0.U(manW.W))

  // default output (overwritten later)
  io.valid := false.B
  io.z     := 0.U

  when(doneEx) {
    val countManNext = countMan - 1.U
    val zmanNext     = (zman << rndW.U)(manW-1, 0) | (io.rnd)(manW.min(rndW)-1, 0)
    when(countManNext === 0.U) {
      // output and reset
      val z = Cat(zsgn, zex, zmanNext)
      assert(z.getWidth == spec.W)
      io.z     := z
      io.valid := true.B

      zex      := (exBias-1).U(exW.W)
      zman     := 0.U
      doneEx   := false.B
      countMan := nRndsForMan.U
    } otherwise {
      zman     := zmanNext
      countMan := countManNext
    }
  } otherwise {
    val rndHas1 = io.rnd =/= 0.U
    val nZeros  = Mux( !rndHas1, rndW.U, PriorityEncoder(io.rnd))
    val zexNext = zex -& nZeros
    val zexNeg  = zexNext.head(1) === 1.U
    zex    := Mux(zexNeg, 0.U(exW.W), zexNext)
    doneEx := rndHas1 || zexNeg
  }
}
