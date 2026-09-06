//% @file crc.scala
//
// Fletcher's checksum
//
// part of RIAL - RIKEN Arithmetic and Logic library
// Copyright (C) Toru Niina RIKEN BDR 2025
//
package rial.ecc

import chisel3._
import chisel3.util._

import rial.util._

// Fletcher-16
//
// Ai = Ai-1 + di
// Bi = Bi-1 + Ai
// checksum = { (Bi mod 255), (Ai mod 255) }
//
// A0 = 0                  | B0 = 0
// A1 = d1                 | B1 = d1
// A2 = d1 + d2            | B2 = d1*2 + d2
// A3 = d1 + d2 + d3       | B3 = d1*3 + d2*2 + d3
// A4 = d1 + d2 + d3 + d4  | B4 = d1*4 + d2*3 + d3*2 + d4
// ...                     | ...
// An = sum({di})          | Bn = d1*n + d2*(n-1) + ... + dn
//
// max(An) = 255 * n
// max(Bn) = 255 * (n + n-1 + ... + 1) = 255 * (n * (n+1) / 2)
//
// (0 <= xi <= 255)
// x = x0 + x1 * 256     + x2 * 256^2     + ... xn * 256^n     (mod255)
//   = x0 + x1 * (1+255) + x2 * (1+255)^2 + ... xn * (1+255)^n (mod255)
//   = x0 + x1           + x2             + ... xn             (mod255)
//
// if 0 <= x <= 254, mod255(x) == x. note that mod255(0xff) == 0
//

/** Fletcher's checksum calculator
 *
 * {{{
 * class Fletcher(...) extends Module {
 *   val io = IO(new Bundle{
 *     val in  = Flipped(Valid(UInt(inputW.W)))
 *     val out = Valid(UInt(checksumW.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * The whole input word is checksummed in one shot: it is split into elements
 * of `checksumW/2` bits, `in.bits` of one cycle is one message.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`  the module input
 *  - `"sum"` after the weighted sum and the sum of the elements
 *  - `"mod"` after the modulo of the weighted sum
 *  - `"out"` the module output
 *
 * [[getStage]] gives the resulting latency, and `out.valid` is `in.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // two registers, placed by the module
 * new Fletcher(16, 64, 2)
 * // hand-picked positions
 * new Fletcher(16, 64, PipelineStageConfig.specify("sum" -> 1, "mod" -> 1))
 * // combinational
 * new Fletcher(16, 64, PipelineStageConfig.none)
 * }}}
 *
 * @constructor construct Fletcher module
 * @param checksumW checksum width. if 16, use Fletcher-16.
 * @param inputW    input width. if it is not a multiple of checksumW/2, the input will be extended with zero (zero-padding)
 * @param stage     pipeline registers of this module.
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class Fletcher(
  checksumW: Int, // if 16, use Fletcher-16.
  inputW:    Int, // if not a multiple of checksumW/2, it will be 0-padded.
  stage:     PipelineStageConfig
) extends Module {

  assert(checksumW % 2 == 0)
  val elemW = checksumW / 2
  val mod = (1 << elemW) - 1

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = Fletcher.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  def splitIntoElems(x: UInt): Seq[UInt] = {
    // 0-pad to make it a multiple of elemW
    val y = x.pad((x.getWidth + elemW - 1) / elemW * elemW)
    Seq.tabulate(y.getWidth / elemW)(i => { y((i+1)*elemW-1, i*elemW) })
  }

  def modsum(lhs: UInt, rhs: UInt): UInt = {
    assert(lhs.getWidth == elemW)
    assert(rhs.getWidth == elemW)

    val sum = lhs +& rhs
    sum(elemW-1,0) + sum(elemW)
  }

  // -------------------------------------------------------------------------

  val io = IO(new Bundle {
    val in  = Flipped(Valid(UInt(inputW.W)))
    val out = Valid(UInt(checksumW.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val xIn = cutIn(io.in.bits)

  //----------------------------------------------------------------------
  // S1 : the weighted sum Bn and, beside it, the sum of the elements An

  val d = VecInit(splitIntoElems(xIn))
  val n = d.length

  // Bn is only ever read modulo 2^elemW - 1, and there 2^elemW is 1, so
  // multiplying an element by 2^k is rotating it left by k bits. Every weight
  // therefore becomes one rotated copy of its element per set bit, which
  // leaves no multiplier and no operand wider than an element plus the depth
  // of the reduction.
  def rotl(x: UInt, k: Int): UInt = {
    val s = k % elemW
    if (s == 0) { x } else { Cat(x(elemW-1-s, 0), x(elemW-1, elemW-s)) }
  }

  val bnTerms = d.zipWithIndex.flatMap{ case(v, i) => {
    val weight = n - i
    (0 until log2Ceil(weight+1)).filter(k => ((weight >> k) & 1) == 1).map(k => rotl(v, k))
  }}
  val bnW    = elemW + log2Ceil(bnTerms.length)
  val bnComb = VecInit(bnTerms.map(_.pad(bnW))).reduceTree(_+_)

  val amod0Comb = d.reduceTree(modsum(_,_))

  //---------------------------- stage "sum" -----------------------------
  private val cutSum = new SharedPipe("sum", cutIn.valid, stages.at("sum"))

  val bn    = cutSum(bnComb)
  val amod0 = cutSum(amod0Comb)

  //----------------------------------------------------------------------
  // S2 : Bn mod 255

  val bmod0Comb = VecInit(splitIntoElems(bn)).reduceTree(modsum(_,_))

  //---------------------------- stage "mod" -----------------------------
  private val cutMod = new SharedPipe("mod", cutSum.valid, stages.at("mod"))

  val bmod0    = cutMod(bmod0Comb)
  val amod0Mod = cutMod(amod0)

  //----------------------------------------------------------------------
  // S3 : mod255(0xff) is 0, and then the two halves are concatenated

  val amodComb = Mux(amod0Mod === mod.U, 0.U(elemW.W), amod0Mod)
  val bmodComb = Mux(bmod0    === mod.U, 0.U(elemW.W), bmod0)
  val checksumComb = Cat(bmodComb, amodComb)

  assert(amodComb.getWidth == elemW)
  assert(bmodComb.getWidth == elemW)
  assert(checksumComb.getWidth == checksumW)

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutMod.valid, stages.at("out"))

  io.out.bits  := cutOut(checksumComb)
  io.out.valid := cutOut.valid
}

/** Pipeline stage points of [[rial.ecc.Fletcher]].
 *
 * @see [[rial.util.StagePoint]]
 */
object Fletcher {

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.ecc.Fletcher]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("sum", 3), StagePoint("mod", 3), StagePoint("out", 1))
}
