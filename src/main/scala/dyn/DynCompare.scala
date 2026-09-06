//% @file DynCompare.scala
//
// Integer comparator whose element width is selected at run time.
//
// part of RIAL - RIKEN Arithmetic and Logic library
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.dyn

import chisel3._
import chisel3.util._

import rial.util._

// All the modes share one subtractor `a - b` built as a chain of laneW-bit
// adders, cut at the element boundaries (see [[rial.dyn.DynLane]]). From the
// chain we get, per lane,
//   - carry-out : 1 iff there was no borrow, i.e. a >= b within the element
//                 up to that lane (unsigned)
//   - the difference, whose being zero means that lane is equal
// and, taking the values at the topmost lane of an element,
//   eq = all the lanes of the element are equal
//   lt = no carry-out
//   gt = carry-out and not eq
// A signed comparison reuses the same unsigned datapath by biasing the topmost
// lane of each element, which costs one XOR gate per lane.

/** Integer comparator with a run-time selectable element width
 *
 * {{{
 * class DynCompare(...) extends Module {
 *   val io = IO(new Bundle{
 *     val in  = Flipped(Valid(new DynCompare.Input(laneW, nLanes)))
 *     val out = Valid(new DynCompare.Output(nLanes))
 *   })
 *   //...
 * }
 * }}}
 *
 * `mode` is the base-2 logarithm of the number of lanes per element, i.e. with
 * the default parameters mode=0 means int8 x 4, mode=1 means int16 x 2 and
 * mode=2 means int32. [[rial.dyn.DynLane.modeOf]] converts an element width in
 * bits into the corresponding mode value.
 *
 * The result of an element is broadcast to all the lanes it occupies, so
 * `io.out.bits.eq(i)` is always the result for the element that contains lane
 * `i` whatever the mode is. To get a per-byte mask, `Cat` of `Fill(laneW,
 * flag)` can be used. `ne`, `le` and `ge` are just `!eq`, `!gt` and `!lt`.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`  the module input
 *  - `"sub"` after the shared carry chain
 *  - `"out"` the module output
 *
 * [[getStage]] gives the resulting latency, and `out.valid` is `in.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // two registers, placed by the module
 * new DynCompare(8, 4, true, 2)
 * // hand-picked positions
 * new DynCompare(8, 4, true, PipelineStageConfig.specify("sub" -> 1, "out" -> 1))
 * // combinational
 * new DynCompare(8, 4, true, PipelineStageConfig.none)
 * }}}
 *
 * @constructor construct DynCompare module
 * @param laneW  width of the narrowest element (and of the carry chain unit)
 * @param nLanes number of lanes. must be a power of 2
 * @param signed if true, elements are compared as signed integers
 * @param stage  pipeline registers of this module
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class DynCompare(
  val laneW  : Int = 8,
  val nLanes : Int = 4,
  val signed : Boolean = true,
  val stage  : PipelineStageConfig = PipelineStageConfig.none
) extends Module {

  require(laneW   >= 1,   s"laneW must be positive, but got ${laneW}")
  require(nLanes  >= 1,   s"nLanes must be positive, but got ${nLanes}")
  require(isPow2(nLanes), s"nLanes must be a power of 2, but got ${nLanes}")

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = DynCompare.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val width = laneW * nLanes
  val modeW = DynLane.modeW(nLanes)

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new DynCompare.Input(laneW, nLanes)))
    val out = Valid(new DynCompare.Output(nLanes))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val aIn    = cutIn(io.in.bits.a)
  val bIn    = cutIn(io.in.bits.b)
  val modeIn = cutIn(io.in.bits.mode)

  //----------------------------------------------------------------------
  // S1 : the shared carry chain of a - b

  val heads = DynLane.heads(modeIn, nLanes)
  val tails = DynLane.tails(heads)

  val a = DynLane.operand(aIn, laneW, nLanes, tails, signed)
  val b = DynLane.operand(bIn, laneW, nLanes, tails, signed)

  val sub      = DynLane.subtract(a, b, heads)
  val coutComb = sub.map(s => s(laneW).asBool)

  // lane equality does not have to wait for the carry chain. At the head of an
  // element the carry-in is the forced 1, so the difference is zero exactly
  // when the two lanes are equal; at any lane above it the accumulate only
  // lets a result through when every lower lane of the element was equal,
  // which is what forces that same carry-in. The bias cancels in an equality,
  // so the raw lanes answer it.
  val eqLaneComb = DynLane.split(aIn, laneW, nLanes)
    .zip(DynLane.split(bIn, laneW, nLanes))
    .map{ case(x, y) => x === y }
  val zaccComb = DynLane.accumulate(eqLaneComb, heads)

  //---------------------------- stage "sub" -----------------------------
  private val cutSub = new SharedPipe("sub", cutIn.valid, stages.at("sub"))

  val cout    = coutComb.map(c => cutSub(c))
  val zacc    = zaccComb.map(z => cutSub(z))
  val modeSub = cutSub(modeIn)

  //----------------------------------------------------------------------
  // S2 : result of the element that ends at lane i. only valid where tails(i).

  val eqAt = zacc
  val ltAt = cout.map(c => !c)
  val gtAt = cout.zip(zacc).map{ case(c, z) => c && !z }

  // give every lane the result of the element that contains it.
  val eqComb = DynLane.broadcast(eqAt, modeSub, nLanes)
  val ltComb = DynLane.broadcast(ltAt, modeSub, nLanes)
  val gtComb = DynLane.broadcast(gtAt, modeSub, nLanes)

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutSub.valid, stages.at("out"))

  for(i <- 0 until nLanes) {
    io.out.bits.eq(i) := cutOut(eqComb(i))
    io.out.bits.lt(i) := cutOut(ltComb(i))
    io.out.bits.gt(i) := cutOut(gtComb(i))
  }
  io.out.valid := cutOut.valid
}

/** Input and output bundles and pipeline stage points of
 *  [[rial.dyn.DynCompare]].
 *
 * @see [[rial.util.StagePoint]]
 */
object DynCompare {

  /** Input of [[rial.dyn.DynCompare]]: the two operands and the element width.
   *
   * @param laneW  width of the narrowest element
   * @param nLanes number of lanes
   */
  class Input(laneW: Int, nLanes: Int) extends Bundle {
    val a    = UInt((laneW * nLanes).W)
    val b    = UInt((laneW * nLanes).W)
    val mode = UInt(DynLane.modeW(nLanes).W)
  }

  /** Output of [[rial.dyn.DynCompare]]: the three flags per lane.
   *
   * @param nLanes number of lanes
   */
  class Output(nLanes: Int) extends Bundle {
    val eq = Vec(nLanes, Bool())
    val lt = Vec(nLanes, Bool())
    val gt = Vec(nLanes, Bool())
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.dyn.DynCompare]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("sub", 3), StagePoint("out", 1))
}
