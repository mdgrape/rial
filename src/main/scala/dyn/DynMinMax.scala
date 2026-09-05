//% @file DynMinMax.scala
//
// Integer min/max whose element width is selected at run time.
//
// part of RIAL - RIKEN Arithmetic and Logic library
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.dyn

import chisel3._
import chisel3.util._

import rial.util._

// The selector comes from the same carry chain as [[rial.dyn.DynCompare]] uses:
// the carry-out at the topmost lane of an element is 1 iff a >= b there. Note
// that the equality chain of DynCompare is not needed here, because when the
// two elements are equal both operands give the same answer; only the carry
// chain is built.
//
// The selector of an element is broadcast to all the lanes it occupies, so the
// output is assembled lane by lane. In the int8x4 mode the four lanes pick
// their operand independently, i.e. the output is a re-mixture of the two
// inputs, while in the int32 mode all the lanes select the same operand.
//
//   a = | a3 | a2 | a1 | a0 |
//   b = | b3 | b2 | b1 | b0 |
//                              mode=0 (int8 x 4) : | a3 | b2 | b1 | a0 | ...
//                              mode=1 (int16 x 2): |   b1:b0 |   a1:a0 | ...
//                              mode=2 (int32)    : |     a3:a0         | ...

/** Integer min/max with a run-time selectable element width
 *
 * {{{
 * class DynMinMax(...) extends Module {
 *   val io = IO(new Bundle{
 *     val in  = Flipped(Valid(new DynMinMax.Input(laneW, nLanes)))
 *     val out = Valid(new DynMinMax.Output(laneW, nLanes))
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
 * The selection is made per element, so in the int8x4 mode each byte of the
 * output comes from whichever input has the smaller (or larger) byte there, and
 * `io.out.bits.min` is generally a mixture of the bytes of the two inputs.
 *
 * `min` and `max` share one comparison, hence a module that outputs both costs
 * almost nothing more than one that outputs only either of them; it is directly
 * usable as the compare-exchange element of a sorting network. Unused output is
 * eliminated by the synthesis tool.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`  the module input
 *  - `"cmp"` after the comparison
 *  - `"out"` the module output
 *
 * [[getStage]] gives the resulting latency, and `out.valid` is `in.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // two registers, placed by the module
 * new DynMinMax(8, 4, true, 2)
 * // hand-picked positions
 * new DynMinMax(8, 4, true, PipelineStageConfig.specify("cmp" -> 1, "out" -> 1))
 * // combinational
 * new DynMinMax(8, 4, true, PipelineStageConfig.none)
 * }}}
 *
 * @constructor construct DynMinMax module
 * @param laneW  width of the narrowest element (and of the carry chain unit)
 * @param nLanes number of lanes. must be a power of 2
 * @param signed if true, elements are compared as signed integers
 * @param stage  pipeline registers of this module
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class DynMinMax(
  val laneW  : Int = 8,
  val nLanes : Int = 4,
  val signed : Boolean = true,
  val stage  : PipelineStageConfig = PipelineStageConfig.none
) extends Module {

  require(laneW   >= 1,   s"laneW must be positive, but got ${laneW}")
  require(nLanes  >= 1,   s"nLanes must be positive, but got ${nLanes}")
  require(isPow2(nLanes), s"nLanes must be a power of 2, but got ${nLanes}")

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = DynMinMax.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val width = laneW * nLanes
  val modeW = DynLane.modeW(nLanes)

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new DynMinMax.Input(laneW, nLanes)))
    val out = Valid(new DynMinMax.Output(laneW, nLanes))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val aIn    = cutIn(io.in.bits.a)
  val bIn    = cutIn(io.in.bits.b)
  val modeIn = cutIn(io.in.bits.mode)

  //----------------------------------------------------------------------
  // S1 : the carry chain of a - b. the carry-out at the topmost lane of an
  // element is 1 iff a >= b there. when a == b either operand gives the same
  // answer, so the equality chain that DynCompare needs is not built here.

  val heads = DynLane.heads(modeIn, nLanes)
  val tails = DynLane.tails(heads)

  val a = DynLane.operand(aIn, laneW, nLanes, tails, signed)
  val b = DynLane.operand(bIn, laneW, nLanes, tails, signed)

  val geAt   = DynLane.subtract(a, b, heads).map(s => s(laneW).asBool)
  val geComb = DynLane.broadcast(geAt, modeIn, nLanes)

  //---------------------------- stage "cmp" -----------------------------
  private val cutCmp = new SharedPipe("cmp", cutIn.valid, stages.at("cmp"))

  val ge   = geComb.map(g => cutCmp(g))
  val aCmp = cutCmp(aIn)
  val bCmp = cutCmp(bIn)

  //----------------------------------------------------------------------
  // S2 : select the original (unbiased) lanes

  val al = DynLane.split(aCmp, laneW, nLanes)
  val bl = DynLane.split(bCmp, laneW, nLanes)

  val minComb = DynLane.concat(Seq.tabulate(nLanes)(i => Mux(ge(i), bl(i), al(i))))
  val maxComb = DynLane.concat(Seq.tabulate(nLanes)(i => Mux(ge(i), al(i), bl(i))))

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutCmp.valid, stages.at("out"))

  io.out.bits.min := cutOut(minComb)
  io.out.bits.max := cutOut(maxComb)
  io.out.valid    := cutOut.valid
}

/** Input and output bundles and pipeline stage points of
 *  [[rial.dyn.DynMinMax]].
 *
 * @see [[rial.util.StagePoint]]
 */
object DynMinMax {

  /** Input of [[rial.dyn.DynMinMax]]: the two operands and the element width.
   *
   * @param laneW  width of the narrowest element
   * @param nLanes number of lanes
   */
  class Input(laneW: Int, nLanes: Int) extends Bundle {
    val a    = UInt((laneW * nLanes).W)
    val b    = UInt((laneW * nLanes).W)
    val mode = UInt(DynLane.modeW(nLanes).W)
  }

  /** Output of [[rial.dyn.DynMinMax]]: the lane-wise minimum and maximum.
   *
   * @param laneW  width of the narrowest element
   * @param nLanes number of lanes
   */
  class Output(laneW: Int, nLanes: Int) extends Bundle {
    val min = UInt((laneW * nLanes).W)
    val max = UInt((laneW * nLanes).W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.dyn.DynMinMax]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("cmp", 3), StagePoint("out", 1))
}
