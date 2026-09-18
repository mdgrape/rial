//% @file DynAdd.scala
//
// Integer adder whose element width is selected at run time.
//
// part of RIAL - RIKEN Arithmetic and Logic library
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.dyn

import chisel3._
import chisel3.util._

import rial.util._

// All the modes share one adder `a + b` built as a chain of laneW-bit adders,
// cut at the element boundaries (see [[rial.dyn.DynLane]]): at a lane that
// starts a new element the carry-in is forced to 0, so the lower elements do
// not reach into the upper ones. What is left of an element is its sum modulo
// its own width, which is all the output keeps.
//
// The whole datapath is that one chain; the mode only picks where it is cut.

/** Integer adder with a run-time selectable element width
 *
 * {{{
 * class DynAdd(...) extends Module {
 *   val io = IO(new Bundle{
 *     val in  = Flipped(Valid(new DynAdd.Input(laneW, nLanes)))
 *     val out = Valid(new DynAdd.Output(laneW, nLanes))
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
 * `sum` is the element-wise sum, each element wrapped around its own width, and
 * it occupies the same bits of the output as its operands do of the inputs. The
 * bits of a wrapped sum are the same whether the elements are read as signed or
 * as unsigned integers, so the module has no `signed` parameter; a carry-out is
 * not reported for the same reason. To detect an overflow, widen the operands
 * (`DynAdd(2*laneW, nLanes)` over sign or zero extended elements).
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`  the module input
 *  - `"out"` the module output
 *
 * [[getStage]] gives the resulting latency, and `out.valid` is `in.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // two registers, placed by the module
 * new DynAdd(8, 4, 2)
 * // hand-picked positions
 * new DynAdd(8, 4, PipelineStageConfig.specify("in" -> 1, "out" -> 1))
 * // combinational
 * new DynAdd(8, 4, PipelineStageConfig.none)
 * }}}
 *
 * @constructor construct DynAdd module
 * @param laneW  width of the narrowest element (and of the carry chain unit)
 * @param nLanes number of lanes. must be a power of 2
 * @param stage  pipeline registers of this module
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class DynAdd(
  val laneW  : Int = 8,
  val nLanes : Int = 4,
  val stage  : PipelineStageConfig = PipelineStageConfig.none
) extends Module {

  require(laneW   >= 1,   s"laneW must be positive, but got ${laneW}")
  require(nLanes  >= 1,   s"nLanes must be positive, but got ${nLanes}")
  require(isPow2(nLanes), s"nLanes must be a power of 2, but got ${nLanes}")

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = DynAdd.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val width = laneW * nLanes
  val modeW = DynLane.modeW(nLanes)

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new DynAdd.Input(laneW, nLanes)))
    val out = Valid(new DynAdd.Output(laneW, nLanes))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val aIn    = cutIn(io.in.bits.a)
  val bIn    = cutIn(io.in.bits.b)
  val modeIn = cutIn(io.in.bits.mode)

  //----------------------------------------------------------------------
  // S1 : the shared carry chain of a + b. the carry-out each lane returns is
  // the bit above that lane, which belongs to the next element or to nothing,
  // so only the sum bits are kept.

  val heads = DynLane.heads(modeIn, nLanes)

  val a = DynLane.split(aIn, laneW, nLanes)
  val b = DynLane.split(bIn, laneW, nLanes)

  val sumComb = DynLane.concat(
    DynLane.add(a, b, heads).map(s => s(laneW-1, 0)))

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutIn.valid, stages.at("out"))

  io.out.bits.sum := cutOut(sumComb)
  io.out.valid    := cutOut.valid
}

/** Input and output bundles and pipeline stage points of [[rial.dyn.DynAdd]].
 *
 * @see [[rial.util.StagePoint]]
 */
object DynAdd {

  /** Input of [[rial.dyn.DynAdd]]: the two operands and the element width.
   *
   * @param laneW  width of the narrowest element
   * @param nLanes number of lanes
   */
  class Input(laneW: Int, nLanes: Int) extends Bundle {
    val a    = UInt((laneW * nLanes).W)
    val b    = UInt((laneW * nLanes).W)
    val mode = UInt(DynLane.modeW(nLanes).W)
  }

  /** Output of [[rial.dyn.DynAdd]]: the element-wise sum.
   *
   * @param laneW  width of the narrowest element
   * @param nLanes number of lanes
   */
  class Output(laneW: Int, nLanes: Int) extends Bundle {
    val sum = UInt((laneW * nLanes).W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.dyn.DynAdd]]. */
  val stagePoints: Seq[StagePoint] = Seq(StagePoint("out", 1))
}
