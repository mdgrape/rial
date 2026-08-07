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
 *     val in  = Input(new Bundle {
 *       val a    = UInt(width.W)
 *       val b    = UInt(width.W)
 *       val mode = UInt(modeW.W)
 *     })
 *     val out = Output(new Bundle {
 *       val min = UInt(width.W)
 *       val max = UInt(width.W)
 *     })
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
 * `io.out.min` is generally a mixture of the bytes of the two inputs.
 *
 * `min` and `max` share one comparison, hence a module that outputs both costs
 * almost nothing more than one that outputs only either of them; it is directly
 * usable as the compare-exchange element of a sorting network. Unused output is
 * eliminated by the synthesis tool.
 *
 * @constructor construct DynMinMax module
 * @param laneW   width of the narrowest element (and of the carry chain unit)
 * @param nLanes  number of lanes. must be a power of 2
 * @param signed  if true, elements are compared as signed integers
 * @param latency latency of this module
 */
class DynMinMax(
  val laneW   : Int = 8,
  val nLanes  : Int = 4,
  val signed  : Boolean = true,
  val latency : Int = 0
) extends Module {

  require(laneW   >= 1,   s"laneW must be positive, but got ${laneW}")
  require(nLanes  >= 1,   s"nLanes must be positive, but got ${nLanes}")
  require(isPow2(nLanes), s"nLanes must be a power of 2, but got ${nLanes}")
  require(latency >= 0,   s"latency must not be negative, but got ${latency}")

  val width = laneW * nLanes
  val modeW = DynLane.modeW(nLanes)

  val io = IO(new Bundle {
    val in  = Input(new Bundle {
      val a    = UInt(width.W)
      val b    = UInt(width.W)
      val mode = UInt(modeW.W)
    })
    val out = Output(new Bundle {
      val min = UInt(width.W)
      val max = UInt(width.W)
    })
  })

  val mode  = io.in.mode
  val heads = DynLane.heads(mode, nLanes)
  val tails = DynLane.tails(heads)

  val a = DynLane.operand(io.in.a, laneW, nLanes, tails, signed)
  val b = DynLane.operand(io.in.b, laneW, nLanes, tails, signed)

  // the carry chain of a - b. the carry-out at the topmost lane of an element
  // is 1 iff a >= b there. when a == b either operand gives the same answer, so
  // the equality chain that DynCompare needs is not built here.
  val geAt = DynLane.subtract(a, b, heads).map(s => s(laneW).asBool)
  val ge   = DynLane.broadcast(geAt, mode, nLanes)

  // select the original (unbiased) lanes
  val al = DynLane.split(io.in.a, laneW, nLanes)
  val bl = DynLane.split(io.in.b, laneW, nLanes)

  val min = DynLane.concat(Seq.tabulate(nLanes)(i => Mux(ge(i), bl(i), al(i))))
  val max = DynLane.concat(Seq.tabulate(nLanes)(i => Mux(ge(i), al(i), bl(i))))

  io.out.min := ShiftRegister(min, latency)
  io.out.max := ShiftRegister(max, latency)
}
