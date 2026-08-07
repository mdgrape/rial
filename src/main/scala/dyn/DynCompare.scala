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
 *     val in  = Input(new Bundle {
 *       val a    = UInt(width.W)
 *       val b    = UInt(width.W)
 *       val mode = UInt(modeW.W)
 *     })
 *     val out = Output(new Bundle {
 *       val eq = Vec(nLanes, Bool())
 *       val lt = Vec(nLanes, Bool())
 *       val gt = Vec(nLanes, Bool())
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
 * The result of an element is broadcast to all the lanes it occupies, so
 * `io.out.eq(i)` is always the result for the element that contains lane `i`
 * whatever the mode is. To get a per-byte mask, `Cat` of `Fill(laneW, flag)`
 * can be used. `ne`, `le` and `ge` are just `!eq`, `!gt` and `!lt`.
 *
 * @constructor construct DynCompare module
 * @param laneW   width of the narrowest element (and of the carry chain unit)
 * @param nLanes  number of lanes. must be a power of 2
 * @param signed  if true, elements are compared as signed integers
 * @param latency latency of this module
 */
class DynCompare(
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
      val eq = Vec(nLanes, Bool())
      val lt = Vec(nLanes, Bool())
      val gt = Vec(nLanes, Bool())
    })
  })

  val mode  = io.in.mode
  val heads = DynLane.heads(mode, nLanes)
  val tails = DynLane.tails(heads)

  val a = DynLane.operand(io.in.a, laneW, nLanes, tails, signed)
  val b = DynLane.operand(io.in.b, laneW, nLanes, tails, signed)

  // the shared carry chain of a - b
  val sub  = DynLane.subtract(a, b, heads)
  val cout = sub.map(s => s(laneW).asBool)
  val zacc = DynLane.accumulate(sub.map(s => s(laneW-1, 0) === 0.U), heads)

  // result of the element that ends at lane i. only valid where tails(i).
  val eqAt = zacc
  val ltAt = cout.map(c => !c)
  val gtAt = cout.zip(zacc).map{ case(c, z) => c && !z }

  // give every lane the result of the element that contains it.
  val eq = DynLane.broadcast(eqAt, mode, nLanes)
  val lt = DynLane.broadcast(ltAt, mode, nLanes)
  val gt = DynLane.broadcast(gtAt, mode, nLanes)

  for(i <- 0 until nLanes) {
    io.out.eq(i) := ShiftRegister(eq(i), latency)
    io.out.lt(i) := ShiftRegister(lt(i), latency)
    io.out.gt(i) := ShiftRegister(gt(i), latency)
  }
}
