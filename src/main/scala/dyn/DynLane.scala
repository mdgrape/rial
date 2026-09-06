//% @file DynLane.scala
//
// Building blocks shared by the run-time selectable element width operators.
//
// part of RIAL - RIKEN Arithmetic and Logic library
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.dyn

import chisel3._
import chisel3.util._

import scala.collection.mutable.ArrayBuffer

/** Building blocks for operators with a run-time selectable element width
 *
 * An input of `laneW * nLanes` bits is split into `nLanes` lanes of `laneW`
 * bits, and an element consists of `2^mode` consecutive lanes. With the default
 * configuration of the modules in this package (laneW=8, nLanes=4) the modes
 * are int8 x 4, int16 x 2 and int32 x 1.
 *
 * {{{
 *   lane        3        2        1        0
 *            +--------+--------+--------+--------+
 *            |        |        |        |        |
 *            +--------+--------+--------+--------+
 *   mode=0   |   e3   |   e2   |   e1   |   e0   |   8bit x 4
 *   mode=1   |       e1        |       e0        |  16bit x 2
 *   mode=2   |                e0                 |  32bit x 1
 * }}}
 *
 * The mode dependent part of an operator is only where the carry chain is cut:
 * at a lane that starts a new element the carry-in is forced instead of being
 * taken from the lane below, which makes the lower lanes invisible to the upper
 * ones. That is exactly what a narrower element width means, so the datapath
 * itself is shared by all the modes.
 *
 * These functions are combinational and take/return plain `Seq`s so that an
 * operator can pick only the parts it needs; the pieces it does not use are
 * simply not built.
 */
object DynLane {

  /** mode value of the widest element */
  def modeMax(nLanes: Int): Int = log2Ceil(nLanes)

  /** width of the mode signal */
  def modeW(nLanes: Int): Int = math.max(1, log2Ceil(modeMax(nLanes)+1))

  /** mode value that makes an operator interpret its input as `elemW` bit integers
   *
   * @param elemW element width in bits
   * @param laneW width of the narrowest element of the operator
   */
  def modeOf(elemW: Int, laneW: Int = 8): Int = {
    require(elemW % laneW == 0, s"elemW(${elemW}) must be a multiple of laneW(${laneW})")
    require(isPow2(elemW / laneW), s"elemW(${elemW}) / laneW(${laneW}) must be a power of 2")
    log2Ceil(elemW / laneW)
  }

  /** `heads(i)` is true iff lane i starts a new element
   *
   * lane i starts a new element iff 2^mode divides i, i.e. iff mode <= ctz(i),
   * which is a comparison against a constant.
   */
  def heads(mode: UInt, nLanes: Int): Seq[Bool] = Seq.tabulate(nLanes)(i =>
    if (i == 0) { true.B } else { mode <= Integer.numberOfTrailingZeros(i).U }
  )

  /** `tails(i)` is true iff lane i ends an element, i.e. iff lane i+1 starts one */
  def tails(heads: Seq[Bool]): Seq[Bool] = Seq.tabulate(heads.length)(i =>
    if (i == heads.length-1) { true.B } else { heads(i+1) }
  )

  /** split x into lanes. `[0]` is the LSB side */
  def split(x: UInt, laneW: Int, nLanes: Int): Seq[UInt] = {
    require(x.getWidth == laneW * nLanes)
    Seq.tabulate(nLanes)(i => x(laneW*(i+1)-1, laneW*i))
  }

  /** concatenate lanes into one value. `[0]` is the LSB side */
  def concat(lanes: Seq[UInt]): UInt = Cat(lanes.reverse)

  /** bias the topmost lane of each element by 2^(laneW-1)
   *
   * It maps the signed order onto the unsigned one, so that the rest of the
   * datapath can be shared between the signed and the unsigned case. Since the
   * same bias is applied to both operands, it does not affect equality.
   */
  def bias(lanes: Seq[UInt], tails: Seq[Bool]): Seq[UInt] = {
    require(lanes.length == tails.length)

    lanes.zip(tails).map{ case(x, t) => {
      val w = x.getWidth
      if (w == 1) { x ^ t.asUInt } else { Cat(x(w-1) ^ t, x(w-2, 0)) }
    }}
  }

  /** split x into lanes, biased if the elements are signed */
  def operand(x: UInt, laneW: Int, nLanes: Int, tails: Seq[Bool], signed: Boolean): Seq[UInt] = {
    val lanes = split(x, laneW, nLanes)
    if (signed) { bias(lanes, tails) } else { lanes }
  }

  /** lane-wise `a - b` as one carry chain over guard interleaved operands
   *
   * The carry-in of a lane that starts a new element is forced to 1 (the `+1`
   * of the two's complement negation) and otherwise comes from the lane below,
   * so the chain is cut at each element boundary. Instead of gating the carry
   * between separate per-lane adders, a guard bit is inserted below every lane
   * but the lowest, holding that lane's head flag on the `a` side against the
   * 1 that the inversion puts on the `b` side. The guard's carry-out is then
   * `maj(head, 1, carry from below) = head | carry`, which is exactly the
   * carry-in the element width asks for, computed inside the chain. The sum
   * bits at the lane positions are unchanged by the guards.
   *
   * @return per-lane result of laneW+1 bits. bit laneW is the carry-out, which
   *         is 1 iff there was no borrow, i.e. iff a >= b (unsigned) within the
   *         element up to that lane. bits laneW-1 to 0 are the difference.
   */
  def subtract(a: Seq[UInt], b: Seq[UInt], heads: Seq[Bool]): Seq[UInt] = {
    require(a.nonEmpty)
    require(a.length == b.length)
    require(a.length == heads.length)

    val laneW = a.head.getWidth
    require(a.forall(_.getWidth == laneW) && b.forall(_.getWidth == laneW))

    val n  = a.length
    val ag = Cat(((n-1) to 1 by -1).flatMap(i => Seq(a(i), heads(i).asUInt)) :+ a(0))
    val bg = Cat(((n-1) to 1 by -1).flatMap(i => Seq(b(i), 0.U(1.W))) :+ b(0))
    val sum = (ag +& (~bg).asUInt) + 1.U // n*(laneW+1) bits, never overflows

    // the guard above lane i holds head(i+1) ^ 1 ^ cout(i), so the carry-out
    // of every lane but the topmost is read back off the chain
    Seq.tabulate(n)(i => {
      val lo   = i * (laneW + 1)
      val cout = if (i == n-1) {
        sum(n * (laneW + 1) - 1)
      } else {
        ~(sum(lo + laneW) ^ heads(i+1).asUInt)
      }
      Cat(cout, sum(lo + laneW - 1, lo))
    })
  }

  /** AND-accumulate x from the head of each element up to each lane */
  def accumulate(x: Seq[Bool], heads: Seq[Bool]): Seq[Bool] = {
    require(x.length == heads.length)

    val acc = new ArrayBuffer[Bool]
    for(i <- 0 until x.length) {
      acc += (if (i == 0) { x(i) } else { x(i) && (heads(i) || acc(i-1)) })
    }
    acc.toSeq
  }

  /** broadcast the value at the tail lane of each element to all its lanes
   *
   * The element that contains lane i ends at lane `i | (2^mode - 1)`, so each
   * lane selects among at most `modeMax+1` candidates, many of which coincide.
   */
  def broadcast[T <: Data](v: Seq[T], mode: UInt, nLanes: Int): Seq[T] = Seq.tabulate(nLanes)(i => {
    val idx = (0 to modeMax(nLanes)).map(m => i | ((1 << m) - 1))
    if (idx.distinct.length == 1) {
      v(idx.head)
    } else {
      MuxLookup(mode, v(idx.last))(idx.zipWithIndex.map{ case(j, m) => (m.U, v(j)) })
    }
  })
}
