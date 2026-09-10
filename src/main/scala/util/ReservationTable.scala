//% @file ReservationTable.scala
//
// Reservation table of the boundaries a pipeline shares
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.util

import chisel3._
import chisel3.util._

/** Says which operations may be issued into a pipeline whose operations take
 *  different numbers of cycles to reach the boundaries they share.
 *
 * One operation at a time enters the pipeline on `io.in`, and `io.in.bits`
 * selects which one it is. An operation of selection value `s` needs
 * `latencies(b)(s)` cycles to reach shared boundary `b`, so two operations
 * issued in different cycles can still reach one of those boundaries in the
 * same cycle. `io.issueOk(s)` is high when an operation of selection value `s`
 * issued this cycle would reach every boundary in a cycle no operation already
 * in flight reaches it in; it has to be honoured, i.e. `io.in.valid` may only
 * be asserted when `io.issueOk(io.in.bits)` is high, which the table does not
 * check.
 *
 * A selection value that a `latencies(b)` does not name never reaches boundary
 * `b` and never reserves it, and a selection value no boundary names is always
 * ready to be issued. Every selection value with the same latency to a
 * boundary shares the register that reserves it, so the table costs one
 * register per cycle of the longest path to each boundary, whatever the number
 * of selection values.
 *
 * {{{
 * // three operations, two shared boundaries: the operation selected by 1
 * // reaches the first boundary after 2 cycles and the second after 5, the one
 * // selected by 2 after 3 and 5, and the one selected by 3 -- which never
 * // passes the first boundary -- after 4 cycles.
 * val table = Module(new ReservationTable(2, 4,
 *   Seq(Map(1 -> 2, 2 -> 3), Map(1 -> 5, 2 -> 5, 3 -> 4))))
 *
 * val go = want.valid && table.io.issueOk(want.bits)
 *
 * table.io.in.valid := go
 * table.io.in.bits  := want.bits
 * }}}
 *
 * @constructor create the reservation registers of every shared boundary.
 * @param selW      width of the selection signal
 * @param nSel      number of selection values `io.issueOk` answers for
 * @param latencies for each shared boundary, the cycles an operation needs to
 *                  reach it, keyed by selection value
 */
class ReservationTable(val selW: Int, val nSel: Int,
  val latencies: Seq[Map[Int, Int]]) extends Module {

  val io = IO(new Bundle {
    /** the operation issued this cycle, `bits` being its selection value. */
    val in = Flipped(Valid(UInt(selW.W)))

    /** whether an operation of each selection value may be issued now. */
    val issueOk = Output(Vec(nSel, Bool()))
  })

  /** whether the operation issued this cycle is one of `sels`. */
  private def issuedOneOf(sels: Seq[Int]): Bool = {
    io.in.valid && sels.map(s => io.in.bits === s.U).reduce(_ || _)
  }

  /** whether an operation in flight reaches boundary `b` in `k` cycles.
   *
   * `reservedAt(b)(k)` reads the register that holds it, `k` counting down as
   * the operation advances; no register is needed for the longest path, which
   * nothing in flight is still that far from.
   */
  private val reservedAt: Seq[Int => Bool] =
    latencies.zipWithIndex.map { case (lat, b) => {
      val n = if (lat.isEmpty) {
        0
      } else {
        lat.values.max
      }
      if (n == 0) {
        (k: Int) => false.B
      } else {
        val reserved = RegInit(VecInit(Seq.fill(n)(false.B)))
        reserved.suggestName(s"reserved$b")

        // one test per distinct latency: every selection value that reaches
        // the boundary in the same number of cycles reserves the same register
        val issuedWith: Map[Int, Bool] =
          lat.values.toSeq.distinct.filter(d => d > 0).map(d => {
            val sels = lat.filter { case (_, e) => e == d }.keys.toSeq.sorted
            d -> issuedOneOf(sels)
          }).toMap

        (0 until n).foreach(k => {
          val ahead = if (k + 1 < n) {
            reserved(k + 1)
          } else {
            false.B
          }
          reserved(k) := ahead || issuedWith.getOrElse(k + 1, false.B)
        })

        (k: Int) => {
          if (k < n) {
            reserved(k)
          } else {
            false.B
          }
        }
      }
    }}

  io.issueOk.zipWithIndex.foreach { case (ok, s) => {
    val knows = latencies.zipWithIndex.filter {
      case (lat, _) => lat.contains(s)
    }
    if (knows.isEmpty) {
      ok := true.B
    } else {
      ok := knows.map { case (lat, b) => !reservedAt(b)(lat(s)) }.reduce(_ && _)
    }
  }}
}
