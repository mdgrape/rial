package rial.tests

import chisel3._
import chisel3.simulator.scalatest.ChiselSim

import org.scalatest.funspec.AnyFunSpec

import rial.util.ReservationTable

import scala.collection.mutable.ArrayBuffer
import scala.util.Random

/** [[rial.util.ReservationTable]]: the cycles in which it lets each selection
 *  value be issued, against a software model of the boundaries in flight.
 */
class ReservationTableTest extends AnyFunSpec with ChiselSim {

  val r = new Random(20260910)

  /** how many cycles each table under test is driven for. */
  val cycles = 500

  /** how often the issuer offers an operation, when it is allowed to. */
  val density = 0.7

  private def runtest(name: String, selW: Int, nSel: Int,
    latencies: Seq[Map[Int, Int]]) = {

    it(name) {
      simulate(new ReservationTable(selW, nSel, latencies)) { c =>

        // the operations still on their way to a boundary: the selection
        // value and the cycle it was issued in
        val inFlight = ArrayBuffer.empty[(Int, Int)]

        /** whether an operation of `s` may be issued in cycle `now`: none in
         *  flight reaches a boundary both of them pass in the same cycle.
         */
        def free(now: Int, s: Int): Boolean = {
          !inFlight.exists { case (g, at) => {
            latencies.exists(lat => {
              lat.contains(s) && lat.contains(g) &&
                at + lat(g) == now + lat(s)
            })
          }}
        }

        /** the last cycle in which an operation of `s` still reserves a
         *  boundary.
         */
        def lastAt(s: Int): Int = {
          latencies.filter(lat => lat.contains(s)).map(lat => lat(s))
            .foldLeft(0)((a, b) => a.max(b))
        }

        var issued = 0

        for (now <- 0 until cycles) {
          inFlight --= inFlight.filter { case (g, at) => at + lastAt(g) < now }

          // what the table answers now is decided by the operations issued in
          // earlier cycles, so ask before issuing one more
          val expectOk = (0 until nSel).map(s => free(now, s))

          (0 until nSel).foreach(s => {
            assert(c.io.issueOk(s).peek().litToBoolean == expectOk(s),
              f"issueOk($s%d) is ${c.io.issueOk(s).peek().litToBoolean}%b at " +
              f"cycle $now%d, expected ${expectOk(s)}%b")
          })

          val s  = r.nextInt(nSel)
          val go = r.nextDouble() < density && expectOk(s)

          c.io.in.valid.poke(go.B)
          c.io.in.bits.poke(s.U(selW.W))
          if (go) {
            inFlight += ((s, now))
            issued   += 1
          }

          c.clock.step(1)
        }

        assert(issued > cycles / 10,
          f"only $issued%d operations were issued in $cycles%d cycles")
      }
    }
  }

  describe("ReservationTable") {

    runtest("two boundaries, every latency of its own", 2, 4,
      Seq(Map(1 -> 1, 2 -> 3, 3 -> 4), Map(1 -> 7, 2 -> 5, 3 -> 6)))

    runtest("two boundaries, latencies shared between selection values", 3, 6,
      Seq(Map(1 -> 2, 2 -> 2, 3 -> 5, 4 -> 5, 5 -> 5),
          Map(1 -> 8, 2 -> 9, 3 -> 8, 4 -> 9, 5 -> 8)))

    runtest("a boundary an operation reaches in no cycle at all", 2, 4,
      Seq(Map(1 -> 0, 2 -> 0, 3 -> 2), Map(1 -> 3, 2 -> 4, 3 -> 4)))

    runtest("a selection value one boundary alone knows", 3, 5,
      Seq(Map(1 -> 2, 3 -> 3), Map(1 -> 6, 2 -> 4, 3 -> 6, 4 -> 5)))

    runtest("one boundary, an empty one beside it", 2, 3,
      Seq(Map.empty[Int, Int], Map(1 -> 3, 2 -> 1)))

    runtest("a single boundary every operation reaches at once", 2, 3,
      Seq(Map(1 -> 0, 2 -> 0)))
  }
}
