//% @file DynMul.scala
//
// Integer multiplier whose element width is selected at run time.
//
// part of RIAL - RIKEN Arithmetic and Logic library
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.dyn

import chisel3._
import chisel3.util._

import rial.util._

import scala.collection.mutable.ArrayBuffer

// All the modes share one bit heap over the whole input width, whatever the
// rows of that heap are. A row may only see the bits of one element, so every
// bit of it is gated by whether the lanes it comes from belong to the same
// element, and a pair that crosses an element boundary contributes nothing:
// what is left of the heap is one independent multiplier per element. Whether
// two lanes belong to the same element is a comparison of the mode against a
// constant, as `heads` is.
//
// The heap is summed in two steps. It is first reduced to two rows by full
// adders, then the two rows are added by one carry chain. Both steps are cut at
// the lanes that start an element: the reduction drops a carry that would leave
// an element and the carry chain forces its carry-in there, so every element's
// product stays inside its own window, i.e. it is taken modulo 2^(2*elemW),
// which is exactly the two's complement product.
//
// A bit of negative weight is never subtracted. It is replaced by its inverse,
// which costs `-2^c` at that weight, and all that is owed is added back as one
// constant row selected by the mode. Which bit is negatively weighted depends
// on the mode - the sign bit of an operand is the topmost bit of a lane that
// ends an element - so the inversions are driven by the lane flags.
//
// Two ways of building the heap are available.
//
// The plain array has one row per bit of the multiplier: the partial product at
// weight 2^(i+j) is `a(i) & b(j)`. Signed elements are handled the
// Baugh-Wooley way, where a partial product with exactly one sign bit among its
// two bits is the negatively weighted one. Only the topmost row and column of a
// lane block can hold such a pair, hence 2*laneW-1 inversions per block.
//
// The radix-4 Booth recoding has one row per two bits of the multiplier: the
// three bits (b(p+1), b(p), b(p-1)) at every even position p give a digit in
// -2..2, and the row of that digit is the digit times the multiplicand, at
// weight 2^p. Half the rows cost, per row, the multiplexing of the digit, a
// `+1` where the digit is negative, and the sign bits of a row that is now
// signed whether the elements are or not. Since a digit reads the bit below its
// position, laneW must be even for a digit to stay inside one lane.
//
// The product of an element is twice as wide as the element, so the output is
// twice as wide as the input and every lane of the input owns two lanes of the
// output. With laneW=8 and nLanes=4:
//
//   operand lane     3       2       1       0
//                 +-------+-------+-------+-------+
//   mode=0        |   a3  |   a2  |   a1  |   a0  |        8 bit x 4
//   mode=1        |       a1      |       a0      |       16 bit x 2
//   mode=2        |               a0              |       32 bit x 1
//
//   product window   3       2       1       0     (each twice as wide)
//                 +-------+-------+-------+-------+
//   mode=0        | a3*b3 | a2*b2 | a1*b1 | a0*b0 |       16 bit x 4
//   mode=1        |     a1*b1     |     a0*b0     |       32 bit x 2
//   mode=2        |             a0*b0             |       64 bit x 1

/** Integer multiplier with a run-time selectable element width
 *
 * {{{
 * class DynMul(...) extends Module {
 *   val io = IO(new Bundle{
 *     val in  = Flipped(Valid(new DynMul.Input(laneW, nLanes)))
 *     val out = Valid(new DynMul.Output(laneW, nLanes))
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
 * `prod` is twice as wide as the operands and holds the full product of each
 * element: the product of an element of `elemW` bits is `2*elemW` bits wide and
 * sits at the window the element itself occupies, doubled. With the default
 * parameters mode=0 gives four 16 bit products, mode=1 two 32 bit products and
 * mode=2 one 64 bit product. A full product never overflows, so the result is
 * exact in every mode.
 *
 * `booth` selects the radix-4 recoding of the multiplier, which takes two bits
 * of it per partial product instead of one and so halves the number of them; it
 * requires an even `laneW`. The result is the same either way.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"reduce"` after the heap has been reduced to two rows
 *  - `"out"`    the module output
 *
 * [[getStage]] gives the resulting latency, and `out.valid` is `in.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // two registers, placed by the module
 * new DynMul(8, 4, true, 2)
 * // hand-picked positions
 * new DynMul(8, 4, true, PipelineStageConfig.specify("reduce" -> 1, "out" -> 1))
 * // combinational
 * new DynMul(8, 4, true, PipelineStageConfig.none)
 * // one partial product per bit of the multiplier
 * new DynMul(8, 4, true, PipelineStageConfig.none, false)
 * }}}
 *
 * @constructor construct DynMul module
 * @param laneW  width of the narrowest element
 * @param nLanes number of lanes. must be a power of 2
 * @param signed if true, elements are multiplied as signed integers
 * @param stage  pipeline registers of this module
 * @param booth  if true, the multiplier is recoded radix-4. needs an even laneW
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class DynMul(
  val laneW  : Int = 8,
  val nLanes : Int = 4,
  val signed : Boolean = true,
  val stage  : PipelineStageConfig = PipelineStageConfig.none,
  val booth  : Boolean = true
) extends Module {

  require(laneW   >= 1,   s"laneW must be positive, but got ${laneW}")
  require(nLanes  >= 1,   s"nLanes must be positive, but got ${nLanes}")
  require(isPow2(nLanes), s"nLanes must be a power of 2, but got ${nLanes}")
  require(!booth || laneW % 2 == 0,
    s"booth needs an even laneW so that a radix-4 digit stays inside a lane, " +
    s"but got ${laneW}")

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = DynMul.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val width = laneW * nLanes
  val modeW = DynLane.modeW(nLanes)
  val nCol  = 2 * width

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new DynMul.Input(laneW, nLanes)))
    val out = Valid(new DynMul.Output(laneW, nLanes))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val aIn    = cutIn(io.in.bits.a)
  val bIn    = cutIn(io.in.bits.b)
  val modeIn = cutIn(io.in.bits.mode)

  //----------------------------------------------------------------------
  // S1 : the heap of partial products, reduced to two rows

  val heads = DynLane.heads(modeIn, nLanes)
  val tails = DynLane.tails(heads)

  // lanes li and lj belong to the same element iff li>>mode == lj>>mode, i.e.
  // iff li^lj < 2^mode, i.e. iff mode > floor(log2(li^lj)), which is a
  // comparison against a constant
  val same = Seq.tabulate(nLanes, nLanes)((li, lj) =>
    if (li == lj) {
      true.B
    } else {
      modeIn > (31 - Integer.numberOfLeadingZeros(li ^ lj)).U
    })

  // bit i-1 of an operand, the one below bit i: nothing below the lowest bit
  // and, at a lane boundary, what the lane below holds only if that lane
  // belongs to the same element as lane i does
  private def below(x: UInt, i: Int): Option[Bool] = {
    if (i == 0) {
      None
    } else if (i % laneW == 0) {
      Some(x(i-1) && !heads(i / laneW))
    } else {
      Some(x(i-1))
    }
  }

  //---------------------------- the plain array -------------------------

  // bit i of an operand can be the sign bit of its element only if it is the
  // topmost bit of its lane; it is one iff that lane also ends the element
  private def signPos(i: Int): Boolean = signed && i % laneW == laneW-1

  // the partial product of bit i of a and bit j of b is negatively weighted
  // iff exactly one of the two bits is a sign bit of its element. None where
  // that cannot happen whatever the mode is.
  private def negative(i: Int, j: Int): Option[Bool] = {
    val li = i / laneW
    val lj = j / laneW
    val sa = if (signPos(i)) { Some(tails(li)) } else { None }
    val sb = if (signPos(j)) { Some(tails(lj)) } else { None }
    // two distinct lanes of one element cannot both end it, so a pair of
    // topmost bits is negatively weighted iff exactly one of its lanes does
    val neg = (sa, sb) match {
      case (Some(x), Some(y)) => { if (li == lj) { None } else { Some(x ^ y) } }
      case (Some(x), None)    => { Some(x) }
      case (None,    Some(y)) => { Some(y) }
      case (None,    None)    => { None }
    }
    neg.map(n => if (li == lj) { n } else { n && same(li)(lj) })
  }

  // the partial product bit at weight 2^(i+j), inverted where the pair is
  // negatively weighted and zero where the two bits belong to different
  // elements
  private def partial(i: Int, j: Int): Bool = {
    val li = i / laneW
    val lj = j / laneW
    val bit = if (li == lj) {
      aIn(i) && bIn(j)
    } else {
      aIn(i) && bIn(j) && same(li)(lj)
    }
    negative(i, j) match {
      case Some(n) => { bit ^ n }
      case None    => { bit }
    }
  }

  // what the inverted partial products of mode m owe: per element the sum of
  // the -2^(i+j) they cost, taken modulo 2^(2*elemW) and placed at the
  // element's own window of the output. For nLanes = 1 this is the usual
  // 2^laneW + 2^(2*laneW-1).
  private def owedByArray(m: Int): BigInt = {
    if (!signed) {
      BigInt(0)
    } else {
      val elemW = laneW << m
      val wrap  = BigInt(1) << (2 * elemW)
      val perElem = (0 until elemW).flatMap(i => (0 until elemW).map(j => (i, j)))
        .filter{ case(i, j) => (i == elemW-1) != (j == elemW-1) }
        .map   { case(i, j) => -(BigInt(1) << (i + j)) }
        .sum
      wrapModes(perElem, m)
    }
  }

  private def arrayHeap: (Seq[Seq[Bool]], Seq[BigInt]) = {
    val columns = Seq.fill(nCol)(new ArrayBuffer[Bool])
    for(i <- 0 until width; j <- 0 until width) {
      columns(i+j) += partial(i, j)
    }
    (columns.map(_.toSeq), owedOfEveryMode(owedByArray))
  }

  //------------------------- the radix-4 recoding -----------------------

  // what the inverted top bits of a recoded mode m owe: per element the
  // -2^(p + elemW) of each of its digits, one weight higher where the elements
  // are unsigned because a row then carries one more bit
  private def owedByBooth(m: Int): BigInt = {
    val elemW = laneW << m
    val top   = if (signed) { elemW } else { elemW + 1 }
    wrapModes((0 until elemW by 2).map(p => -(BigInt(1) << (p + top))).sum, m)
  }

  private def boothHeap: (Seq[Seq[Bool]], Seq[BigInt]) = {
    val columns = Seq.fill(nCol)(new ArrayBuffer[Bool])

    // first(lh)(lp) is true iff lane lh is the lane that starts the element
    // containing lane lp, i.e. it is a head and no lane between the two is.
    // Exactly one lane of an element satisfies it, lane 0 being a head always.
    val first = Seq.tabulate(nLanes, nLanes)((lh, lp) =>
      if (lh > lp) {
        false.B
      } else {
        ((lh+1) to lp).foldLeft(heads(lh))((acc, l) => acc && !heads(l))
      })

    for(p <- 0 until width by 2) {
      val lp = p / laneW
      val b0 = bIn(p)
      val b1 = bIn(p+1)

      // the digit of the three bits (b(p+1), b(p), b(p-1)) is -2*b1 + b0 +
      // b(p-1): `one` and `two` select its magnitude, a or 2a, and `neg` asks
      // for the opposite sign. The digit 111 is a zero, which never negates.
      val (neg, one, two) = below(bIn, p) match {
        case Some(bm) => {
          (b1 && !(b0 && bm), b0 ^ bm,
            (b1 && !b0 && !bm) || (!b1 && b0 && bm))
        }
        case None => { (b1, b0, b1 && !b0) }
      }

      // the row is the digit times the multiplicand of the same element, so a
      // bit of a that belongs to another element contributes nothing. A
      // negative digit takes the inverse of the magnitude, whose `+1` is added
      // at the lowest bit of the element the row belongs to.
      for(i <- 0 until width) {
        val la  = i / laneW
        val src = below(aIn, i) match {
          case Some(ap) => { (one && aIn(i)) || (two && ap) }
          case None     => { one && aIn(i) }
        }
        val bit = src ^ neg
        columns(i + p) += (if (la == lp) { bit } else { bit && same(la)(lp) })
      }
      for(lh <- 0 to lp) {
        columns(laneW * lh + p) += (neg && first(lh)(lp))
      }

      // the topmost bits of the row sit above the element, so which weight
      // they take depends on where the element ends
      for(lt <- lp until nLanes) {
        val hi = laneW * (lt + 1) - 1
        val ends = if (lt == lp) {
          tails(lt)
        } else {
          tails(lt) && same(lt)(lp)
        }
        if (signed) {
          // the row of a signed multiplicand reaches bit hi+2 and its two
          // topmost bits are both the sign, so the weights they carry add up
          // to the single negative weight -2^(p+hi+1)
          val sign = ((one || two) && aIn(hi)) ^ neg
          columns(p + hi + 1) += (!sign && ends)
        } else {
          // an unsigned multiplicand has no sign bit, so 2a reaches bit hi+1
          // and the row is negative only where the digit asks for it
          columns(p + hi + 1) += (((two && aIn(hi)) ^ neg) && ends)
          columns(p + hi + 2) += (!neg && ends)
        }
      }
    }

    // the recoding reads the multiplier as a signed number, so an unsigned one
    // needs one more digit, of value b(hi) at the weight just above the
    // element. Its row is never negative.
    if (!signed) {
      for(lt <- 0 until nLanes) {
        val hi = laneW * (lt + 1) - 1
        for(i <- 0 to hi) {
          val la = i / laneW
          val ends = if (la == lt) {
            tails(lt)
          } else {
            tails(lt) && same(la)(lt)
          }
          columns(i + hi + 1) += (bIn(hi) && aIn(i) && ends)
        }
      }
    }

    (columns.map(_.toSeq), owedOfEveryMode(owedByBooth))
  }

  //--------------------------- the shared part --------------------------

  // one element's worth of what is owed, taken modulo the element's own width
  // and repeated at the window of every element of mode m
  private def wrapModes(perElem: BigInt, m: Int): BigInt = {
    val elemW   = laneW << m
    val wrap    = BigInt(1) << (2 * elemW)
    val wrapped = ((perElem % wrap) + wrap) % wrap
    (0 until (nLanes >> m)).map(e => wrapped << (2 * elemW * e)).sum
  }

  private def owedOfEveryMode(owed: Int => BigInt): Seq[BigInt] = {
    Seq.tabulate(DynLane.modeMax(nLanes)+1)(m => owed(m))
  }

  private val (heap, owedOfMode) = if (booth) { boothHeap } else { arrayHeap }

  // what is owed is negative at every weight it covers, so it is added back as
  // one more row of the heap, selected by the mode
  val columns = if (owedOfMode.forall(_ == 0)) {
    heap
  } else {
    val owedRow = if (owedOfMode.length == 1) {
      owedOfMode.head.U(nCol.W)
    } else {
      MuxLookup(modeIn, 0.U(nCol.W))(
        owedOfMode.zipWithIndex.map{ case(v, m) => (m.U, v.U(nCol.W)) })
    }
    Seq.tabulate(nCol)(c => {
      if (owedOfMode.exists(_.testBit(c))) {
        heap(c) :+ owedRow(c)
      } else {
        heap(c)
      }
    })
  }

  // the heap is as wide as the output, so its lanes are twice as wide as the
  // operand lanes while the lanes that start an element are the same ones
  val (row0Comb, row1Comb) = DynLane.compress(columns, 2 * laneW, heads)

  //--------------------------- stage "reduce" ---------------------------
  private val cutReduce = new SharedPipe("reduce", cutIn.valid, stages.at("reduce"))

  val row0       = cutReduce(row0Comb)
  val row1       = cutReduce(row1Comb)
  val modeReduce = cutReduce(modeIn)

  //----------------------------------------------------------------------
  // S2 : the carry chain that adds the two rows, cut at the same lanes. the
  // carry-out of a lane belongs to the next element, so only the sum is kept.

  val headsReduce = DynLane.heads(modeReduce, nLanes)

  val prodComb = DynLane.concat(DynLane.add(
    DynLane.split(row0, 2 * laneW, nLanes),
    DynLane.split(row1, 2 * laneW, nLanes),
    headsReduce).map(s => s(2*laneW-1, 0)))

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutReduce.valid, stages.at("out"))

  io.out.bits.prod := cutOut(prodComb)
  io.out.valid     := cutOut.valid
}

/** Input and output bundles and pipeline stage points of [[rial.dyn.DynMul]].
 *
 * @see [[rial.util.StagePoint]]
 */
object DynMul {

  /** Input of [[rial.dyn.DynMul]]: the two operands and the element width.
   *
   * @param laneW  width of the narrowest element
   * @param nLanes number of lanes
   */
  class Input(laneW: Int, nLanes: Int) extends Bundle {
    val a    = UInt((laneW * nLanes).W)
    val b    = UInt((laneW * nLanes).W)
    val mode = UInt(DynLane.modeW(nLanes).W)
  }

  /** Output of [[rial.dyn.DynMul]]: the full product of each element, twice as
   *  wide as the operands.
   *
   * @param laneW  width of the narrowest element
   * @param nLanes number of lanes
   */
  class Output(laneW: Int, nLanes: Int) extends Bundle {
    val prod = UInt((2 * laneW * nLanes).W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.dyn.DynMul]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("reduce", 4), StagePoint("out", 2))
}
