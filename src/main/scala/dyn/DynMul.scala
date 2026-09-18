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

// All the modes share one array multiplier over the whole input width. Its
// partial product at weight 2^(i+j) is `a(i) & b(j)`, gated by whether bit i
// and bit j belong to the same element: a pair that crosses an element
// boundary contributes nothing, so what is left of the array is one
// independent multiplier per element. Whether two lanes belong to the same
// element is a comparison of the mode against a constant, as `heads` is.
//
// The array is summed in two steps. The bit heap is first reduced to two rows
// by full adders, then the two rows are added by one carry chain. Both steps
// are cut at the lanes that start an element: the reduction drops a carry that
// would leave an element and the carry chain forces its carry-in there, so
// every element's product stays inside its own window, i.e. it is taken modulo
// 2^(2*elemW), which is exactly the two's complement product.
//
// Signed elements are handled the Baugh-Wooley way. A partial product with
// exactly one sign bit among its two bits is negatively weighted; it is
// replaced by its inverse, which costs `-2^(i+j)` per such weight, and all
// that is owed is added back as one constant row. Which bit of an operand is a
// sign bit depends on the mode - it is the topmost bit of a lane that ends an
// element - so the inversions are driven by the lane flags and the constant is
// selected by the mode. Only the topmost row and column of a lane block can
// hold a negative partial product, hence 2*laneW-1 inversions per block.
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
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"reduce"` after the array has been reduced to two rows
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
 * }}}
 *
 * @constructor construct DynMul module
 * @param laneW  width of the narrowest element
 * @param nLanes number of lanes. must be a power of 2
 * @param signed if true, elements are multiplied as signed integers
 * @param stage  pipeline registers of this module
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class DynMul(
  val laneW  : Int = 8,
  val nLanes : Int = 4,
  val signed : Boolean = true,
  val stage  : PipelineStageConfig = PipelineStageConfig.none
) extends Module {

  require(laneW   >= 1,   s"laneW must be positive, but got ${laneW}")
  require(nLanes  >= 1,   s"nLanes must be positive, but got ${nLanes}")
  require(isPow2(nLanes), s"nLanes must be a power of 2, but got ${nLanes}")

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = DynMul.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val width = laneW * nLanes
  val modeW = DynLane.modeW(nLanes)

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
  // S1 : the array of partial products, reduced to two rows

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
  private def owed(m: Int): BigInt = {
    if (!signed) {
      BigInt(0)
    } else {
      val elemW = laneW << m
      val wrap  = BigInt(1) << (2 * elemW)
      val perElem = (0 until elemW).flatMap(i => (0 until elemW).map(j => (i, j)))
        .filter{ case(i, j) => (i == elemW-1) != (j == elemW-1) }
        .map   { case(i, j) => -(BigInt(1) << (i + j)) }
        .sum
      val wrapped = ((perElem % wrap) + wrap) % wrap
      (0 until (nLanes >> m)).map(e => wrapped << (2 * elemW * e)).sum
    }
  }

  val nCol = 2 * width

  val columns = Seq.fill(nCol)(new ArrayBuffer[Bool])
  for(i <- 0 until width; j <- 0 until width) {
    columns(i+j) += partial(i, j)
  }

  val owedOfMode = Seq.tabulate(DynLane.modeMax(nLanes)+1)(m => owed(m))
  if (owedOfMode.exists(_ != 0)) {
    val constant = if (owedOfMode.length == 1) {
      owedOfMode.head.U(nCol.W)
    } else {
      MuxLookup(modeIn, 0.U(nCol.W))(
        owedOfMode.zipWithIndex.map{ case(v, m) => (m.U, v.U(nCol.W)) })
    }
    for(c <- 0 until nCol) {
      if (owedOfMode.exists(_.testBit(c))) { columns(c) += constant(c) }
    }
  }

  // the heap is as wide as the output, so its lanes are twice as wide as the
  // operand lanes while the lanes that start an element are the same ones
  val (row0Comb, row1Comb) = DynLane.compress(
    columns.map(_.toSeq), 2 * laneW, heads)

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
