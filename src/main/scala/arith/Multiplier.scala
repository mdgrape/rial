//% @file Multiplier.scala
//
// Integer multiplier with a selectable internal structure, shared by every
// module of the library that needs a product.
//
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.arith

import chisel3._
import chisel3.util._
import chisel3.util.experimental.InlineInstance

import rial.util._

/** Product of two integers, with the structure [[rial.arith.MulImpl]] names.
 *
 * `io.out.bits` is the product of `io.in.bits.x` and `io.in.bits.y`, unsigned
 * by default and two's complement when `signed` is set (the operands, the
 * addend and the result are then read with `asSInt`).
 *
 * With `shift` the product is divided by `2^shift`, rounded down; every
 * structure takes it, and a module that wants no shift leaves it out.
 *
 * `MulImpl.macc` makes the product a multiply-add: the input then carries an
 * addend `io.in.bits.c` of `cW` bits, and the result is
 *
 *   `io.out.bits = (x * y >> shift) + c`
 *
 * so the addend is at the weight of the shifted result and a caller hands over
 * the number it has. [[Multiplier.outWidth]] is the width of the result. The
 * other structures have no addend and must not be given one.
 *
 * A module that needs a product instantiates this one and does nothing else
 * about it: what a wide operand costs, what a two's complement one costs once
 * it is sliced, and where the addend has to sit for the product to stay one
 * multiply-add, is this module's business. `MulImpl.plain` is the plain
 * product and asks nothing of the operands.
 *
 * A multiply-add is one multiplication where the narrower operand fits a
 * slice, and one per slice of the wider operand where it does not; the addend
 * is one more term of the sum that adds them up. The partial products and the
 * addend are what the boundary `"pp"` holds, and the sum is behind it, so a
 * register at `"pp"` leaves each multiplication on the near side of it, which
 * is where a hard multiplier keeps its own.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"pp"`     after the partial products, only with `MulImpl.split` and
 *               `MulImpl.macc`; with `macc` the addend crosses it beside them,
 *               and the sum that adds them all is behind it
 *  - `"chainK"` after the k-th partial sum, `K` = 1 ... one less than the
 *               number of hops, so only when a plain unsigned product is split
 *               into more than two slices
 *  - `"out"`    the module output, after the sum of the partial products
 *
 * [[getStage]] gives the resulting latency, and `io.out.valid` is
 * `io.in.valid` delayed by that many cycles. The module is inlined into the
 * one that instantiates it, so a product is the same netlist wherever it is
 * written and a hierarchy boundary never sits between it and its neighbours.
 *
 * {{{
 * // one product, two registers at the output
 * new Multiplier(24, 24, MulImpl.plain, PipelineStageConfig.atOut(2))
 * // one partial product per 17 bit slice of x, registered, then summed
 * new Multiplier(24, 24, MulImpl.split(17),
 *   PipelineStageConfig.specify("pp" -> 1, "out" -> 1))
 * // a two's complement product, combinational
 * new Multiplier(19, 24, MulImpl.split(17), PipelineStageConfig.none, true)
 * // a two's complement multiply-add, its product divided by 2^18 and its
 * // addend at the weight of that, its partial products and its addend
 * // registered
 * new Multiplier(19, 24, MulImpl.macc(17),
 *   PipelineStageConfig.specify("pp" -> 1), true, Some(25), Some(18))
 * }}}
 *
 * @constructor create a chisel Module Multiplier
 * @param xW     width of the left hand operand
 * @param yW     width of the right hand operand
 * @param impl   structure of the product
 * @param stage  pipeline registers of this module
 * @param signed read both operands and the result as two's complement
 * @param cW     width of the addend, which `MulImpl.macc` has and the other
 *               structures do not; the addend is at the weight of the result,
 *               after the shift
 * @param shift  bits the product drops at its bottom, when it drops any
 *
 * @see [[rial.arith.MulImpl]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class Multiplier(
  xW:     Int,
  yW:     Int,
  impl:   MulImpl = MulImpl.plain,
  stage:  PipelineStageConfig = PipelineStageConfig.none,
  signed: Boolean = false,
  cW:     Option[Int] = None,
  shift:  Option[Int] = None,
) extends Module with InlineInstance {

  /** Number of boundaries the chain of partial sums offers here. */
  val nChain = Multiplier.chainPoints(impl, xW, yW, signed)

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = Multiplier.stagePoints(impl, nChain)

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  /** Bits the result drops at its bottom. */
  val shiftW = shift.getOrElse(0)

  /** True when the narrower operand fits a slice, so that a multiply-add is
   *  one multiplication and the wider operand is not sliced.
   */
  val oneMultiply = Multiplier.fitsSlice(xW.min(yW), impl)

  val outW = Multiplier.outWidth(xW, yW, cW, shiftW)

  /** Width of the product on its own. */
  private val prodW = xW + yW

  /** Width at which a multiply-add sums its terms: one bit above both the
   *  product and the addend, so that the sum is exact before the shift.
   */
  private val sumW = outW + shiftW

  require(!impl.isMacc || cW.isDefined,
    "a macc multiplier is a multiply-add and needs the width of its addend")
  require(impl.isMacc || cW.isEmpty,
    s"a ${impl.kind} multiplier has no addend; ask for MulImpl.macc to get one")
  require(shiftW >= 0,
    s"a multiplier drops no negative number of bits, but shift is $shiftW")
  require(outW > 0,
    s"a multiplier of $xW by $yW bits with an addend of ${cW.getOrElse(0)} " +
    s"bits has nothing left after a shift of $shiftW")

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new Multiplier.Input(xW, yW, cW)))
    val out = Valid(UInt(outW.W))
  })

  // --------------------------- stage "in" ---------------------------
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val xIn = cutIn(io.in.bits.x)
  val yIn = cutIn(io.in.bits.y)
  val cIn = io.in.bits.c.map(c => cutIn(c))

  // --------------------------- stage "pp" ---------------------------
  private val nPp = if (impl.isSliced) {
    stages.at("pp")
  } else {
    0
  }

  private val cutPp = new SharedPipe("pp", cutIn.valid, nPp)

  /** The cut the sum of the partial products ends at. */
  private var cutLast: SharedPipe = cutPp

  /** The product of the partial products [[Multiplier.partials]] returns, each
   *  adding its own to what came before it, shifted down by the slice width;
   *  the boundary `"chainK"` sits between two of those additions.
   */
  private def chained(parts: Seq[UInt]): UInt = {
    if (parts.length == 1) {
      Multiplier.fit(parts.head, prodW)
    } else {
      val sw = impl.sliceW.get
      var carried = parts.head
      var rest    = parts.tail
      var low     = Seq.empty[UInt]
      var k       = 1
      while (rest.nonEmpty) {
        val next = rest.head +& (carried >> sw)
        val done = low :+ carried(sw - 1, 0)
        rest = rest.tail
        if (rest.nonEmpty) {
          val cut = new SharedPipe(s"chain$k", cutLast.valid, stages.at(s"chain$k"))
          carried = cut(next)
          low     = done.map(l => cut(l))
          rest    = rest.map(r => cut(r))
          cutLast = cut
          k      += 1
        } else {
          carried = next
          low     = done
        }
      }
      Multiplier.fit(Cat(carried, Cat(low.reverse)), prodW)
    }
  }

  /** `v` divided by `2^shift`, rounded down the way the module reads it. */
  private def dropped(v: UInt): UInt = {
    if (shiftW == 0) {
      v
    } else if (signed) {
      (v.asSInt >> shiftW).asUInt
    } else {
      v >> shiftW
    }
  }

  /** What it costs to read two's complement operands unsigned in the slices,
   *  `sumW` bits wide and to be subtracted from the sum of the slices.
   *
   * With `xu` and `yu` the operands read unsigned, `x = xu - x[xW-1] * 2^xW`
   * and `y = yu - y[yW-1] * 2^yW` give
   *
   *   x * y = xu*yu - x[xW-1] * y * 2^xW - y[yW-1] * xu * 2^yW
   *
   * exactly, at every width: the left hand term keeps the sign of `y` and so
   * carries the corner where both operands are negative. Both terms are
   * selections of an operand, so their sum is formed while the partial
   * products are.
   */
  private def signCost(): UInt = {
    val m = xW.min(yW)
    val cost = (Mux(xIn(xW-1), yIn.asSInt, 0.S(yW.W)) << (xW - m)) +&
               (Mux(yIn(yW-1), xIn, 0.U).zext << (yW - m))
    Multiplier.fit((cost.pad(sumW) << m).asUInt, sumW)
  }

  /** The addend at the weight of the product, `sumW` bits wide, which is where
   *  it joins the product: `(x * y + c * 2^shift) >> shift` is the same number
   *  as `(x * y >> shift) + c`, and the addend reaches the sum of the partial
   *  products without an adder of its own.
   */
  private def atProduct(c: UInt): UInt = {
    if (signed) {
      Multiplier.fit((c.asSInt.pad(sumW) << shiftW).asUInt, sumW)
    } else {
      Multiplier.fit(c << shiftW, sumW)
    }
  }

  /** The product of the two operands, in the structure `impl` names.
   *
   * A slice of an operand is what does not fit one multiplication, and the
   * slices are unsigned, so a two's complement operand is read unsigned in
   * them and what that costs comes back: `x = xu - x[xW-1] * 2^xW` gives
   *
   *   x * y = xu*yu - (x < 0 ? yu : 0) * 2^xW - (y < 0 ? xu : 0) * 2^yW
   *
   * modulo 2^(xW+yW). The two terms beside the product are selections of an
   * operand, so their sum is formed while the partial products are and then
   * rides on the addend input of the least significant multiplier, which is
   * where a hard multiplier subtracts a term for free.
   */
  private def product(): UInt = {
    if (!impl.isSliced) {
      // the product as it would be written by hand: `MulImpl.plain` asks
      // nothing of the operands, so there is nothing to shape
      dropped(if (signed) {
        (xIn.asSInt * yIn.asSInt).asUInt
      } else {
        xIn * yIn
      })
    } else if (!signed) {
      // one partial product per slice, each added to what came before it: an
      // unsigned `MulImpl.split`
      dropped(chained(Multiplier.partials(xIn, yIn, impl).map(p => cutPp(p))))
    } else {
      // the partial products and the cost of their sign in one sum: a two's
      // complement `MulImpl.split`
      val parts = Multiplier.partials(xIn, yIn, impl)
      val m   = xW.min(yW)
      val sub = (Mux(xIn(xW-1), yIn, 0.U) << (xW - m)) +&
                (Mux(yIn(yW-1), xIn, 0.U) << (yW - m))
      val head = Multiplier.fit(parts.head, prodW) - Cat(
        Multiplier.fit(sub, prodW - m), 0.U(m.W))
      dropped(Multiplier.sumWide(
        (head +: Multiplier.placed(parts.tail, impl, prodW, 1)).map(t => cutPp(t)),
        prodW))
    }
  }

  /** The product of the two operands with the addend `c`, which is at the
   *  weight of the shifted result: `(x * y >> shift) + c`.
   */
  private def multiplyAdd(c: UInt): UInt = {
    // the operands are read unsigned, so that what a multiplication is given
    // is a slice and what it costs to read a two's complement operand that way
    // is one term of the sum; the addend is one more term of it, and the sum
    // is exact at `sumW` bits
    val parts = if (oneMultiply) {
      Seq(xIn * yIn)
    } else {
      Multiplier.partials(xIn, yIn, impl)
    }
    val atWeight = Multiplier.fit(parts.head, sumW) +:
      Multiplier.placed(parts.tail, impl, sumW, 1)
    val terms = if (signed) {
      (atWeight.head - signCost()) +: atWeight.tail
    } else {
      atWeight
    }
    dropped(Multiplier.sumWide(
      (terms :+ atProduct(c)).map(t => cutPp(t)), sumW))
  }

  val prod = cIn match {
    case Some(c) => {
      multiplyAdd(c)
    }
    case None => {
      product()
    }
  }

  require(prod.getWidth == outW,
    s"product is ${prod.getWidth} bits, expected $outW")

  // --------------------------- stage "out" --------------------------
  io.out := Pipe(cutLast.valid, prod, stages.at("out"))
}

/** Input bundle and pipeline stage points of [[rial.arith.Multiplier]].
 *
 * @see [[rial.util.StagePoint]]
 */
object Multiplier {

  /** Input of [[rial.arith.Multiplier]]: the two factors, and the addend of a
   *  multiply-add when the module was asked for one.
   *
   * @param xW width of the left hand operand
   * @param yW width of the right hand operand
   * @param cW width of the addend, when there is one
   */
  class Input(xW: Int, yW: Int, cW: Option[Int] = None) extends Bundle {
    val x = UInt(xW.W)
    val y = UInt(yW.W)
    val c = cW.map(w => UInt(w.W))
  }

  /** The `cW` a module that has no addend of its own passes, which is the
   *  narrowest one `MulImpl.macc` accepts and nothing at all for a structure
   *  that has no addend. Such a module drives `io.in.bits.c` with zero and
   *  reads the low `xW + yW` bits of the result, which are the product.
   *
   * @param impl structure of the product
   */
  def zeroAddendW(impl: MulImpl): Option[Int] = {
    if (impl.isMacc) {
      Some(1)
    } else {
      None
    }
  }

  /** Width of the result of [[rial.arith.Multiplier]] for these parameters,
   *  which is one bit above both the shifted product and the addend.
   *
   * @param xW    width of the left hand operand
   * @param yW    width of the right hand operand
   * @param cW    width of the addend, when there is one
   * @param shift bits the product drops at its bottom
   */
  def outWidth(xW: Int, yW: Int, cW: Option[Int], shift: Int): Int = {
    cW match {
      case Some(w) => {
        (xW + yW - shift).max(w) + 1
      }
      case None => {
        xW + yW - shift
      }
    }
  }

  /** True when an operand of `w` bits fits one slice of `impl`, which is
   *  `MulImpl.sliceW` bits, so that a multiply-add of such an operand is one
   *  multiplication. A structure that slices nothing fits every operand.
   *
   * @param w    width of the operand
   * @param impl structure of the product
   */
  def fitsSlice(w: Int, impl: MulImpl): Boolean = {
    impl.sliceW match {
      case Some(sw) if impl.isSliced => {
        w <= sw
      }
      case _ => {
        true
      }
    }
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.arith.Multiplier]]. `"pp"` exists only with the structures that
   *  slice, which are `MulImpl.split` and `MulImpl.macc`, and the `nChain`
   *  boundaries `"chain1"` ... `"chain<nChain>"` only when the sum of the
   *  partial products is a chain of more than one addition; [[chainPoints]] is
   *  how many that is.
   *
   * @param impl   structure of the product
   * @param nChain number of boundaries inside the sum
   */
  def stagePoints(impl: MulImpl, nChain: Int = 0): Seq[StagePoint] = {
    if (impl.isSliced) {
      (StagePoint("pp", 3.0) +:
        (1 to nChain).map(k => StagePoint(s"chain$k", 2.0))) :+
        StagePoint("out", 3.0)
    } else {
      Seq(StagePoint("out", 4.0))
    }
  }

  /** Number of boundaries the sum of the partial products offers for these
   *  parameters, which is one less than the number of additions it takes.
   *
   * A product that is not split, or split in two, is one addition and has
   * none; a two's complement product and a multiply-add are summed all at once
   * rather than in a chain, and have none either.
   *
   * @param impl   structure of the product
   * @param xW     width of the left hand operand
   * @param yW     width of the right hand operand
   * @param signed the operands are two's complement
   */
  def chainPoints(impl: MulImpl, xW: Int, yW: Int,
    signed: Boolean = false): Int = {

    if (signed || impl.isMacc) {
      0
    } else {
      (sliceBits(xW.max(yW), impl).length - 2).max(0)
    }
  }

  /** The slices `impl` splits an operand of `w` bits into, least significant
   *  first, as (high bit, low bit) pairs.
   */
  private def sliceBits(w: Int, impl: MulImpl): Seq[(Int, Int)] = {
    impl.sliceW match {
      case Some(sw) if impl.isSliced && sw < w => {
        val n = (w + sw - 1) / sw
        (0 until n).map(i => ((w - 1).min((i + 1) * sw - 1), i * sw))
      }
      case _ => {
        Seq((w - 1, 0))
      }
    }
  }

  /** `x * y` as one unsigned partial product per slice of the wider operand
   *  (the left hand one on a tie), least significant slice first, so that each
   *  of them is one hard multiplier of `MulImpl.sliceW` bits on its narrow
   *  port and a pipeline register on it lands inside that multiplier.
   *
   * When the wider operand fits one slice this is the single product `x * y`.
   */
  private def partials(x: UInt, y: UInt, impl: MulImpl): Seq[UInt] = {
    if (x.getWidth >= y.getWidth) {
      sliceBits(x.getWidth, impl).map({ case (hi, lo) => x(hi, lo) * y })
    } else {
      sliceBits(y.getWidth, impl).map({ case (hi, lo) => x * y(hi, lo) })
    }
  }

  /** The partial products [[partials]] returns, each at its own weight and
   *  `outW` bits wide, for a sum that has other terms in it. `first` is the
   *  slice index of `parts.head`.
   */
  private def placed(parts: Seq[UInt], impl: MulImpl, outW: Int,
    first: Int): Seq[UInt] = {

    parts.zipWithIndex.map({ case (p, i) =>
      val sh = (first + i) * impl.sliceW.get
      Cat(fit(p, outW - sh), 0.U(sh.W))
    })
  }

  /** The low `w` bits of `x`, zero extended if it is narrower. */
  private def fit(x: UInt, w: Int): UInt = {
    x.pad(w)(w - 1, 0)
  }

  /** The low `w` bits of the sum of `ops`. Three operands at a time become two
   *  by a carry-save step, so however many there are the sum ends in one carry
   *  chain.
   */
  private def sumWide(ops: Seq[UInt], w: Int): UInt = {
    var cur = ops.map(o => fit(o, w))
    while (cur.length > 2) {
      val a = cur(0)
      val b = cur(1)
      val c = cur(2)
      val s = a ^ b ^ c
      val carry = (((a & b) | (b & c) | (c & a)) << 1)(w - 1, 0)
      cur = s +: carry +: cur.drop(3)
    }
    cur.reduce(_ + _)
  }
}
