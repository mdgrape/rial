//% @file DSP.scala
//
// The products that are written to fit the hard multipliers (the DSP tiles)
// of a device: one partial product per slice of an operand, and a
// multiply-add whose addend rides on the multiplier that forms the product,
// together with the structure a user of such a module names.
//
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.arith

import chisel3._
import chisel3.util._
import chisel3.util.experimental.InlineInstance

import rial.util._

/** Structure of a product, chosen by the user of the module that builds it.
 *
 * `kind` names the structure and `sliceW` carries its parameter, so that a
 * configuration read from a file can be turned into one of these without
 * knowing the set of structures:
 *
 *  - `MulImpl.plain`    one product, written as it is. It asks nothing of the
 *                       operands and the module that needs it writes it in
 *                       place.
 *  - `MulImpl.split(w)` one partial product per `w` bit slice of the wider
 *                       operand, summed afterwards, which is
 *                       [[rial.arith.SplitMultiplier]]. A pipeline register
 *                       may be placed on the partial products, which is where
 *                       a hard multiplier of at most `w` bits on its narrow
 *                       port has one.
 *  - `MulImpl.splitSmall(w)` `split`, but for two's complement operands
 *                       the slices of the wider operand are read signed, each
 *                       in its own register at `"in"`. It uses less logic and
 *                       suits a wider product, and maps onto a hard multiplier
 *                       only where `"in"` holds a register.
 *  - `MulImpl.macc(w)`  a multiply-add, which is [[rial.arith.MultiplyAdd]]:
 *                       the product carries an addend, and the addend joins
 *                       the sum of the partial products rather than an adder
 *                       of its own, which is where a hard multiplier of at
 *                       most `w` bits on its narrow port has an addend input.
 *                       Only a module that has an addend to fold, such as the
 *                       coefficient of a Horner step, takes this one.
 *  - `MulImpl.maccSmall(w)` `macc` with its slices read as `splitSmall`
 *                       reads them.
 *
 * @param kind   `"plain"`, `"split"`, `"splitSmall"`, `"macc"` or
 *               `"maccSmall"`
 * @param sliceW width of an operand slice, for every kind but `"plain"`
 *
 * @see [[rial.arith.SplitMultiplier]]
 * @see [[rial.arith.MultiplyAdd]]
 */
case class MulImpl(kind: String, sliceW: Option[Int] = None) {

  require(MulImpl.kinds.contains(kind),
    s"unknown multiplier structure '$kind'; " +
    s"known ones are ${MulImpl.kinds.mkString(", ")}")
  require(!isSliced || sliceW.exists(_ > 0),
    s"a $kind multiplier needs a positive slice width")

  /** True for the structures that slice one operand. */
  def isSplit: Boolean = {
    kind == "split" || kind == "splitSmall"
  }

  /** True for the structures that carry an addend. */
  def isMacc: Boolean = {
    kind == "macc" || kind == "maccSmall"
  }

  /** True for the structures that read the slices of a two's complement
   *  operand signed.
   */
  def isSmall: Boolean = {
    kind == "splitSmall" || kind == "maccSmall"
  }

  /** True for the structures that slice an operand. */
  def isSliced: Boolean = {
    isSplit || isMacc
  }

  /** Short name, for a generated module or file name, e.g. `"split17"`. */
  def tag: String = {
    kind + sliceW.map(_.toString).getOrElse("")
  }
}

/** Constructor shorthands for [[rial.arith.MulImpl]]. */
object MulImpl {

  /** Every structure [[rial.arith.MulImpl]] names. */
  val kinds: Seq[String] = Seq("plain", "split", "splitSmall", "macc",
    "maccSmall")

  /** One product. */
  val plain: MulImpl = MulImpl("plain")

  /** One partial product per `w` bit slice of the wider operand. */
  def split(w: Int): MulImpl = {
    MulImpl("split", Some(w))
  }

  /** [[split]] with the slices of a two's complement operand read signed. */
  def splitSmall(w: Int): MulImpl = {
    MulImpl("splitSmall", Some(w))
  }

  /** A multiply-add whose addend rides on a hard multiplier of `w` bits. */
  def macc(w: Int): MulImpl = {
    MulImpl("macc", Some(w))
  }

  /** [[macc]] with the slices of a two's complement operand read signed. */
  def maccSmall(w: Int): MulImpl = {
    MulImpl("maccSmall", Some(w))
  }

  /** Reads `"plain"`, `"split:<w>"`, `"splitSmall:<w>"`, `"macc:<w>"` or
   *  `"maccSmall:<w>"`.
   *
   * @throws IllegalArgumentException if `s` names none of them.
   */
  def parse(s: String): MulImpl = {
    val t = s.trim
    val sliced = Seq("split" -> (split _), "splitSmall" -> (splitSmall _),
      "macc" -> (macc _), "maccSmall" -> (maccSmall _)).collectFirst {
      case (k, make) if t.startsWith(k + ":") => {
        val v = t.drop(k.length + 1)
        val w = try {
          v.toInt
        } catch {
          case _: NumberFormatException => {
            throw new IllegalArgumentException(
              s"multiplier structure '$k:$v': '$v' is not an int")
          }
        }
        make(w)
      }
    }
    if (t == "plain") {
      plain
    } else {
      sliced.getOrElse(throw new IllegalArgumentException(
        s"unknown multiplier structure '$s'; " +
        s"expected plain, split:<w>, splitSmall:<w>, macc:<w> or " +
        s"maccSmall:<w>"))
    }
  }
}

/** The slices of an operand and the partial products made of them.
 *
 * A product whose operands are wider than one hard multiplier is written as
 * one multiplication per slice of an operand; these are the pieces of such a
 * product, and they are pure functions on hardware values.
 *
 * @see [[rial.arith.SplitMultiplier]]
 * @see [[rial.arith.MultiplyAdd]]
 */
object Slice {

  /** The slices of at most `sliceW` bits an operand of `w` bits is cut into,
   *  least significant first, as (high bit, low bit) pairs. An operand that
   *  already fits `sliceW` bits is one slice.
   */
  def bits(w: Int, sliceW: Int): Seq[(Int, Int)] = {
    if (sliceW < w) {
      val n = (w + sliceW - 1) / sliceW
      (0 until n).map(i => ((w - 1).min((i + 1) * sliceW - 1), i * sliceW))
    } else {
      Seq((w - 1, 0))
    }
  }

  /** `x * y` as one unsigned partial product per slice of the wider operand
   *  (the left hand one on a tie), least significant slice first, so that each
   *  of them is one hard multiplier of `sliceW` bits on its narrow port and a
   *  pipeline register on it lands inside that multiplier.
   *
   * When the wider operand fits one slice this is the single product `x * y`.
   */
  def partials(x: UInt, y: UInt, sliceW: Int): Seq[UInt] = {
    if (x.getWidth >= y.getWidth) {
      bits(x.getWidth, sliceW).map({ case (hi, lo) => x(hi, lo) * y })
    } else {
      bits(y.getWidth, sliceW).map({ case (hi, lo) => x * y(hi, lo) })
    }
  }

  /** The partial products [[partials]] returns, each at its own weight and
   *  `outW` bits wide, for a sum that has other terms in it. `first` is the
   *  slice index of `parts.head`.
   */
  def placed(parts: Seq[UInt], sliceW: Int, outW: Int,
    first: Int): Seq[UInt] = {

    parts.zipWithIndex.map({ case (p, i) =>
      val sh = (first + i) * sliceW
      Cat(fit(p, outW - sh), 0.U(sh.W))
    })
  }

  /** `x * y` of two's complement operands as one signed partial product per
   *  slice of the wider operand (the left hand one on a tie), least
   *  significant slice first, each the two's complement bits of its product.
   *  Every slice but the top one is read with a zero bit on top of it, so it
   *  is a non-negative signed number. Each slice and the narrower operand go
   *  through `hold` on their own, so that a register there holds a slice
   *  whole and its sign extension is that of a register.
   */
  def signedPartials(x: UInt, y: UInt, sliceW: Int,
    hold: UInt => UInt): Seq[UInt] = {

    val (wide, narrow) = if (x.getWidth >= y.getWidth) {
      (x, y)
    } else {
      (y, x)
    }
    val cuts = bits(wide.getWidth, sliceW)
    val narrowHeld = hold(narrow)
    cuts.zipWithIndex.map({ case ((hi, lo), i) =>
      val slice = if (i == cuts.length - 1) {
        hold(wide(hi, lo))
      } else {
        hold(Cat(0.U(1.W), wide(hi, lo)))
      }
      (slice.asSInt * narrowHeld.asSInt).asUInt
    })
  }

  /** The partial products [[signedPartials]] returns, each sign extended to
   *  `outW` bits at its own weight.
   */
  def placedSigned(parts: Seq[UInt], sliceW: Int, outW: Int): Seq[UInt] = {
    parts.zipWithIndex.map({ case (p, i) =>
      fit(p.asSInt.pad(outW).asUInt << (i * sliceW), outW)
    })
  }

  /** The low `w` bits of `x`, zero extended if it is narrower. */
  def fit(x: UInt, w: Int): UInt = {
    x.pad(w)(w - 1, 0)
  }

  /** The low `w` bits of the sum of `ops`. Three operands at a time become two
   *  by a carry-save step, so however many there are the sum ends in one carry
   *  chain.
   */
  def sumWide(ops: Seq[UInt], w: Int): UInt = {
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

/** Product of two integers, one partial product per `sliceW` bit slice.
 *
 * `io.out.bits` is the product of `io.in.bits.x` and `io.in.bits.y`, `xW + yW`
 * bits wide, unsigned by default and two's complement when `signed` is set
 * (the operands and the result are then read with `asSInt`).
 *
 * The wider operand is cut into slices of `sliceW` bits and each of them
 * multiplies the narrower one, so that every multiplication is one hard
 * multiplier of `sliceW` bits on its narrow port; the partial products are
 * summed afterwards.
 *
 * A two's complement product has two forms. By default the operands are read
 * unsigned in the slices and what that costs is subtracted from the first
 * partial product, which is this module's business and not its caller's;
 * where the product fits a hard multiplier, its post adder absorbs the
 * subtraction. With `small` the slices of the wider operand are read signed,
 * each in its own register at `"in"`; that uses less logic and suits a wider
 * product, and maps onto a hard multiplier only when `"in"` holds a register.
 * An unsigned product has one form and ignores `small`.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"pp"`     after the partial products, which is where a hard multiplier
 *               keeps its own register
 *  - `"chainK"` after the k-th partial sum, `K` = 1 ... [[chainPoints]], so
 *               only when an unsigned product is cut into more than two slices
 *               and its partial products are summed in a chain
 *  - `"out"`    the module output, after the sum of the partial products
 *
 * [[getStage]] gives the resulting latency, and `io.out.valid` is
 * `io.in.valid` delayed by that many cycles. The module is inlined into the
 * one that instantiates it, so a product is the same netlist wherever it is
 * written and a hierarchy boundary never sits between it and its neighbours.
 *
 * {{{
 * // one partial product per 17 bit slice of x, registered, then summed
 * new SplitMultiplier(24, 24, 17,
 *   PipelineStageConfig.specify("pp" -> 1, "out" -> 1))
 * }}}
 *
 * @constructor create a chisel Module SplitMultiplier
 * @param xW     width of the left hand operand
 * @param yW     width of the right hand operand
 * @param sliceW width of an operand slice
 * @param stage  pipeline registers of this module
 * @param signed read both operands and the result as two's complement
 * @param small  read the slices of a two's complement operand signed
 *
 * @see [[rial.arith.MulImpl]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class SplitMultiplier(
  xW:     Int,
  yW:     Int,
  sliceW: Int,
  stage:  PipelineStageConfig = PipelineStageConfig.none,
  signed: Boolean = false,
  small:  Boolean = false,
) extends Module with InlineInstance {

  /** Number of boundaries the chain of partial sums offers here. */
  val nChain = SplitMultiplier.chainPoints(xW, yW, sliceW, signed)

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = SplitMultiplier.stagePoints(nChain)

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  /** Width of the product. */
  val outW = xW + yW

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new SplitMultiplier.Input(xW, yW)))
    val out = Valid(UInt(outW.W))
  })

  // --------------------------- stage "in" ---------------------------
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  // --------------------------- stage "pp" ---------------------------
  private val cutPp = new SharedPipe("pp", cutIn.valid, stages.at("pp"))

  /** The cut the sum of the partial products ends at. */
  private var cutLast: SharedPipe = cutPp

  /** The sum of the partial products [[Slice.partials]] returns, each adding
   *  its own to what came before it, shifted down by the slice width; the
   *  boundary `"chainK"` sits between two of those additions.
   */
  private def chained(parts: Seq[UInt]): UInt = {
    if (parts.length == 1) {
      Slice.fit(parts.head, outW)
    } else {
      var carried = parts.head
      var rest    = parts.tail
      var low     = Seq.empty[UInt]
      var k       = 1
      while (rest.nonEmpty) {
        val next = rest.head +& (carried >> sliceW)
        val done = low :+ carried(sliceW - 1, 0)
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
      Slice.fit(Cat(carried, Cat(low.reverse)), outW)
    }
  }

  val prod = if (!signed) {
    // one partial product per slice, each added to what came before it
    val xIn = cutIn(io.in.bits.x)
    val yIn = cutIn(io.in.bits.y)
    chained(Slice.partials(xIn, yIn, sliceW).map(p => cutPp(p)))
  } else if (!small) {
    // the partial products and the cost of their sign in one sum: with `xu`
    // and `yu` the operands read unsigned, `x = xu - x[xW-1] * 2^xW` and
    // `y = yu - y[yW-1] * 2^yW` give
    //
    //   x * y = xu*yu - (x < 0 ? yu : 0) * 2^xW - (y < 0 ? xu : 0) * 2^yW
    //
    // modulo 2^(xW+yW), and both terms beside the product are selections of
    // an operand, so their sum is formed while the partial products are
    val xIn   = cutIn(io.in.bits.x)
    val yIn   = cutIn(io.in.bits.y)
    val parts = Slice.partials(xIn, yIn, sliceW)
    val m   = xW.min(yW)
    val sub = (Mux(xIn(xW-1), yIn, 0.U) << (xW - m)) +&
              (Mux(yIn(yW-1), xIn, 0.U) << (yW - m))
    val head = Slice.fit(parts.head, outW) - Cat(
      Slice.fit(sub, outW - m), 0.U(m.W))
    Slice.sumWide(
      (head +: Slice.placed(parts.tail, sliceW, outW, 1)).map(t => cutPp(t)),
      outW)
  } else {
    // signed slices and signed products, so the sum has no term but them
    val parts = Slice.signedPartials(io.in.bits.x, io.in.bits.y, sliceW,
      v => cutIn(v))
    Slice.sumWide(
      Slice.placedSigned(parts.map(p => cutPp(p)), sliceW, outW), outW)
  }

  require(prod.getWidth == outW,
    s"product is ${prod.getWidth} bits, expected $outW")

  // --------------------------- stage "out" --------------------------
  io.out := Pipe(cutLast.valid, prod, stages.at("out"))
}

/** Input bundle and pipeline stage points of [[rial.arith.SplitMultiplier]].
 *
 * @see [[rial.util.StagePoint]]
 */
object SplitMultiplier {

  /** Input of [[rial.arith.SplitMultiplier]]: the two factors.
   *
   * @param xW width of the left hand operand
   * @param yW width of the right hand operand
   */
  class Input(xW: Int, yW: Int) extends Bundle {
    val x = UInt(xW.W)
    val y = UInt(yW.W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.arith.SplitMultiplier]]. The `nChain` boundaries `"chain1"` ...
   *  `"chain<nChain>"` exist only when the sum of the partial products is a
   *  chain of more than one addition; [[chainPoints]] is how many that is.
   *
   * @param nChain number of boundaries inside the sum
   */
  def stagePoints(nChain: Int = 0): Seq[StagePoint] = {
    (StagePoint("pp", 3.0) +:
      (1 to nChain).map(k => StagePoint(s"chain$k", 2.0))) :+
      StagePoint("out", 3.0)
  }

  /** Number of boundaries the sum of the partial products offers for these
   *  parameters, which is one less than the number of additions it takes.
   *
   * A product that is not cut, or cut in two, is one addition and has none; a
   * two's complement product is summed all at once rather than in a chain,
   * in either of its forms, and has none either.
   *
   * @param xW     width of the left hand operand
   * @param yW     width of the right hand operand
   * @param sliceW width of an operand slice
   * @param signed the operands are two's complement
   */
  def chainPoints(xW: Int, yW: Int, sliceW: Int,
    signed: Boolean = false): Int = {

    if (signed) {
      0
    } else {
      (Slice.bits(xW.max(yW), sliceW).length - 2).max(0)
    }
  }
}

/** A multiply-add: the product of two integers, shifted, with an addend.
 *
 * `io.out.bits` is
 *
 *   `io.out.bits = (x * y >> shift) + c`
 *
 * with `x`, `y` and `c` the fields of `io.in.bits`, unsigned by default and
 * two's complement when `signed` is set (the operands, the addend and the
 * result are then read with `asSInt`). The product is divided by `2^shift`,
 * rounded down, and the addend is at the weight of that, so a caller hands
 * over the number it has. [[MultiplyAdd.outWidth]] is the width of the result.
 *
 * The module is one multiplication where the narrower operand fits a slice of
 * `sliceW` bits, and one per slice of the wider operand where it does not;
 * the addend is one more term of the sum that adds them up. The partial
 * products and the addend are what the boundary `"pp"` holds and the sum is
 * behind it, so a register at `"pp"` leaves each multiplication and its addend
 * on the near side of it, which is where a hard multiplier keeps its own.
 *
 * A two's complement multiply-add has the two forms of
 * [[rial.arith.SplitMultiplier]]: by default the operands are read unsigned
 * and what that costs is subtracted from the first partial product; with
 * `small` the slices of the wider operand are read signed, each in its own
 * register at `"in"`, which uses less logic and suits a wider product, and
 * maps onto a hard multiplier only when `"in"` holds a register.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`  the module input
 *  - `"pp"`  after the partial products, which the addend crosses beside them
 *  - `"out"` the module output, after the sum and the shift
 *
 * [[getStage]] gives the resulting latency, and `io.out.valid` is
 * `io.in.valid` delayed by that many cycles. The module is inlined into the
 * one that instantiates it, so a multiply-add is the same netlist wherever it
 * is written and a hierarchy boundary never sits between it and its
 * neighbours.
 *
 * {{{
 * // a two's complement multiply-add on 17 bit slices, its product divided by
 * // 2^18 and its addend of 25 bits at the weight of that, its partial
 * // products and its addend registered
 * new MultiplyAdd(19, 24, 25, 18, 17,
 *   PipelineStageConfig.specify("pp" -> 1), signed = true)
 * }}}
 *
 * @constructor create a chisel Module MultiplyAdd
 * @param xW     width of the left hand operand
 * @param yW     width of the right hand operand
 * @param cW     width of the addend, which is at the weight of the result
 * @param shift  bits the product drops at its bottom
 * @param sliceW width of an operand slice
 * @param stage  pipeline registers of this module
 * @param signed read the operands, the addend and the result as two's
 *               complement
 * @param small  read the slices of a two's complement operand signed
 *
 * @see [[rial.arith.MulImpl]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class MultiplyAdd(
  xW:     Int,
  yW:     Int,
  cW:     Int,
  shift:  Int,
  sliceW: Int,
  stage:  PipelineStageConfig = PipelineStageConfig.none,
  signed: Boolean = false,
  small:  Boolean = false,
) extends Module with InlineInstance {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = MultiplyAdd.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  /** True when the narrower operand fits a slice, so that the multiply-add is
   *  one multiplication and the wider operand is not cut.
   */
  val oneMultiply = xW.min(yW) <= sliceW

  val outW = MultiplyAdd.outWidth(xW, yW, cW, shift)

  /** Width at which the terms are summed: one bit above both the product and
   *  the addend, so that the sum is exact before the shift.
   */
  private val sumW = outW + shift

  require(shift >= 0,
    s"a multiply-add drops no negative number of bits, but shift is $shift")
  require(outW > 0,
    s"a multiply-add of $xW by $yW bits with an addend of $cW bits has " +
    s"nothing left after a shift of $shift")

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new MultiplyAdd.Input(xW, yW, cW)))
    val out = Valid(UInt(outW.W))
  })

  // --------------------------- stage "in" ---------------------------
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val cIn = cutIn(io.in.bits.c)

  // --------------------------- stage "pp" ---------------------------
  private val cutPp = new SharedPipe("pp", cutIn.valid, stages.at("pp"))

  /** `v` divided by `2^shift`, rounded down the way the module reads it. */
  private def dropped(v: UInt): UInt = {
    if (shift == 0) {
      v
    } else if (signed) {
      (v.asSInt >> shift).asUInt
    } else {
      v >> shift
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
  private def signCost(xIn: UInt, yIn: UInt): UInt = {
    val m = xW.min(yW)
    val cost = (Mux(xIn(xW-1), yIn.asSInt, 0.S(yW.W)) << (xW - m)) +&
               (Mux(yIn(yW-1), xIn, 0.U).zext << (yW - m))
    Slice.fit((cost.pad(sumW) << m).asUInt, sumW)
  }

  /** The addend at the weight of the product, `sumW` bits wide, which is where
   *  it joins the product: `(x * y + c * 2^shift) >> shift` is the same number
   *  as `(x * y >> shift) + c`, and the addend reaches the sum of the partial
   *  products without an adder of its own.
   */
  private def atProduct(c: UInt): UInt = {
    if (signed) {
      Slice.fit((c.asSInt.pad(sumW) << shift).asUInt, sumW)
    } else {
      Slice.fit(c << shift, sumW)
    }
  }

  // the products and what they need besides, each registered at "pp" and at
  // its weight in a sum that is exact at `sumW` bits; the addend is one more
  // term of it
  private val terms = if (signed && small) {
    // signed slices and signed products, so nothing is added for the sign
    val parts = if (oneMultiply) {
      Seq((cutIn(io.in.bits.x).asSInt * cutIn(io.in.bits.y).asSInt).asUInt)
    } else {
      Slice.signedPartials(io.in.bits.x, io.in.bits.y, sliceW, v => cutIn(v))
    }
    Slice.placedSigned(parts.map(p => cutPp(p)), sliceW, sumW)
  } else {
    // the operands are read unsigned, so that what a multiplication is given
    // is a slice, and what it costs to read a two's complement operand that
    // way is subtracted from the first partial product
    val xIn = cutIn(io.in.bits.x)
    val yIn = cutIn(io.in.bits.y)
    val parts = if (oneMultiply) {
      Seq(xIn * yIn)
    } else {
      Slice.partials(xIn, yIn, sliceW)
    }
    val atWeight = Slice.fit(parts.head, sumW) +:
      Slice.placed(parts.tail, sliceW, sumW, 1)
    (if (signed) {
      (atWeight.head - signCost(xIn, yIn)) +: atWeight.tail
    } else {
      atWeight
    }).map(t => cutPp(t))
  }

  val prod = dropped(Slice.sumWide(terms :+ cutPp(atProduct(cIn)), sumW))

  require(prod.getWidth == outW,
    s"multiply-add is ${prod.getWidth} bits, expected $outW")

  // --------------------------- stage "out" --------------------------
  io.out := Pipe(cutPp.valid, prod, stages.at("out"))
}

/** Input bundle and pipeline stage points of [[rial.arith.MultiplyAdd]].
 *
 * @see [[rial.util.StagePoint]]
 */
object MultiplyAdd {

  /** Input of [[rial.arith.MultiplyAdd]]: the two factors and the addend.
   *
   * @param xW width of the left hand operand
   * @param yW width of the right hand operand
   * @param cW width of the addend
   */
  class Input(xW: Int, yW: Int, cW: Int) extends Bundle {
    val x = UInt(xW.W)
    val y = UInt(yW.W)
    val c = UInt(cW.W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.arith.MultiplyAdd]].
   */
  def stagePoints: Seq[StagePoint] = {
    Seq(StagePoint("pp", 3.0), StagePoint("out", 3.0))
  }

  /** Width of the result of [[rial.arith.MultiplyAdd]] for these parameters,
   *  which is one bit above both the shifted product and the addend.
   *
   * @param xW    width of the left hand operand
   * @param yW    width of the right hand operand
   * @param cW    width of the addend
   * @param shift bits the product drops at its bottom
   */
  def outWidth(xW: Int, yW: Int, cW: Int, shift: Int): Int = {
    (xW + yW - shift).max(cW) + 1
  }
}
