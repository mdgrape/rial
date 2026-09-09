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

/** Structure of a product, chosen by the user of the module that builds it.
 *
 * `kind` names the structure and `sliceW` carries its parameter, so that a
 * configuration read from a file can be turned into one of these without
 * knowing the set of structures:
 *
 *  - `MulImpl.plain`    one product.
 *  - `MulImpl.split(w)` one partial product per `w` bit slice of the wider
 *                       operand, summed afterwards. A pipeline register may be
 *                       placed on the partial products, which is where a hard
 *                       multiplier of at most `w` bits on its narrow port has
 *                       one.
 *  - `MulImpl.macc(w)`  a multiply-add: the product carries an addend, and it
 *                       is shaped so that the addend rides on the addend input
 *                       a hard multiplier of at most `w` bits on its narrow
 *                       port has, rather than on an adder of its own.
 *                       [[rial.arith.Multiplier]] asks its caller for the
 *                       addend, and only with this structure; a module that
 *                       has none of its own gives zero.
 *
 * @param kind   `"plain"`, `"split"` or `"macc"`
 * @param sliceW width of an operand slice, for `"split"` and `"macc"`
 *
 * @see [[rial.arith.Multiplier]]
 */
case class MulImpl(kind: String, sliceW: Option[Int] = None) {

  require(MulImpl.kinds.contains(kind),
    s"unknown multiplier structure '$kind'; " +
    s"known ones are ${MulImpl.kinds.mkString(", ")}")
  require(!isSliced || sliceW.exists(_ > 0),
    s"a $kind multiplier needs a positive slice width")

  /** True for the structure that slices one operand. */
  def isSplit: Boolean = {
    kind == "split"
  }

  /** True for the structure that carries an addend. */
  def isMacc: Boolean = {
    kind == "macc"
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
  val kinds: Seq[String] = Seq("plain", "split", "macc")

  /** One product. */
  val plain: MulImpl = MulImpl("plain")

  /** One partial product per `w` bit slice of the wider operand. */
  def split(w: Int): MulImpl = {
    MulImpl("split", Some(w))
  }

  /** A multiply-add whose addend rides on a hard multiplier of `w` bits. */
  def macc(w: Int): MulImpl = {
    MulImpl("macc", Some(w))
  }

  /** Reads `"plain"`, `"split:<w>"` or `"macc:<w>"`.
   *
   * @throws IllegalArgumentException if `s` names none of them.
   */
  def parse(s: String): MulImpl = {
    val t = s.trim
    val sliced = Seq("split" -> (split _), "macc" -> (macc _)).collectFirst {
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
        s"expected plain, split:<w> or macc:<w>"))
    }
  }
}

/** Product of two integers, with the structure [[rial.arith.MulImpl]] names.
 *
 * `io.out.bits` is the product of `io.in.bits.x` and `io.in.bits.y`, unsigned
 * by default and two's complement when `signed` is set (the operands, the
 * addend and the result are then read with `asSInt`).
 *
 * `MulImpl.macc` makes the product a multiply-add: the input then carries an
 * addend `io.in.bits.c` at the weight of the product, and `cW` gives its
 * width. The other structures have no addend and must not be given one.
 *
 * With `shift` the result is divided by `2^shift`, rounded down, so that
 * `io.out.bits` is `(x * y + c) >> shift`; every structure takes it, and a
 * module that wants no shift leaves it out.
 *
 * A module that needs a product instantiates this one and does nothing else
 * about it: what a wide operand costs, what a two's complement one costs once
 * it is sliced, and where the addend has to sit for the product to stay one
 * multiply-add, is this module's business. `MulImpl.plain` is the plain
 * product and asks nothing of the operands.
 *
 * A two's complement multiply-add splits the right hand operand in two when
 * neither operand fits a slice, so it then needs `yW` to be at most twice the
 * slice width; the left hand operand is not sliced. Where the narrower
 * operand already fits a slice the product is one multiplication and the
 * addend rides on it. Both shapes put the multiply and the addend of a tile
 * behind the same boundary, so a register at `"pp"` and none at `"out"`
 * leaves each tile a multiply-add of one hard multiplier.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"pp"`     after the partial products, only with `MulImpl.split` and
 *               `MulImpl.macc`; with `macc` these are the tiles of the
 *               multiply-add
 *  - `"chainK"` after the k-th partial sum, `K` = 1 ... one less than the
 *               number of hops, so only when a plain product is split into
 *               more than two slices
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
 * // a two's complement multiply-add, its result divided by 2^18, its tiles
 * // registered
 * new Multiplier(19, 24, MulImpl.macc(17),
 *   PipelineStageConfig.specify("pp" -> 1), true, Some(41), Some(18))
 * }}}
 *
 * @constructor create a chisel Module Multiplier
 * @param xW     width of the left hand operand
 * @param yW     width of the right hand operand
 * @param impl   structure of the product
 * @param stage  pipeline registers of this module
 * @param signed read both operands and the result as two's complement
 * @param cW     width of the addend, which `MulImpl.macc` has and the other
 *               structures do not
 * @param shift  bits the result drops at its bottom, when it drops any
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

  /** Bit at which the right hand operand is split for a two's complement
   *  multiply-add, or zero when the product is one multiplication.
   */
  val loW = Multiplier.loSliceW(xW, yW, impl, signed, shiftW)

  val outW = Multiplier.outWidth(xW, yW, cW, shiftW)

  /** Width of the product on its own. */
  private val prodW = xW + yW

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

  /** The cuts a partial product crosses before it reaches the sum, so that an
   *  addend that joins the sum can be made to cross the same ones.
   */
  private var cutsToSum: Seq[SharedPipe] = Seq(cutPp)

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
          cutLast   = cut
          cutsToSum = cutsToSum :+ cut
          k        += 1
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

  /** A two's complement multiply-add as two multiplications that each take one
   *  addend, so that each of them is one multiply-add of a hard multiplier.
   *
   * The left hand operand enters in offset binary, `x = xOff - 2^(xW-1)`, and
   * the right hand one unsigned and split at `loW` into `yHi` and `yLo`, so that
   * with `off = c - y * 2^(xW-1)` gathering the offset of `xOff` and the
   * addend, and `off = offHi * 2^loW + offLo`,
   *
   *   x * y + c = 2^loW * (xOff*yHi + offHi) - neg + 2^loW - 1
   *
   * where `neg = (y < 0 ? xOff : 0) * 2^yW + (2^loW - 1 - offLo) - xOff*yLo`.
   * The `2^loW - 1` is what makes the shift of a subtracted term exact, and
   * `2^loW - 1 - offLo` is the bitwise complement of `offLo`, so the two
   * corrections stay one addend each.
   */
  private def offsetTiles(c: UInt): UInt = {
    val xOff = Cat(~xIn(xW-1), xIn(xW-2, 0))
    val off  = c.asSInt -& (yIn.asSInt << (xW - 1))

    val prodHi = (xOff * yIn(yW-1, loW)).zext +& (off >> loW)
    val negBase = Cat(Mux(yIn(yW-1), xOff, 0.U(xW.W)),
      0.U((yW - loW).W), ~off(loW-1, 0))
    val prodNeg = negBase.zext -& (xOff * yIn(loW-1, 0)).zext

    require(prodHi.getWidth == Multiplier.hiWidth(xW, yW, cW.get, loW),
      s"multiply-add tile is ${prodHi.getWidth} bits, " +
      s"expected ${Multiplier.hiWidth(xW, yW, cW.get, loW)}")
    require(prodNeg.getWidth == Multiplier.negWidth(xW, yW, loW),
      s"multiply-add tile is ${prodNeg.getWidth} bits, " +
      s"expected ${Multiplier.negWidth(xW, yW, loW)}")

    val hiPp  = cutPp(prodHi)
    val negPp = cutPp(prodNeg)

    // the low `loW` bits of the sum are dropped anyway when the shift reaches
    // past them, and then the two tiles meet over the high bits alone
    val sum = if (loW <= shiftW) {
      (hiPp -& (negPp >> loW)) >> (shiftW - loW)
    } else {
      (Cat(hiPp.asUInt, Fill(loW, 1.U(1.W))).asSInt -& negPp) >> shiftW
    }

    require(sum.getWidth >= outW,
      s"multiply-add result is ${sum.getWidth} bits, expected at least $outW")

    sum.asUInt(outW-1, 0)
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

  /** The product of the two operands with the addend `c` on it, shaped so that
   *  each multiplication and the addend it carries stay one multiply-add.
   */
  private def multiplyAdd(c: UInt): UInt = {
    if (signed && loW == 0) {
      // one multiplication and its addend: the operands fit one slice, so
      // there is nothing to slice and nothing the addend has to be shaped
      // around
      dropped(cutPp(
        ((xIn.asSInt * yIn.asSInt) +& c.asSInt).asUInt))
    } else if (!signed) {
      // one partial product per slice chained into the sum, the addend on the
      // final sum across the cuts the partial products cross: an unsigned
      // multiply-add
      dropped(chained(Multiplier.partials(xIn, yIn, impl).map(p => cutPp(p)))
        +& cutsToSum.foldLeft(c)((d, cut) => cut(d)))
    } else {
      // two tiles in offset binary, one multiply-add each: a two's complement
      // multiply-add whose operands do not fit one slice
      offsetTiles(c)
    }
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

  /** Width of the result of [[rial.arith.Multiplier]] for these parameters.
   *
   * @param xW    width of the left hand operand
   * @param yW    width of the right hand operand
   * @param cW    width of the addend, when there is one
   * @param shift bits the result drops at its bottom
   */
  def outWidth(xW: Int, yW: Int, cW: Option[Int], shift: Int): Int = {
    cW match {
      case Some(w) => {
        (xW + yW).max(w) + 1 - shift
      }
      case None => {
        xW + yW - shift
      }
    }
  }

  /** Bit at which a two's complement multiply-add splits its right hand operand,
   *  or zero when it forms the product as one multiplication.
   *
   * There is no split when the product already fits one multiplication, that is
   * when the narrower operand is within a slice. Otherwise the split leaves both
   * halves within one slice and, where it can, at or below the shift, so that
   * the two halves meet again over the high bits alone.
   *
   * @param xW     width of the left hand operand
   * @param yW     width of the right hand operand
   * @param impl   structure of the product
   * @param signed the operands are two's complement
   * @param shift  bits the result drops at its bottom
   *
   * @throws IllegalArgumentException if the right hand operand is wider than
   *         two slices, which this form does not split.
   */
  def loSliceW(xW: Int, yW: Int, impl: MulImpl, signed: Boolean,
    shift: Int): Int = {

    impl.sliceW match {
      case Some(w) if impl.isMacc && signed && xW.min(yW) > w => {
        require(yW <= 2 * w,
          s"a two's complement multiply-add splits its right hand operand in " +
          s"two, so $yW bits do not fit slices of $w bits")
        val lower = (yW - w).max(1)
        val upper = w.min(shift)
        if (upper >= lower) {
          (yW / 2).max(lower).min(upper)
        } else {
          (yW / 2).max(lower).min(w)
        }
      }
      case _ => {
        0
      }
    }
  }

  /** Width of the tile of a two's complement multiply-add that carries the
   *  addend; see [[rial.arith.Multiplier]].
   *
   * @param xW  width of the left hand operand
   * @param yW  width of the right hand operand
   * @param cW  width of the addend
   * @param loW bit at which the right hand operand is split
   */
  def hiWidth(xW: Int, yW: Int, cW: Int, loW: Int): Int = {
    val mulW = xW + yW - loW + 1
    val addW = cW.max(xW + yW - 1) + 1 - loW
    mulW.max(addW) + 1
  }

  /** Width of the tile of a two's complement multiply-add that carries the
   *  offsets; see [[rial.arith.Multiplier]].
   *
   * @param xW  width of the left hand operand
   * @param yW  width of the right hand operand
   * @param loW bit at which the right hand operand is split
   */
  def negWidth(xW: Int, yW: Int, loW: Int): Int = {
    (xW + yW + 1).max(xW + loW + 1) + 1
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.arith.Multiplier]]. `"pp"` exists only with `MulImpl.split`, and
   *  the `nChain` boundaries `"chain1"` ... `"chain<nChain>"` only when the
   *  sum of the partial products is a chain of more than one addition;
   *  [[chainPoints]] is how many that is.
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
   * none; a two's complement product is summed together with the terms its
   * sign costs rather than in a chain, and has none either.
   *
   * @param xW     width of the left hand operand
   * @param yW     width of the right hand operand
   * @param signed the operands are two's complement
   */
  def chainPoints(impl: MulImpl, xW: Int, yW: Int,
    signed: Boolean = false): Int = {

    if (signed) {
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
