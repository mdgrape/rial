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
 *
 * @param kind   `"plain"` or `"split"`
 * @param sliceW width of an operand slice, for `"split"`
 *
 * @see [[rial.arith.Multiplier]]
 */
case class MulImpl(kind: String, sliceW: Option[Int] = None) {

  require(MulImpl.kinds.contains(kind),
    s"unknown multiplier structure '$kind'; " +
    s"known ones are ${MulImpl.kinds.mkString(", ")}")
  require(!isSplit || sliceW.exists(_ > 0),
    "a split multiplier needs a positive slice width")

  /** True for the structure that slices one operand. */
  def isSplit: Boolean = {
    kind == "split"
  }

  /** Short name, for a generated module or file name, e.g. `"split17"`. */
  def tag: String = {
    kind + sliceW.map(_.toString).getOrElse("")
  }
}

/** Constructor shorthands for [[rial.arith.MulImpl]]. */
object MulImpl {

  /** Every structure [[rial.arith.MulImpl]] names. */
  val kinds: Seq[String] = Seq("plain", "split")

  /** One product. */
  val plain: MulImpl = MulImpl("plain")

  /** One partial product per `w` bit slice of the wider operand. */
  def split(w: Int): MulImpl = {
    MulImpl("split", Some(w))
  }

  /** Reads `"plain"` or `"split:<w>"`.
   *
   * @throws IllegalArgumentException if `s` names neither.
   */
  def parse(s: String): MulImpl = {
    val t = s.trim
    if (t == "plain") {
      plain
    } else if (t.startsWith("split:")) {
      val v = t.drop("split:".length)
      val w = try {
        v.toInt
      } catch {
        case _: NumberFormatException => {
          throw new IllegalArgumentException(
            s"multiplier structure 'split:$v': '$v' is not an int")
        }
      }
      split(w)
    } else {
      throw new IllegalArgumentException(
        s"unknown multiplier structure '$s'; expected plain or split:<w>")
    }
  }
}

/** Product of two integers, with the structure [[rial.arith.MulImpl]] names.
 *
 * `io.out.bits` is the `xW + yW` bit product of `io.in.bits.x` and
 * `io.in.bits.y`, unsigned by default and two's complement when `signed` is
 * set (the operands and the result are then read with `asSInt`). A module that
 * needs a product instantiates this one and does nothing else about it: what
 * a wide operand costs, and what a two's complement one costs once it is
 * sliced, is this module's business. `MulImpl.plain` is the plain product and
 * asks nothing of the operands.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`  the module input
 *  - `"pp"`  after the partial products, only with `MulImpl.split`
 *  - `"out"` the module output, after the sum of the partial products
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
 * }}}
 *
 * @constructor create a chisel Module Multiplier
 * @param xW     width of the left hand operand
 * @param yW     width of the right hand operand
 * @param impl   structure of the product
 * @param stage  pipeline registers of this module
 * @param signed read both operands and the result as two's complement
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
) extends Module with InlineInstance {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = Multiplier.stagePoints(impl)

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val outW = xW + yW

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new Multiplier.Input(xW, yW)))
    val out = Valid(UInt(outW.W))
  })

  // --------------------------- stage "in" ---------------------------
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val xIn = cutIn(io.in.bits.x)
  val yIn = cutIn(io.in.bits.y)

  // --------------------------- stage "pp" ---------------------------
  private val nPp = if (impl.isSplit) {
    stages.at("pp")
  } else {
    0
  }

  private val cutPp = new SharedPipe("pp", cutIn.valid, nPp)

  // `plain` is the product as it would be written by hand and nothing else: a
  // hard multiplier takes both operands as they are, so there is nothing to
  // shape. A slice of an operand is what does not fit one, and the slices are
  // unsigned, so a two's complement operand is read unsigned there and what
  // that costs comes back: `x = xu - x[xW-1] * 2^xW` gives
  //
  //   x * y = xu*yu - (x < 0 ? yu : 0) * 2^xW - (y < 0 ? xu : 0) * 2^yW
  //
  // modulo 2^(xW+yW). The two terms beside the product are selections of an
  // operand, so their sum is formed while the partial products are and then
  // rides on the addend input of the least significant multiplier, which is
  // where a hard multiplier subtracts a term for free.
  val prod = if (!impl.isSplit) {
    if (signed) {
      (xIn.asSInt * yIn.asSInt).asUInt
    } else {
      xIn * yIn
    }
  } else {
    val raw = Multiplier.partials(xIn, yIn, impl)
    if (signed) {
      val m   = xW.min(yW)
      val sub = (Mux(xIn(xW-1), yIn, 0.U) << (xW - m)) +&
                (Mux(yIn(yW-1), xIn, 0.U) << (yW - m))
      val head  = Multiplier.fit(raw.head, outW) - Cat(
        Multiplier.fit(sub, outW - m), 0.U(m.W))
      val terms = head +: Multiplier.placed(raw.tail, impl, outW, 1)
      Multiplier.sumWide(terms.map(t => cutPp(t)), outW)
    } else {
      Multiplier.sum(raw.map(p => cutPp(p)), impl, outW)
    }
  }

  // --------------------------- stage "out" --------------------------
  io.out := Pipe(cutPp.valid, prod, stages.at("out"))
}

/** Input bundle and pipeline stage points of [[rial.arith.Multiplier]].
 *
 * @see [[rial.util.StagePoint]]
 */
object Multiplier {

  /** Input of [[rial.arith.Multiplier]]: the two factors.
   *
   * @param xW width of the left hand operand
   * @param yW width of the right hand operand
   */
  class Input(xW: Int, yW: Int) extends Bundle {
    val x = UInt(xW.W)
    val y = UInt(yW.W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.arith.Multiplier]]. `"pp"` exists only with `MulImpl.split`.
   */
  def stagePoints(impl: MulImpl): Seq[StagePoint] = {
    if (impl.isSplit) {
      Seq(StagePoint("pp", 3.0), StagePoint("out", 3.0))
    } else {
      Seq(StagePoint("out", 4.0))
    }
  }

  /** The slices `impl` cuts an operand of `w` bits into, least significant
   *  first, as (high bit, low bit) pairs.
   */
  private def sliceBits(w: Int, impl: MulImpl): Seq[(Int, Int)] = {
    impl.sliceW match {
      case Some(sw) if impl.isSplit && sw < w => {
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

  /** The `outW` bit product of the partial products [[partials]] returns.
   *
   * The low slice of the least significant partial product is already a slice
   * of the product, and the rest is one flat sum, so the whole summation is
   * one carry chain.
   */
  private def sum(parts: Seq[UInt], impl: MulImpl, outW: Int): UInt = {
    if (parts.length == 1) {
      fit(parts.head, outW)
    } else {
      val sw   = impl.sliceW.get
      val hiW  = outW - sw
      val ends = parts.zipWithIndex.map({ case (p, i) =>
        if (i == 0) {
          fit(p >> sw, hiW)
        } else if (i == 1) {
          fit(p, hiW)
        } else {
          Cat(fit(p, hiW - (i - 1) * sw), 0.U(((i - 1) * sw).W))
        }
      })
      Cat(sumWide(ends, hiW), parts.head(sw - 1, 0))
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
