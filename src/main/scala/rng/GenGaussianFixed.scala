//% @file GenGaussianFixed.scala
//
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.rng

import scala.math._

import chisel3._
import chisel3.util._

import rial.arith._
import rial.fixed._
import rial.table._
import rial.util._

import org.apache.commons.math3.special.Erf

//
// |z| of the gaussian, in units of 2^zIntW, the power of two just above the
// largest |z| of the table (at most 2^(intW-1)).
//
// Rows, each split into 2^splitW sub-intervals of x:
//   - row 0, the central page: |z| = sqrt(2) * erfInv(x), x = t/2 in [0, 0.5).
//   - row 1+i, octave i of the tail: |z| = sqrt(2) * erfcInv(x), where the
//     position t in [0, 1) runs from x = 2^-(i+1) down to 2^-(i+2), so that
//     |z| grows with t in every row.
//   - the last row: x from 2^-(nOctaves+1) down to 2^-(nOctaves+1+lumpOctaves),
//     linearly in t, with |z| capped at its value at the deep end.
// No row goes beyond zCap, the largest magnitude the output holds.
//
private[rial] class ProbitFixedTable(
  val zSpec:     FixedSpec,
  val nOrder:    Int,
  val nOctaves:  Int,
  val posW:      Int,
  val splitW:    Int,
  val guardBits: Int,
) {
  private val sqrt2 = sqrt(2.0)

  val zCap = GenRandomGaussianFixed.zCap(zSpec)

  val lastBegin = pow(2.0, -(nOctaves+1))
  val lastEnd   = pow(2.0, -(nOctaves+1+GenRandomGaussianFixed.lumpOctaves))
  val lastCap   = min(ProbitTable.probitHalf(lastEnd), zCap)

  // the evaluator needs a table that reaches the upper half of its output, so
  // the output covers no more integer bits than the tail reaches
  val zIntW  = (0 until zSpec.intW - 1).find(m => pow(2.0, m) > lastCap).getOrElse(zSpec.intW - 1)
  val zScale = pow(2.0, zIntW)

  // the top splitW bits of the position choose the sub-interval, the rest is
  // the polynomial input; the output is |z| in units of 2^-(fracW+guardBits)
  val inBitW  = posW - splitW
  val outBitW = zIntW + zSpec.fracW + guardBits

  val ft = new FuncTable(nOrder)

  val nSub = 1 << splitW

  // sub-interval 0 holds the xBegin end, the one the position 0 reads
  private def addSplit(f: Double => Double, xBegin: Double, xEnd: Double): Unit = {
    for(s <- 0 until nSub) {
      ft.addInterval(f, xBegin + (xEnd - xBegin) * s / nSub,
                        xBegin + (xEnd - xBegin) * (s+1) / nSub)
    }
  }

  addSplit((x: Double) => { min(sqrt2 * Erf.erfInv(x), zCap) / zScale }, 0.0, 0.5)

  // the octave that crosses zCap is fitted as it is, not capped: the kink of
  // the cap would widen the coefficients of every row, and the evaluator and
  // the rounding saturate the result at zCap anyway
  for(i <- 0 until nOctaves) {
    addSplit((x: Double) => { ProbitTable.probitHalf(x) / zScale },
      pow(2.0, -(i+1)), pow(2.0, -(i+2)))
  }

  addSplit((x: Double) => { min(ProbitTable.probitHalf(x), lastCap) / zScale },
    lastBegin, lastEnd)

  val nRows = 1 + nOctaves + 1

  assert(ft.intervals.length == nRows * nSub)

  /** The function row `row` (of the `nRows`, not of the sub-intervals)
   *  approximates, as |z|, at x. */
  def exact(row: Int, x: Double): Double = {
    if (row == 0) {
      min(sqrt2 * Erf.erfInv(x), zCap)
    } else if (row <= nOctaves) {
      min(ProbitTable.probitHalf(x), zCap)
    } else {
      min(ProbitTable.probitHalf(x), lastCap)
    }
  }

  // the probit is shallow at the top end of a row, where its polynomial may
  // reach slightly below zero; the evaluator clamps a negative result to zero
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW,
    negativeIsClamped = true)

  val bitW   = gen.getBitWidth()
  val coeffs = gen.generateTable(bitW)
}

/** Generates a gaussian random number in fixed point out of `rndW` uniform
 * random bits.
 *
 * The result `io.z` is a two's complement number of `zSpec`, symmetric about
 * zero, and held at `2^(intW-1) - 2^-fracW`, the largest magnitude `zSpec`
 * holds. It is the inverse CDF of the gaussian, approximated piecewise by
 * polynomials evaluated by [[rial.table.PolynomialEval]], rounded to the
 * nearest multiple of `2^-fracW`.
 *
 * The random word is read from the top down as the sign, one bit that picks
 * the central half or the tail, a scan field of `rndW - 2 - posW` bits and the
 * position of `posW` bits. Each bit of the scan field takes the tail one
 * octave of probability deeper, down to the depth at which |z| reaches the
 * largest magnitude of `zSpec`. [[GenRandomGaussianFixed.octaves]] gives the
 * number of bits that are read; the bits of the field above them are ignored.
 * So `rndW = 2 + posW + octaves` is the shortest word with the full tail, and
 * a shorter word gives a shallower one.
 *
 * The tail reaches |z| = [[zReach]], which [[GenRandomGaussianFixed.reach]]
 * computes without building the module, and no larger |z| occurs. With the
 * full tail, zReach is the largest magnitude of `zSpec`. With a shallower one
 * of `n` octaves, the words that set none of the `n` bits, with probability
 * `2^-(n+1)`, give |z| in the octave just below zReach in place of the tail
 * beyond it.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"decode"` after the row decode
 *  - `"table"`  after the coefficient table lookup
 *  - `"poly.*"` the boundaries of the [[rial.table.PolynomialEval]] submodule
 *               that evaluates the polynomial of the row, see its doc for the
 *               names; the group name `"poly"` stands for `"poly.out"`
 *  - `"round"`  after the magnitude is rounded to `zSpec`
 *  - `"out"`    the module output, after the sign is applied
 *
 * [[getStage]] gives the resulting latency, and `z.valid` is `rnd.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // Q3.13 (16 bits, 13 fractional), three registers placed by the module
 * new GenRandomGaussianFixed(32, FixedSpec(16, 13), stage = 3)
 * // hand-picked positions, one of them inside the polynomial
 * new GenRandomGaussianFixed(32, FixedSpec(16, 13),
 *   stage = PipelineStageConfig.specify("table" -> 1, "poly.prod2" -> 1, "out" -> 1))
 * // combinational
 * new GenRandomGaussianFixed(32, FixedSpec(16, 13), stage = PipelineStageConfig.none)
 * }}}
 *
 * @param rndW        input random number bit width.
 * @param zSpec       output fixed point spec; must be signed.
 * @param nOrder      order of the polynomial of a row.
 * @param stage       pipeline registers of this module.
 * @param splitW      every row of the table is divided into `2^splitW`
 *                    sub-intervals, each with its own polynomial.
 * @param guardBits   bits below `2^-fracW` that the polynomial is computed
 *                    with, rounded away at the end.
 * @param posW        width of the position inside a row; `zSpec.fracW + 4`
 *                    when `None`.
 * @param polyMulImpl structure of the products of the polynomial.
 *
 * @see [[rial.arith.MulImpl]]
 * @see [[rial.util.PipelineStageConfig]]
 */
class GenRandomGaussianFixed(
  rndW:        Int,
  zSpec:       FixedSpec,
  nOrder:      Int = 3,
  stage:       PipelineStageConfig,
  splitW:      Int = 3,
  guardBits:   Int = 3,
  posW:        Option[Int] = None,
  polyMulImpl: MulImpl = MulImpl.plain,
) extends Module {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = GenRandomGaussianFixed.stagePoints(nOrder, polyMulImpl)

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  require(zSpec.signed, "the result is a signed fixed point number")
  require(zSpec.intW >= 1, "the sign bit is part of intW")
  require(guardBits >= 0)

  val positionW = GenRandomGaussianFixed.positionWidth(zSpec, posW)
  val scanW     = rndW - 2 - positionW
  require(scanW >= 1, s"rndW=${rndW} leaves no scan field next to the ${positionW}-bit position")

  /** The number of octaves of the tail the scan field reads. */
  val nOctaves = GenRandomGaussianFixed.octaves(rndW, zSpec, posW)
  require(nOctaves >= 1, s"${zSpec} holds no magnitude beyond the central half")

  require(splitW >= 0 && positionW - splitW >= 2,
    "the position inside a sub-interval needs at least two bits")

  val table = GenRandomGaussianFixed.table(zSpec, nOrder, nOctaves, positionW,
    splitW, guardBits)
  val bitW   = table.bitW
  val coeffs = table.coeffs

  /** The largest |z| the module produces; see the class doc. */
  lazy val zReach = GenRandomGaussianFixed.reach(rndW, zSpec, nOrder, splitW,
    guardBits, posW)

  val io = IO(new Bundle {
    val rnd = Flipped(Valid(UInt(rndW.W)))
    val z   = Valid(UInt(zSpec.W.W))
  })

  val inBitW = table.inBitW

  val packed = VecInit(PolynomialGenerator.packTable(coeffs).map(c => {
    c.U(coeffs.bitW.coefBits.reduce(_+_).W)
  }))

  // ---------------------------------------------------------------
  // S0 : input registers ("in")

  private val cutIn = new SharedPipe("in", io.rnd.valid, stages.at("in"))

  val rndIn   = cutIn(io.rnd.bits)
  val validIn = cutIn.valid

  // ---------------------------------------------------------------
  // S1 : row decode

  val sgn  = rndIn(rndW-1)
  val page = rndIn(rndW-2)
  val scan = rndIn(positionW + nOctaves - 1, positionW)
  val pos  = rndIn(positionW-1, 0)

  // row 0 is the central page, row 1+i the octave whose bit i is the lowest
  // set bit of the scan, the last row the one with no bit set; the page sits
  // below the scan as one more, lowest, bit, so the row needs no selection of
  // its own
  val lastRow = (table.nRows - 1).U(log2Ceil(table.nRows).W)
  val row = GenRandomGaussian.firstSetSelect(Cat(scan, ~page), lastRow)(i => {
    i.U(log2Ceil(table.nRows).W)
  })

  // --------------------------- stage "decode" --------------------------
  private val cutDecode = new SharedPipe("decode", validIn, stages.at("decode"))

  val rowDecode   = cutDecode(row)
  val posDecode   = cutDecode(pos)
  val sgnDecode   = cutDecode(sgn)
  val validDecode = cutDecode.valid

  // the sub-interval the position falls in sits below the row
  val romAddr = if(splitW != 0) {
    Cat(rowDecode, posDecode(positionW-1, positionW-splitW))
  } else {
    rowDecode
  }

  // the coefficients are read as one word, the shape a block RAM can register
  val coefWord = packed(romAddr)

  // position inside the sub-interval, as a signed value in [-1, 1)
  val sdx = if(nOrder != 0) {
    val dx = posDecode(inBitW-1, 0)
    Some(Cat(~dx(inBitW-1), dx(inBitW-2, 0)).asSInt)
  } else {
    None
  }

  // --------------------------- stage "table" ---------------------------
  private val cutTable = new SharedPipe("table", validDecode, stages.at("table"))

  val coefTable  = cutTable(coefWord).asTypeOf(MixedVec(bitW.coefBits.map(w => SInt(w.W))))
  val sdxTable   = sdx.map(x => cutTable(x))
  val sgnTable   = cutTable(sgnDecode)
  val validTable = cutTable.valid

  // ---------------------------------------------------------------
  // S2 : polynomial of the row ("poly.*")

  val polynomial = Module(new PolynomialEval(bitW, stages.sub("poly"),
    clampNegToZero = true, mulImpl = polyMulImpl))

  polynomial.io.in.valid       := validTable
  polynomial.io.in.bits.coeffs := coefTable

  if(nOrder != 0) {
    polynomial.io.in.bits.x.get := sdxTable.get
  }

  // the sign bypasses the polynomial, delayed by its latency
  private val cutPoly = new SharedPipe("poly", validTable, stages.sub("poly").total)

  val sgnPoly = cutPoly(sgnTable)

  // -----------------------------------------------------------------
  // S3 : the magnitude, rounded to zSpec

  // the guard bits are rounded away (half up); a carry beyond the bits of
  // zSpec saturates to the largest magnitude
  val magW = zSpec.W - 1
  val rounded = if(guardBits == 0) {
    polynomial.io.out.bits
  } else {
    val res = polynomial.io.out.bits
    (res >> guardBits) +& res(guardBits-1)
  }
  val mag = if(rounded.getWidth > magW) {
    Mux(rounded(magW) === 1.U, Fill(magW, 1.U(1.W)), rounded(magW-1, 0))
  } else {
    rounded.pad(magW)
  }

  // --------------------------- stage "round" ---------------------------
  private val cutRound = new SharedPipe("round", polynomial.io.out.valid, stages.at("round"))

  val magRound   = cutRound(mag)
  val sgnRound   = cutRound(sgnPoly)
  val validRound = cutRound.valid

  // -----------------------------------------------------------------
  // S4 : the sign; the magnitude stays below 2^(intW-1), so its negation fits

  val magZ = Cat(0.U(1.W), magRound)
  val z    = Mux(sgnRound === 1.U, 0.U(zSpec.W.W) - magZ, magZ)
  require(z.getWidth == zSpec.W)

  // ---------------------------- stage "out" ----------------------------
  io.z := Pipe(validRound, z, stages.at("out"))
}

/** Pipeline stage points and the tail of [[rial.rng.GenRandomGaussianFixed]].
 *
 * @see [[rial.util.StagePoint]]
 * @see [[rial.util.StageGroup]]
 */
object GenRandomGaussianFixed {

  /** The octaves of the tail the last row spans. One more would reach
   *  further at the cost of the coefficient widths of the whole table. */
  private[rial] val lumpOctaves = 1

  /** `2^(intW-1) - 2^-fracW`, the largest magnitude `zSpec` holds. */
  private[rial] def zCap(zSpec: FixedSpec): Double = {
    pow(2.0, zSpec.intW - 1) - pow(2.0, -zSpec.fracW)
  }

  /** The width of the position, `zSpec.fracW + 4` unless `posW` gives one. */
  def positionWidth(zSpec: FixedSpec, posW: Option[Int] = None): Int = {
    posW.getOrElse(zSpec.fracW + 4)
  }

  /** The number of octaves of the tail the module reads out of a `rndW`-bit
   * word: one per bit of the scan field, up to the octave at which the result
   * reaches the largest magnitude of `zSpec`.
   *
   * @param rndW  input random number bit width.
   * @param zSpec output fixed point spec.
   * @param posW  width of the position, as given to the module.
   */
  def octaves(rndW: Int, zSpec: FixedSpec, posW: Option[Int] = None): Int = {
    val scanW = rndW - 2 - positionWidth(zSpec, posW)
    val zcap  = zCap(zSpec)
    // every x below 2^-(n+1) gives |z| >= zCap once this does
    (0 until scanW).find(n => {
      ProbitTable.probitHalf(pow(2.0, -(n+1))) >= zcap
    }).getOrElse(scanW)
  }

  /** The largest |z| the module of these parameters produces. No result
   * beyond it occurs.
   *
   * @param rndW      input random number bit width.
   * @param zSpec     output fixed point spec.
   * @param nOrder    order of the polynomial of a row.
   * @param splitW    sub-intervals of a row, as given to the module.
   * @param guardBits guard bits of the polynomial, as given to the module.
   * @param posW      width of the position, as given to the module.
   */
  def reach(rndW: Int, zSpec: FixedSpec, nOrder: Int = 3, splitW: Int = 3,
    guardBits: Int = 3, posW: Option[Int] = None): Double = {
    val t = table(zSpec, nOrder, octaves(rndW, zSpec, posW),
      positionWidth(zSpec, posW), splitW, guardBits)
    largestMagnitude(t).toDouble * pow(2.0, -zSpec.fracW)
  }

  // building a table solves erfc over every row, so it is made once per
  // configuration
  private val tables =
    scala.collection.mutable.Map.empty[(FixedSpec, Int, Int, Int, Int, Int), ProbitFixedTable]

  private[rial] def table(zSpec: FixedSpec, nOrder: Int, nOctaves: Int,
    posW: Int, splitW: Int, guardBits: Int): ProbitFixedTable = {
    tables.synchronized {
      tables.getOrElseUpdate((zSpec, nOrder, nOctaves, posW, splitW, guardBits), {
        new ProbitFixedTable(zSpec, nOrder, nOctaves, posW, splitW, guardBits)
      })
    }
  }

  /** The magnitude, in units of `2^-fracW`, of a polynomial output: the guard
   *  bits rounded away half up, the carry saturated to the largest magnitude.
   */
  private[rial] def roundGuard(polyOut: BigInt, zSpec: FixedSpec, guardBits: Int): BigInt = {
    val magMax = (BigInt(1) << (zSpec.W - 1)) - 1
    if (guardBits == 0) { polyOut } else {
      val rounded = (polyOut >> guardBits) + ((polyOut >> (guardBits - 1)) & 1)
      if (rounded > magMax) { magMax } else { rounded }
    }
  }

  /** The largest magnitude, in units of `2^-fracW`, that `t` gives at any
   *  sub-interval and position.
   *
   *  A Horner step multiplies by a position in [-1, 1) and truncates, so no
   *  sub-interval exceeds its constant coefficient plus the magnitudes of the
   *  others, plus one LSB per step. The sub-intervals are read at every
   *  position in the order of that bound, until the bound of the next one is
   *  no more than the largest result found.
   */
  private[rial] def largestMagnitude(t: ProbitFixedTable): BigInt = {
    val bitW   = t.coeffs.bitW
    val nOrder = bitW.nOrder
    val nPos   = if (nOrder == 0) { 1 } else { 1 << bitW.inBitW }
    val magMax = (BigInt(1) << (t.zSpec.W - 1)) - 1

    def bound(cs: Seq[Long]): BigInt = {
      val b = BigInt(cs.head) + cs.tail.map(c => BigInt(c).abs).sum + nOrder
      roundGuard(b.max(BigInt(0)).min((BigInt(1) << bitW.outBitW) - 1),
        t.zSpec, t.guardBits)
    }

    val rows = t.coeffs.coeffs.indices.sortBy(i => -bound(t.coeffs.coeffs(i)))
    var best = BigInt(0)
    val it   = rows.iterator
    var done = false
    while (!done && it.hasNext) {
      val i = it.next()
      if (best >= magMax || bound(t.coeffs.coeffs(i)) <= best) {
        done = true
      } else {
        val res = (0 until nPos).iterator.map(dx => evalRow(i, dx, t.coeffs)).max
        best = best.max(roundGuard(
          PolynomialEval.clamp(res, bitW, clampNegToZero = true), t.zSpec, t.guardBits))
      }
    }
    best
  }

  /** [[rial.table.PolynomialEval.eval]] of one row, in Long where the widths
   *  allow it. */
  private def evalRow(row: Int, dx: Int, coeffs: PolynomialCoefficients): BigInt = {
    val bitW = coeffs.bitW
    if (bitW.nOrder == 0 || bitW.inBitW + bitW.calcBits.max > 62) {
      PolynomialEval.eval(row, BigInt(dx), coeffs)
    } else {
      val cs = coeffs.coeffs(row)
      val sx = dx.toLong - (1L << bitW.inPrec)
      var calcW = bitW.calcBits(bitW.nOrder)
      var sum   = cs(bitW.nOrder)
      for (k <- (bitW.nOrder - 1) to 0 by -1) {
        val xminW = min(calcW, bitW.inBitW - 1)
        sum   = (((sx >> (bitW.inBitW - 1 - xminW)) * sum) >> xminW) + cs(k)
        calcW = bitW.calcBits(k)
      }
      BigInt(sum)
    }
  }

  /** Boundaries at which pipeline registers can be placed; see
   * [[rial.rng.GenRandomGaussianFixed]]. The boundaries of the polynomial
   * appear as the group `"poly"`.
   *
   * @param nOrder      order of the polynomial of a row.
   * @param polyMulImpl structure of the products of the polynomial.
   */
  def stagePoints(nOrder: Int, polyMulImpl: MulImpl = MulImpl.plain): Seq[StageElem] = {
    Seq(
      StagePoint("decode", 2),
      StagePoint("table", 3),
      StageGroup("poly", PolynomialEval.stagePoints(nOrder, polyMulImpl)),
      StagePoint("round", 1),
      StagePoint("out", 1)
    )
  }
}
