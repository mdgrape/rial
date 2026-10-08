//% @file GenRandomGaussian.scala
//
// Copyright (C) Toru Niina RIKEN BDR 2026
//
package rial.rng

import scala.math._
import scala.collection.mutable.ArrayBuffer

import chisel3._
import chisel3.util._

import rial.arith._
import rial.table._
import rial.util._

import org.apache.commons.math3.special.Erf

import org.hipparchus.analysis.UnivariateFunction
import org.hipparchus.analysis.solvers.BrentSolver

//
// inverse CDF of gaussian
//
private[rial] class ProbitTable(
  val spec:        RealSpec,
  val nOrder:      Int,
  val nIntervals:  Int,
  val inBitW:      Int,
  val outBitW:     Int,
  val lumpOctaves: Int = 2, // deepest octaves merged into one (capped) interval
  val splitW:      Int = 3, // every interval split into 2^splitW sub-intervals
) {
  private val sqrt2 = scala.math.sqrt(2)

  val scaleTable = ArrayBuffer.empty[BigInt]

  val ft = new FuncTable(nOrder)

  val nSub = 1 << splitW

  // Every interval is divided into `nSub` equal sub-intervals of x, each with
  // its own polynomial, in the order of the position inside the interval: the
  // sub-interval 0 holds the xMin end, the one the position 0 reads.
  private def addSplit(f: Double => Double, xMin: Double, xMax: Double): Unit = {
    for(s <- 0 until nSub) {
      ft.addInterval(f, xMin + (xMax - xMin) * s / nSub,
                        xMin + (xMax - xMin) * (s+1) / nSub)
    }
  }

  addSplit((x: Double) => { sqrt2 * Erf.erfInv(x) / x - 1.0 }, 0.0, 0.5)
  scaleTable += BigInt(0)

  require(nIntervals > 2)
  // interval 0 is the central page, 1..nOctaves the octaves, the last one the
  // lumped deep tail
  val nOctaves = nIntervals - 2

  // erfcInv(x) is calculated as erfInv(1-x), and it becomes Inf if x is small
  // enough. instead, we will explicitly solve erfc - x.
  // since erfc(x) -> 0 if x -> 0, solver mistakenly judges its already solved
  // when x -> 0. we avoid this misjudge by applying log to both side.
  private val erfcSolver = new BrentSolver(1e-14, 1e-13)
  val probitHalf = (x: Double) => {
    val logx = log(x)
    val f: UnivariateFunction = (t: Double) => log(Erf.erfc(t)) - logx
    sqrt2 * erfcSolver.solve(200, f, 0.0, 25.0)
  }
  for(i <- 1 to nOctaves) {
    val xMin  = pow(2.0, -i)
    val xMax  = pow(2.0, -(i+1))
    val fMin  = probitHalf(xMin)
    val fMax  = probitHalf(xMax)
    assert(fMax / fMin < 2.0)

    val f = (x: Double) => { probitHalf(x) / fMin - 1.0 }
    addSplit(f, xMin, xMax)
    scaleTable += (new RealGeneric(spec, fMin)).value.toBigInt
  }

  // The deepest `lumpOctaves` octaves (below the individual ones) are merged
  // into a single capped interval, x in [2^-(nOctaves+lumpOctaves+1),
  // 2^-(nOctaves+1)]. A single cubic over a few octaves is a slight
  // approximation (see GenGaussianTailProbe); lumpOctaves trades that error for
  // table size. The cap (min with fMax) bounds |z| at the deepest x.
  val xMin  = pow(2.0, -(nOctaves+1))
  val xMax  = pow(2.0, -(nOctaves+lumpOctaves+1))
  val fMin  = probitHalf(xMin)
  val fMax  = probitHalf(xMax)
  assert(fMax / fMin < 2.0, f"lump xMin=${xMin}, xMax=${xMax}, fMin=${fMin}, fMax=${fMax}, fMax/fMin=${fMax/fMin}")

  val f = (x: Double) => { min(probitHalf(x), fMax) / fMin - 1.0 }
  addSplit(f, xMin, xMax)
  scaleTable += (new RealGeneric(spec, fMin)).value.toBigInt

  assert(ft.intervals.length == nIntervals * nSub)

  // ------------------------------------------------------------------------

  // the probit is shallow at the top end, where a row of its table reaches
  // slightly below zero; GenRandomGaussian reads it with an evaluator that
  // clamps a negative result to zero
  private val gen = new PolynomialGenerator(ft, inBitW, outBitW,
    negativeIsClamped = true)

  def getBitWidth(): PolynomialBitWidth = {
    gen.getBitWidth()
  }
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients = {
    gen.generateTable(bitW)
  }
}

// recommended rndW = 2*manW (64 for FP32)
/** Generates a gaussian random number out of `rndW` uniform random bits.
 *
 * The generator is based on a piecewise polynomial approximation of the
 * inverse CDF of the gaussian, evaluated by [[rial.table.PolynomialEval]].
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"decode"` after the interval decode
 *  - `"table"`  after the coefficient table lookup
 *  - `"poly.*"` the boundaries of the [[rial.table.PolynomialEval]] submodule
 *               that evaluates the polynomial of the interval, see its doc for
 *               the names; the group name `"poly"` stands for `"poly.out"`
 *  - `"scale"`  after the scale lookup
 *  - `"mul.*"`  the boundaries of the mantissa product: with `MulImpl.split`
 *               those of [[rial.arith.SplitMultiplier]], and `"mul.out"`
 *               alone for a plain product; the group name `"mul"` stands for
 *               `"mul.out"`
 *  - `"out"`    the module output
 *
 * [[getStage]] gives the resulting latency, and `z.valid` is `rnd.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // four registers, placed by the module
 * new GenRandomGaussian(rndW, spec, nOrder, 4)
 * // hand-picked positions, one of them inside the polynomial
 * new GenRandomGaussian(rndW, spec, nOrder,
 *   PipelineStageConfig.specify("in" -> 1, "poly.prod1" -> 1, "mul" -> 1))
 * // combinational
 * new GenRandomGaussian(rndW, spec, nOrder, PipelineStageConfig.none)
 * }}}
 *
 * @param rndW         input random number bit width.
 * @param spec         output floating-point spec.
 * @param nOrder       order of the polynomial of an interval.
 * @param stage        pipeline registers of this module.
 * @param lumpOctaves  deepest tail octaves merged into one interval.
 * @param splitW       every interval of the table is divided into `2^splitW`
 *                     sub-intervals, each with its own polynomial.
 * @param guardBits    bits below the mantissa that the polynomial is computed
 *                     with, rounded away before the scale multiply.
 * @param polyMulImpl  structure of the products of the polynomial. Only this
 *                     one has an addend to give a multiply-add, the
 *                     coefficient of a Horner step.
 * @param scaleMulImpl structure of the product of the mantissas.
 *
 * @see [[rial.arith.MulImpl]]
 * @see [[rial.util.PipelineStageConfig]]
 */
class GenRandomGaussian(
  rndW: Int,
  spec: RealSpec,
  nOrder:      Int = 3,
  stage:       PipelineStageConfig,
  lumpOctaves: Int = 2,
  splitW:      Int = 3,
  guardBits:   Int = 3,
  polyMulImpl:  MulImpl = MulImpl.plain,
  scaleMulImpl: MulImpl = MulImpl.plain,
) extends Module {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints =
    GenRandomGaussian.stagePoints(nOrder, polyMulImpl, scaleMulImpl)

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  require(spec.W+1 < rndW)

  val io = IO(new Bundle {
    val rnd = Flipped(Valid(UInt(rndW.W)))
    val z   = Valid(UInt(spec.W.W))
  })

  // ---------------------------------------------------------------
  // S0 : input registers ("in")

  private val cutIn = new SharedPipe("in", io.rnd.valid, stages.at("in"))

  val rndIn   = cutIn(io.rnd.bits)
  val validIn = cutIn.valid

  // ---------------------------------------------------------------
  // S1 : interval decode, coefficient table and position in the interval

  val manW = spec.manW

  val sgn  = rndIn(rndW-1)
  val page = rndIn(rndW-2)
  val xex0 = rndIn(rndW-3, manW)
  val xman = rndIn(manW-1, 0)

  // the normalised position inside an interval is manW bits wide: its top
  // splitW bits choose the sub-interval, the rest is the polynomial input
  val inBitW  = manW - splitW
  val outBitW = manW + guardBits
  require(splitW >= 0 && inBitW >= 2, "the position inside a sub-interval needs at least two bits")
  require(guardBits >= 0)

  // Intervals: page-0 (idx 0) + nOctaves individual octaves (idx 1..nOctaves) +
  // 1 lumped deep interval (idx nOctaves+1). The main range provides
  // xex0.getWidth octaves; the deep tail provides manW octaves, of which
  // the deepest `lumpOctaves` are merged, leaving (manW - lumpOctaves)
  // individual deep octaves.
  val mainOctaves = xex0.getWidth                          // e.g. 39
  val deepIndiv   = manW - lumpOctaves                     // e.g. 21
  val nOctaves    = mainOctaves + deepIndiv                // e.g. 60
  val nIntervals  = 1 + nOctaves + 1                       // e.g. 62
  val addrW       = log2Ceil(nIntervals)                   // e.g. 6
  val lumpAddr    = (nIntervals - 1).U(addrW.W)

  val table  = new ProbitTable(spec, nOrder, nIntervals, inBitW, outBitW,
    lumpOctaves, splitW)
  val bitW   = table.getBitWidth()
  val coeffs = table.generateTable(bitW)

  val packed = VecInit(PolynomialGenerator.packTable(coeffs).map(c => {
    c.U(coeffs.bitW.coefBits.reduce(_+_).W)
  }))

  // octave i (x in [2^-(i+1), 2^-i]) lives at table index i:
  //   page 0                  -> idx 0
  //   page 1, xex0 != 0       -> idx xex+1            (main octaves 1..mainOctaves)
  //   page 1, xex0 == 0, deep -> idx mainOctaves+1+clz (deep octaves), clamped to
  //                              the lumped interval once clz reaches deepIndiv.
  //
  // The main range takes the octave from the lowest set bit of xex0
  //   0bxxxxxx1 -> 0          : 1/ 4~1/2 : exBias-2
  //   0bxxxxx10 -> 1          : 1/ 8~1/4 : exBias-3
  //   ...
  // and the deep tail from the leading zeros of xman, so one scan for the lowest
  // set bit of the two of them, laid end to end, tells the octave of either.
  // The scans are grouped, and they select the table index and the shift of the
  // position rather than the octave number, so that the deepest paths of this
  // stage are two selections rather than a priority chain over the whole word
  // followed by an address decode and a shifter control.
  val scan = Cat(Reverse(xman), xex0)

  // the deep tail, saturated at the lumped interval
  val deepAddr = GenRandomGaussian.firstSetSelect(Reverse(xman), lumpAddr)(i => {
    min(mainOctaves + 1 + i, nIntervals - 1).U(addrW.W)
  })
  // the main range, with the central page as one more, lowest, position and the
  // deep tail as the result when no octave bit of it is set, so that the
  // address needs no selection of its own
  val addr = GenRandomGaussian.firstSetSelect(Cat(xex0, ~page), deepAddr)(i => {
    i.U(addrW.W)
  })

  // Shift of the position inside the interval. `dxBits` is 1.xman followed by
  // manW zeros, and every interval takes its position as the manW-bit window
  // of it that starts at `dxShift`:
  //   page 0            : 1.xman >> 1 >> xex, i.e. the window at manW+1+xex
  //   page 1, main      : xman itself, the window at manW
  //   page 1, deep      : xman without its leading zeros and leading one, the
  //                       window at manW-1-clz
  //   page 1, lumped    : the lowest lumpOctaves bits of xman, spanning the
  //                       interval linearly, the window at manW-deepIndiv
  val dxShiftW  = log2Ceil(manW + 2 + mainOctaves)
  val centralShift = (i: Int) => manW + 1 + min(i, mainOctaves)
  val octaveShift  = (i: Int) => manW - min(deepIndiv, max(0, i - mainOctaves + 1))
  val dxShiftOf = (i: Int) => {
    Mux(page === 0.U, centralShift(i).U(dxShiftW.W), octaveShift(i).U(dxShiftW.W))
  }
  val dxShift = GenRandomGaussian.firstSetSelect(scan, dxShiftOf(scan.getWidth))(dxShiftOf)

  // --------------------------- stage "decode" --------------------------
  private val cutDecode = new SharedPipe("decode", validIn, stages.at("decode"))

  val addrDecode  = cutDecode(addr)
  val shiftDecode = cutDecode(dxShift)
  val pageDecode  = cutDecode(page)
  val xmanDecode  = cutDecode(xman)
  val sgnDecode   = cutDecode(sgn)
  val validDecode = cutDecode.valid

  // the normalised position inside the interval
  val dxBits = Cat(1.U(1.W), xmanDecode, 0.U(manW.W))
  val dxFull = (dxBits >> shiftDecode)(manW-1, 0)

  // the sub-interval the position falls in sits below the interval index
  val romAddr = if(splitW != 0) {
    Cat(addrDecode, dxFull(manW-1, manW-splitW))
  } else {
    addrDecode
  }

  // the row is cut as one word, the shape a block RAM can register
  val coefWord = packed(romAddr)

  // position inside the sub-interval, as a signed value in [-1, 1)
  val sdx = if(nOrder != 0) {
    val dx = dxFull(inBitW-1, 0)
    Some(Cat(~dx(inBitW-1), dx(inBitW-2, 0)).asSInt)
  } else {
    None
  }

  // --------------------------- stage "table" ---------------------------
  private val cutTable = new SharedPipe("table", validDecode, stages.at("table"))

  val coefTable  = cutTable(coefWord).asTypeOf(MixedVec(bitW.coefBits.map(w => SInt(w.W))))
  val sdxTable   = sdx.map(x => cutTable(x))
  val pageTable  = cutTable(pageDecode)
  val addrTable  = cutTable(addrDecode)
  val shiftTable = cutTable(shiftDecode)
  val xmanTable  = cutTable(xmanDecode)
  val sgnTable   = cutTable(sgnDecode)
  val validTable = cutTable.valid

  // ---------------------------------------------------------------
  // S2 : polynomial of the interval ("poly.*")

  // we need to clamp negative value. since the last part of the probit is steep,
  // the 3rd order polynomial sometimes fail to fit around x~0.
  val polynomial = Module(new PolynomialEval(bitW, stages.sub("poly"),
    clampNegToZero = true, mulImpl = polyMulImpl))

  polynomial.io.in.valid       := validTable
  polynomial.io.in.bits.coeffs := coefTable

  if(nOrder != 0) {
    polynomial.io.in.bits.x.get := sdxTable.get
  }

  // latency of the polynomial; the signals that bypass it are delayed by the
  // same number of cycles, so they stay aligned with polynomial.io.out
  val nPoly = stages.sub("poly").total

  private val cutPoly = new SharedPipe("poly", validTable, nPoly)

  val postPage  = cutPoly(pageTable)
  val postAddr  = cutPoly(addrTable)
  val postShift = cutPoly(shiftTable)
  val postMan   = cutPoly(xmanTable)
  val postSgn   = cutPoly(sgnTable)

  // -----------------------------------------------------------------
  // S3 : scale table and the operands of the mantissa product

  // the guard bits of the polynomial are rounded away (half up), with the
  // carry into 1.0 clamped, as the math functions round theirs
  val polyMan = if(guardBits == 0) {
    polynomial.io.out.bits
  } else {
    val res   = polynomial.io.out.bits
    val zman0 = (res >> guardBits) +& res(guardBits-1)
    Mux(zman0(manW) === 1.U, Fill(manW, 1.U(1.W)), zman0(manW-1, 0))
  }

  // 1.0~2.0
  val yManW1 = Cat(1.U(1.W), polyMan)

  // page=1
  val scaleTable = VecInit(table.scaleTable.map(_.U(spec.W.W)).toSeq)
  val scale      = scaleTable(postAddr)
  val scaleEx    = scale(manW+spec.exW-1, manW)
  val scaleMan   = scale(manW-1, 0)

  // page=0; the shift of the position carries the octave, so the exponent
  // (exBias-2) - xex is the one of the window at manW+1+xex
  val mulEx    = Mux(postPage === 0.U,
    (spec.exBias - 1 + manW).U(spec.exW.W) - postShift, scaleEx)
  val mulMan   = Mux(postPage === 0.U, postMan, scaleMan)
  val mulManW1 = Cat(1.U(1.W), mulMan)

  // --------------------------- stage "scale" ---------------------------
  private val cutScale = new SharedPipe("scale", polynomial.io.out.valid, stages.at("scale"))

  val yManW1Scale   = cutScale(yManW1)
  val mulExScale    = cutScale(mulEx)
  // the exponent of a product that came out >= 2.0, ready so that the
  // normalisation of the next stage is a selection and not an increment
  val mulExUpScale  = cutScale(mulEx + 1.U)
  val mulManW1Scale = cutScale(mulManW1)
  val sgnScale      = cutScale(postSgn)
  val validScale    = cutScale.valid

  // -----------------------------------------------------------------
  // S4 : mantissa product ("mul.*")

  val prodW = (1+manW)*2

  require(!scaleMulImpl.isMacc,
    "the product of the mantissas has no addend to fold into a multiply-add; " +
    "ask for MulImpl.split to slice it")

  val prodMul = if (scaleMulImpl.isSplit) {
    val mult = Module(new SplitMultiplier(
      yManW1Scale.getWidth, mulManW1Scale.getWidth, scaleMulImpl.sliceW.get,
      stages.sub("mul")))

    mult.io.in.valid  := validScale
    mult.io.in.bits.x := yManW1Scale
    mult.io.in.bits.y := mulManW1Scale
    mult.io.out.bits
  } else {
    // the product as it is written by hand: `MulImpl.plain` asks nothing of
    // the operands, so there is nothing to shape
    val cutMulIn  = new SharedPipe("mulIn",  validScale,     stages.at("mul.in"))
    val cutMulOut = new SharedPipe("mulOut", cutMulIn.valid, stages.at("mul.out"))
    cutMulOut(cutMulIn(yManW1Scale) * cutMulIn(mulManW1Scale))
  }

  require(prodMul.getWidth == prodW,
    s"the product of the mantissas is ${prodMul.getWidth} bits, expected $prodW")

  // ---------------------------- stage "mul" ----------------------------
  // latency of the product; the signals that bypass it are delayed by the
  // same number of cycles, so they stay aligned with it
  private val cutMul = new SharedPipe("mul", validScale, stages.latency("mul"))

  val mulExMul   = cutMul(mulExScale)
  val mulExUpMul = cutMul(mulExUpScale)
  val sgnMul     = cutMul(sgnScale)
  val validMul   = cutMul.valid

  // -----------------------------------------------------------------
  // S5 : normalise the product and assemble the result

  val prodMoreThan2 = prodMul.head(1) === 1.U
  val prodManW1     = Mux(prodMoreThan2,
    prodMul(prodW-1, prodW-1-manW),
    prodMul(prodW-2, prodW-2-manW)
  )

  val zsgn = sgnMul
  val zex  = Mux(prodMoreThan2, mulExUpMul, mulExMul)
  val zman = prodManW1(manW-1, 0)

  val z = Cat(zsgn, zex, zman)
  require(z.getWidth == spec.W)

  // ---------------------------- stage "out" ----------------------------
  io.z := Pipe(validMul, z, stages.at("out"))
}

/** Pipeline stage points of [[rial.rng.GenRandomGaussian]].
 *
 * @see [[rial.util.StagePoint]]
 * @see [[rial.util.StageGroup]]
 */
object GenRandomGaussian {

  /** Bits of a lowest-set-bit scan that share one candidate value. */
  private val scanGroupW = 6

  /** Selects the value that the lowest set bit of `bits` stands for: bit `i`
   * gives `value(i)`, and `none` is the result when no bit is set.
   *
   * The bits are scanned in groups: every group offers the value of its own
   * lowest set bit, and the first non-empty group picks one of them, so that
   * the result is two selections deep rather than one priority chain over the
   * whole word.
   *
   * @param bits  the word to scan.
   * @param none  the result when `bits` is zero.
   * @param value the value of each bit position.
   */
  private[rial] def firstSetSelect(bits: UInt, none: UInt)(value: Int => UInt): UInt = {
    val nBits = bits.getWidth
    if (nBits <= scanGroupW) {
      PriorityMux(bits.asBools :+ true.B, (0 until nBits).map(value) :+ none)
    } else {
      val nGroups = (nBits + scanGroupW - 1) / scanGroupW
      val groups  = (0 until nGroups).map(g => {
        bits(min(nBits, (g + 1) * scanGroupW) - 1, g * scanGroupW)
      })
      val chosen = groups.zipWithIndex.map({ case (grp, g) => {
        PriorityMux(grp.asBools, (0 until grp.getWidth).map(j => {
          value(g * scanGroupW + j)
        }))
      }})
      firstSetSelect(VecInit(groups.map(_.orR)).asUInt, none)(g => chosen(g))
    }
  }

  /** Boundaries at which pipeline registers can be placed; see
   * [[rial.rng.GenRandomGaussian]], for a plain product.
   *
   * @param nOrder order of the polynomial of an interval.
   */
  def stagePoints(nOrder: Int): Seq[StageElem] = {
    stagePoints(nOrder, MulImpl.plain, MulImpl.plain)
  }

  /** Boundaries at which pipeline registers can be placed; see
   * [[rial.rng.GenRandomGaussian]]. The boundaries of the polynomial and of
   * the mantissa product appear as the groups `"poly"` and `"mul"`.
   *
   * @param nOrder       order of the polynomial of an interval.
   * @param polyMulImpl  structure of the products of the polynomial.
   * @param scaleMulImpl structure of the product of the mantissas.
   */
  def stagePoints(
    nOrder: Int, polyMulImpl: MulImpl, scaleMulImpl: MulImpl
  ): Seq[StageElem] = {
    Seq(
      StagePoint("decode", 3),
      StagePoint("table", 3),
      StageGroup("poly", PolynomialEval.stagePoints(nOrder, polyMulImpl)),
      StagePoint("scale", 2),
      StageGroup("mul", if (scaleMulImpl.isSplit) {
        SplitMultiplier.stagePoints()
      } else {
        Seq(StagePoint("out", 4.0))
      }),
      StagePoint("out", 1)
    )
  }
}
