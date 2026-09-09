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

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

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
) {
  private val sqrt2 = scala.math.sqrt(2)
  private val log2 = (a: Double) => log(a) / log(2.0)

  val scaleTable = ArrayBuffer.empty[BigInt]

  val ft = new FuncTable(nOrder)

  ft.addInterval((x: Double) => { sqrt2 * Erf.erfInv(x) / x - 1.0 }, 0.0, 0.5)
  scaleTable += 0.toBigInt

  require(nIntervals > 2)
  // 0: center part
  // 1~nOctaves-1: normal polynomial
  // nOctaves: the last subnormal

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
    ft.addInterval(f, xMin, xMax)
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
  ft.addInterval(f, xMin, xMax)
  scaleTable += (new RealGeneric(spec, fMin)).value.toBigInt

  val maxValue = new RealGeneric(spec, fMax)

  assert(ft.intervals.length == nIntervals)

  // ------------------------------------------------------------------------

  private val gen = new PolynomialGenerator(ft, inBitW, outBitW)

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
 *  - `"mul.*"`  the boundaries of the [[rial.arith.Multiplier]] that forms the
 *               mantissa product, see its doc for the names; the group name
 *               `"mul"` stands for `"mul.out"`
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
 * @constructor        create a chisel Module.
 * @param rndW         input random number bit width.
 * @param spec         output floating-point spec.
 * @param nOrder       order of the polynomial of an interval.
 * @param stage        pipeline registers of this module.
 * @param lumpOctaves  deepest tail octaves merged into one interval.
 * @param polyMulImpl  structure of the products of the polynomial. Only this
 *                     one has an addend to give a multiply-add, the
 *                     coefficient of a Horner step.
 * @param scaleMulImpl structure of the product of the mantissas.
 *
 * @see [[rial.arith.MulImpl]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class GenRandomGaussian(
  rndW: Int,
  spec: RealSpec,
  nOrder:      Int,
  stage:       PipelineStageConfig,
  lumpOctaves: Int = 2, // deepest tail octaves merged into one interval
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

  val sgn  = rndIn(rndW-1)
  val page = rndIn(rndW-2)
  val xex0 = rndIn(rndW-3, spec.manW)
  val xman = rndIn(spec.manW-1, 0)

  val inBitW  = spec.manW
  val outBitW = spec.manW

  // Intervals: page-0 (idx 0) + nOctaves individual octaves (idx 1..nOctaves) +
  // 1 lumped deep interval (idx nOctaves+1). The main range provides
  // xex0.getWidth octaves; the deep tail provides spec.manW octaves, of which
  // the deepest `lumpOctaves` are merged, leaving (manW - lumpOctaves)
  // individual deep octaves.
  val mainOctaves = xex0.getWidth                          // e.g. 39
  val deepIndiv   = spec.manW - lumpOctaves                // e.g. 21
  val nOctaves    = mainOctaves + deepIndiv                // e.g. 60
  val nIntervals  = 1 + nOctaves + 1                       // e.g. 62
  val addrW       = log2Ceil(nIntervals)                   // e.g. 6
  val lumpAddr    = (nIntervals - 1).U(addrW.W)

  val table  = new ProbitTable(spec, nOrder, nIntervals, inBitW, outBitW, lumpOctaves)
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
  // inBitW zeros, and every interval takes its position as the inBitW-bit
  // window of it that starts at `dxShift`:
  //   page 0            : 1.xman >> 1 >> xex, i.e. the window at inBitW+1+xex
  //   page 1, main      : xman itself, the window at inBitW
  //   page 1, deep      : xman without its leading zeros and leading one, the
  //                       window at inBitW-1-clz
  //   page 1, lumped    : the lowest lumpOctaves bits of xman, spanning the
  //                       interval linearly, the window at inBitW-deepIndiv
  val dxShiftW  = log2Ceil(inBitW + 2 + mainOctaves)
  val centralShift = (i: Int) => inBitW + 1 + min(i, mainOctaves)
  val octaveShift  = (i: Int) => inBitW - min(deepIndiv, max(0, i - mainOctaves + 1))
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

  val coef = packed(addrDecode).asTypeOf( MixedVec(bitW.coefBits.map(w => SInt(w.W))) )

  // position inside the interval, as a signed value in [-1, 1)
  val sdx = if(nOrder != 0) {
    val dxBits = Cat(1.U(1.W), xmanDecode, 0.U(inBitW.W))
    val dx = (dxBits >> shiftDecode)(inBitW-1, 0)
    Some(Cat(~dx(inBitW-1), dx(inBitW-2, 0)).asSInt)
  } else {
    None
  }

  // --------------------------- stage "table" ---------------------------
  private val cutTable = new SharedPipe("table", validDecode, stages.at("table"))

  val coefTable  = cutTable(coef)
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

  // 1.0~2.0
  val yManW1 = Cat(1.U(1.W), polynomial.io.out.bits)

  // page=1
  val scaleTable = VecInit(table.scaleTable.map(_.U(spec.W.W)).toSeq)
  val scale      = scaleTable(postAddr)
  val scaleEx    = scale(spec.manW+spec.exW-1, spec.manW)
  val scaleMan   = scale(spec.manW-1, 0)

  // page=0; the shift of the position carries the octave, so the exponent
  // (exBias-2) - xex is the one of the window at inBitW+1+xex
  val mulEx    = Mux(postPage === 0.U,
    (spec.exBias - 1 + inBitW).U(spec.exW.W) - postShift, scaleEx)
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

  val prodW = (1+spec.manW)*2

  // the product of the mantissas has no addend of its own, so a multiply-add
  // gets zero for one
  val mult = Module(new Multiplier(
    yManW1Scale.getWidth, mulManW1Scale.getWidth, scaleMulImpl,
    stages.sub("mul"), cW = Multiplier.zeroAddendW(scaleMulImpl)))

  mult.io.in.valid  := validScale
  mult.io.in.bits.x := yManW1Scale
  mult.io.in.bits.y := mulManW1Scale
  mult.io.in.bits.c.foreach(c => c := 0.U)

  // ---------------------------- stage "mul" ----------------------------
  // latency of the product; the signals that bypass it are delayed by the
  // same number of cycles, so they stay aligned with mult.io.out
  private val cutMul = new SharedPipe("mul", validScale, stages.latency("mul"))

  val prodMul    = mult.io.out.bits(prodW-1, 0)
  val mulExMul   = cutMul(mulExScale)
  val mulExUpMul = cutMul(mulExUpScale)
  val sgnMul     = cutMul(sgnScale)
  val validMul   = cutMul.valid

  // -----------------------------------------------------------------
  // S5 : normalise the product and assemble the result

  val prodMoreThan2 = prodMul.head(1) === 1.U
  val prodManW1     = Mux(prodMoreThan2,
    prodMul(prodW-1, prodW-1-spec.manW),
    prodMul(prodW-2, prodW-2-spec.manW)
  )

  val zsgn = sgnMul
  val zex  = Mux(prodMoreThan2, mulExUpMul, mulExMul)
  val zman = prodManW1(spec.manW-1, 0)

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
  def stagePoints(nOrder: Int, polyMulImpl: MulImpl,
    scaleMulImpl: MulImpl): Seq[StageElem] = Seq(
    StagePoint("decode", 3),
    StagePoint("table", 3),
    StageGroup("poly", PolynomialEval.stagePoints(nOrder, polyMulImpl)),
    StagePoint("scale", 2),
    StageGroup("mul", Multiplier.stagePoints(scaleMulImpl)),
    StagePoint("out", 1))
}
