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
 * The uniform input is split into a sign, a page bit and a uniform value `x`;
 * the octave of `x` selects an interval of the probit (inverse gaussian CDF),
 * whose polynomial is evaluated by [[rial.table.PolynomialEval]] and scaled by
 * the value of the probit at the shallow end of the interval.
 *
 * The latency of this module can be configured via
 * [[rial.util.PipelineStageConfig]]. The generator declares the boundaries of
 * its datapath, so the registers can be placed inside it - either
 * automatically or by name:
 *
 *  - `"in"`     the module input
 *  - `"table"`  after the interval decode, the coefficient table read and the
 *               position `dx` inside the interval, i.e. everything that feeds
 *               the polynomial
 *  - `"poly.*"` the boundaries of [[rial.table.PolynomialEval]], which
 *               evaluates the polynomial of the interval: `"poly.step1"` ...
 *               `"poly.stepN"` and `"poly.out"`; the group name `"poly"`
 *               stands for `"poly.out"`
 *  - `"mul"`    after the scale table read, the exponent / mantissa select and
 *               the mantissa product
 *  - `"out"`    the module output; the normalisation of the product and the
 *               assembly of the result belong to the logic just before it
 *
 * The signals that bypass the polynomial - the page, the table address and the
 * exponent and the mantissa of the uniform value - are delayed by the latency
 * of the `"poly"` group, so that they stay aligned with its result.
 *
 * {{{
 * // four registers, placed by the module
 * new GenRandomGaussian(rndW, spec, nOrder, 4)
 * // hand-picked positions, one of them inside the polynomial
 * new GenRandomGaussian(rndW, spec, nOrder,
 *   PipelineStageConfig.specify("in" -> 1, "poly.step1" -> 1, "mul" -> 1))
 * // combinational
 * new GenRandomGaussian(rndW, spec, nOrder, PipelineStageConfig.none)
 * }}}
 *
 * @constructor        create a chisel Module.
 * @param rndW         input random number bit width.
 * @param spec         output floating-point spec.
 * @param nOrder       order of the polynomial of an interval.
 * @param stage        pipeline configuration of this module.
 * @param lumpOctaves  deepest tail octaves merged into one interval.
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class GenRandomGaussian(
  rndW: Int,
  spec: RealSpec,
  nOrder:      Int,
  stage:       PipelineStageConfig,
  lumpOctaves: Int = 2, // deepest tail octaves merged into one interval
) extends Module {

  /** Boundaries of the datapath at which pipeline registers can be placed.
   * The weights are rough, target-independent estimates of the delay of the
   * logic that ends at each boundary. */
  val stagePoints = GenRandomGaussian.stagePoints(nOrder)

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  require(spec.W+1 < rndW)

  val io = IO(new Bundle {
    val rnd = Flipped(Valid(UInt(rndW.W)))
    val z   = Valid(UInt(spec.W.W))
  })

  /** The valid flag follows the data through the same cuts. */
  private def cutValid(v: Bool, n: Int): Bool = {
    ShiftRegister(v, n, false.B, true.B)
  }

  // ---------------------------------------------------------------
  // S0 : input registers ("in")

  val rndIn   = ShiftRegister(io.rnd.bits, stages.at("in"))
  val validIn = cutValid(io.rnd.valid, stages.at("in"))

  // ---------------------------------------------------------------
  // S1 : interval decode, coefficient table and position in the interval

  val sgn  = rndIn(rndW-1)
  val page = rndIn(rndW-2)
  val xex0 = rndIn(rndW-3, spec.manW)
  val xman = rndIn(spec.manW-1, 0)

  // Main range (xex0 != 0): the lowest set bit of xex0 selects the octave.
  // 0bxxxxxx1 -> 0          : 1/ 4~1/2 : exBias-2
  // 0bxxxxx10 -> 1          : 1/ 8~1/4 : exBias-3
  // ...
  val xex = Mux(xex0 === 0.U, xex0.getWidth.U, PriorityEncoder(xex0))
  // mantissa part for xex0===0.U
  val clz = Mux(xman === 0.U, spec.manW.U, PriorityEncoder(Reverse(xman)))

  // ---------------------------------------------------------------

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
  val lumpAddr    = (nIntervals - 1).U

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
  val addr = WireDefault(0.U(log2Ceil(nIntervals).W))
  when(page === 0.U) {
    addr := 0.U
  }.elsewhen(xex0 =/= 0.U) {
    addr := xex +& 1.U
  }.elsewhen(clz < deepIndiv.U) {
    addr := (mainOctaves + 1).U +& clz
  }.otherwise {
    addr := lumpAddr
  }
  val coef = packed(addr).asTypeOf( MixedVec(bitW.coefBits.map(w => SInt(w.W))) )

  // position inside the interval, as a signed value in [-1, 1)
  val sdx = if(nOrder != 0) {
    val dx = WireDefault(0.U(inBitW.W))
    when(page === 0.U) {
      dx := (Cat(1.U(1.W), xman(spec.manW-1, 1)) >> xex)
    }.elsewhen(xex0 =/= 0.U) {
      dx := xman
    }.elsewhen(clz < deepIndiv.U) {
      // deep octave: drop xman's leading zeros and leading 1, MSB-align the rest
      dx := (xman << (clz +& 1.U))(inBitW-1, 0)
    }.otherwise {
      // lumped interval: the lowest `lumpOctaves` bits of xman span it linearly
      dx := (xman << (spec.manW - lumpOctaves))(inBitW-1, 0)
    }
    Some(Cat(~dx(inBitW-1), dx(inBitW-2, 0)).asSInt)
  } else {
    None
  }

  // --------------------------- stage "table" ---------------------------
  private def cutTable[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("table"))
  }

  val coefTable  = cutTable(coef)
  val sdxTable   = sdx.map(x => cutTable(x))
  val pageTable  = cutTable(page)
  val addrTable  = cutTable(addr)
  val xexTable   = cutTable(xex)
  val xmanTable  = cutTable(xman)
  val sgnTable   = cutTable(sgn)
  val validTable = cutValid(validIn, stages.at("table"))

  // ---------------------------------------------------------------
  // S2 : polynomial of the interval ("poly.*")

  // we need to clamp negative value. since the last part of the probit is steep,
  // the 3rd order polynomial sometimes fail to fit around x~0.
  val polynomial = Module(new PolynomialEval(bitW, stages.sub("poly"), clampNegToZero = true))

  polynomial.io.in.valid       := validTable
  polynomial.io.in.bits.coeffs := coefTable

  if(nOrder != 0) {
    polynomial.io.in.bits.x.get := sdxTable.get
  }

  // latency of the polynomial; the signals that bypass it are delayed by the
  // same number of cycles, so they stay aligned with polynomial.io.out
  val nPoly = stages.sub("poly").total

  val postPage = ShiftRegister(pageTable, nPoly)
  val postAddr = ShiftRegister(addrTable, nPoly)
  val postEx   = ShiftRegister(xexTable,  nPoly)
  val postMan  = ShiftRegister(xmanTable, nPoly)
  val postSgn  = ShiftRegister(sgnTable,  nPoly)

  // -----------------------------------------------------------------
  // S3 : postproc (mul)

  // 1.0~2.0
  val yManW1 = Cat(1.U(1.W), polynomial.io.out.bits)

  // page=1
  val scaleTable = VecInit(table.scaleTable.map(_.U(spec.W.W)).toSeq)
  val scale      = scaleTable(postAddr)
  val scaleEx    = scale(spec.manW+spec.exW-1, spec.manW)
  val scaleMan   = scale(spec.manW-1, 0)

  // page=0
  val mulEx    = Mux(postPage === 0.U, (spec.exBias-2).U(spec.exW.W) - postEx, scaleEx)
  val mulMan   = Mux(postPage === 0.U, postMan, scaleMan)
  val mulManW1 = Cat(1.U(1.W), mulMan)

  val prodW = (1+spec.manW)*2
  val prod  = yManW1 * mulManW1

  // ---------------------------- stage "mul" ----------------------------
  private def cutMul[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("mul"))
  }

  val prodMul  = cutMul(prod)
  val mulExMul = cutMul(mulEx)
  val sgnMul   = cutMul(postSgn)
  val validMul = cutValid(polynomial.io.out.valid, stages.at("mul"))

  // -----------------------------------------------------------------
  // S4 : normalise the product and assemble the result

  val prodMoreThan2 = prodMul.head(1) === 1.U
  val prodManW1     = Mux(prodMoreThan2,
    prodMul(prodW-1, prodW-1-spec.manW),
    prodMul(prodW-2, prodW-2-spec.manW)
  )

  val zsgn = sgnMul
  val zex  = mulExMul + prodMoreThan2
  val zman = prodManW1(spec.manW-1, 0)

  val z = Cat(zsgn, zex, zman)
  require(z.getWidth == spec.W)

  // ---------------------------- stage "out" ----------------------------
  io.z.bits  := ShiftRegister(z, stages.at("out"))
  io.z.valid := cutValid(validMul, stages.at("out"))
}

/** Pipeline stage points of [[rial.rng.GenRandomGaussian]].
 *
 * @see [[rial.util.StagePoint]]
 * @see [[rial.util.StageGroup]]
 */
object GenRandomGaussian {

  /** Ordered boundaries of the datapath, with a rough, target-independent
   * estimate of the delay of the logic ending at each. The polynomial of the
   * interval is a submodule, so its own boundaries appear as the group
   * `"poly"`, i.e. as `"poly.step1"` ... `"poly.stepN"` and `"poly.out"`.
   *
   * @param nOrder order of the polynomial of an interval.
   */
  def stagePoints(nOrder: Int): Seq[StageElem] = Seq(
    StagePoint("table", 3),
    StageGroup("poly", PolynomialEval.stagePoints(nOrder)),
    StagePoint("mul", 4),
    StagePoint("out", 1))
}
