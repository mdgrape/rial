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
class GenRandomGaussian(
  rndW: Int,
  spec: RealSpec,
  nOrder:          Int,
  polynomialStage: Int,
  mulStage:        Int,
  lumpOctaves:     Int = 2, // deepest tail octaves merged into one interval
) extends Module {

  val nStage = 1 + polynomialStage + mulStage

  require(spec.W+1 < rndW)

  val io = IO(new Bundle {
    val rnd = Flipped(Valid(UInt(rndW.W)))
    val z   = Valid(UInt(spec.W.W))
  })

  val rnd = Pipe(io.rnd, 1)

  val page = rnd.bits(rndW-2)
  val xex0 = rnd.bits(rndW-3, spec.manW)
  val xman = rnd.bits(spec.manW-1, 0)

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

  // we need to clamp negative value. since the last part of the probit is steep,
  // the 3rd order polynomial sometimes fail to fit around x~0.
  val polynomial = Module(new PolynomialEval(bitW, polynomialStage, clampNegToZero = true))

  polynomial.io.in.valid       := rnd.valid
  polynomial.io.in.bits.coeffs := coef

  if(nOrder != 0) {
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
    val sdx = Cat(~dx(inBitW-1), dx(inBitW-2, 0))
    polynomial.io.in.bits.x.get := sdx.asSInt
  }

  // -----------------------------------------------------------------
  // postproc (mul)

  // 1.0~2.0
  val yManW1 = Cat(1.U(1.W), polynomial.io.out.bits)

  // page=0
  val postPage = Pipe(rnd.valid, page, polynomialStage)
  val postAddr = Pipe(rnd.valid, addr, polynomialStage)
  val postEx   = Pipe(rnd.valid, xex , polynomialStage)
  val postMan  = Pipe(rnd.valid, xman, polynomialStage)

  // page=1
  val scaleTable = VecInit(table.scaleTable.map(_.U(spec.W.W)).toSeq)
  val scale      = scaleTable(postAddr.bits)
  val scaleEx    = scale(spec.manW+spec.exW-1, spec.manW)
  val scaleMan   = scale(spec.manW-1, 0)

  val mulEx    = Mux(postPage.bits === 0.U, (spec.exBias-2).U(spec.exW.W) - postEx.bits, scaleEx)
  val mulMan   = Mux(postPage.bits === 0.U, postMan.bits, scaleMan)
  val mulManW1 = Cat(1.U(1.W), mulMan)

  val prodW         = (1+spec.manW)*2
  val prod          = yManW1 * mulManW1
  val prodMoreThan2 = prod.head(1) === 1.U
  val prodManW1     = Mux(prodMoreThan2,
    prod(prodW-1, prodW-1-spec.manW),
    prod(prodW-2, prodW-2-spec.manW)
  )

  val zsgn = ShiftRegister(rnd.bits(rndW-1), polynomialStage + mulStage, 0.U, true.B)
  val zex  = ShiftRegister(mulEx + prodMoreThan2,   , mulStage, 0.U, true.B)
  val zman = ShiftRegister(prodManW1(spec.manW-1, 0), mulStage, 0.U, true.B)

  val z = Cat(zsgn, zex, zman)

  io.z.valid := ShiftRegister(polynomial.io.out.valid, mulStage, false.B, true.B)
  io.z.bits  := z
}
