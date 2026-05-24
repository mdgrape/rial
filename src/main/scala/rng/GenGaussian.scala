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

//
// inverse CDF of gaussian
//
private[rial] class ProbitTable(
  val spec:       RealSpec,
  val nOrder:     Int,
  val nIntervals: Int,
  val inBitW:     Int,
  val outBitW:    Int,
) {
  private val sqrt2 = scala.math.sqrt(2)
  private val log2 = (a: Double) => log(a) / log(2.0)

  val probitHalf = (x: Double) => { sqrt2 * Erf.erfInv(x) }

  val scaleTable = ArrayBuffer.empty[BigInt]

  val ft = new FuncTable(nOrder)

  // probitHalf(0~0.5) -> 1.2 ~ 1.4
  println("pribitHalf(x)/x-1.0")
  ft.addInterval((x: Double) => { probitHalf(x) / x - 1.0 }, 0.0, 0.5)
  scaleTable += 0.toBigInt

  for(i <- 1 until nIntervals) {
    val xMin  = 1.0 - pow(2.0, -i)
    val xMax  = 1.0 - pow(2.0, -(i+1))
    val fMin  = probitHalf(xMin)
    val fMax  = probitHalf(xMax)
    assert(fMax / fMin < 2.0)

    val f = (x: Double) => { probitHalf(x) / fMin - 1.0 }
    ft.addInterval(f, xMin, xMax)
    scaleTable += (new RealGeneric(spec, fMin)).value.toBigInt

    // println(f"pribitHalf(x) for ${xMin} to ${xMax}: ${fMin}, ${fMax}")
    // println(f"pribitHalf(x)/fMin-1.0 for ${xMin} to ${xMax}: ${f(xMin)}, ${f(xMax)}")
  }
  assert(ft.intervals.length == nIntervals)

  val maxValue = new RealGeneric(spec, probitHalf(1.0 - pow(2.0, -(nIntervals+1))))

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

  // 0bxxxxxx1 -> 0          : 1/ 4~1/2 : exBias-2
  // 0bxxxxx10 -> 1          : 1/ 8~1/4 : exBias-3
  // 0bxxxx100 -> 2          : 1/16~1/8 : exBias-4
  // ...
  // 0b1000000 -> getWidth-1 :
  // 0b0000000 -> ?
  val xex = Mux(xex0 === 0.U, xex0.getWidth.U, PriorityEncoder(xex0))

  // ---------------------------------------------------------------

  val inBitW  = spec.manW
  val outBitW = spec.manW

  val nIntervals = xex0.getWidth + 1
  val table  = new ProbitTable(spec, nOrder, nIntervals, inBitW, outBitW)
  val bitW   = table.getBitWidth()
  val coeffs = table.generateTable(bitW)

  val packed = VecInit(PolynomialGenerator.packTable(coeffs).map(c => {
    c.U(coeffs.bitW.coefBits.reduce(_+_).W)
  }))

  val addrRaw = Mux(page === 0.U, 0.U, xex +& 1.U)
  val addr    = Mux(addrRaw >= nIntervals.U, (nIntervals - 1).U, addrRaw)
  val coef = packed(addr).asTypeOf( MixedVec(bitW.coefBits.map(w => SInt(w.W))) )

  val polynomial = Module(new PolynomialEval(bitW, polynomialStage))

  polynomial.io.in.valid       := rnd.valid
  polynomial.io.in.bits.coeffs := coef

  if(nOrder != 0) {
    val dx = WireDefault(0.U(inBitW.W))
    when(page === 0.U) {
      dx := (Cat(1.U(1.W), xman(spec.manW-1, 1)) >> xex)
    }.otherwise {
      dx := xman
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
  val zex  = mulEx + prodMoreThan2
  val zman = prodManW1(spec.manW-1, 0)

  val z = Cat(zsgn, zex, zman)

  io.z.valid := ShiftRegister(polynomial.io.out.valid, mulStage, false.B, true.B)
  io.z.bits  := z
}
