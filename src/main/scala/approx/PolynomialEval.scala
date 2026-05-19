package rial.approx

import scala.math._

import rial.util._

import chisel3._
import chisel3.util._

class PolynomialEval(
  bitW: PolynomialBitWidth,
  stage: Int,
) extends Module {

  val order = bitW.nOrder

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
     val coeffs = Input(MixedVec(bitW.coefBits.map{w => SInt(w.W)}))
     val x = if (order != 0) { Some(SInt(bitW.inBitW.W)) } else { None }
    }))
    val out = Valid(UInt(bitW.outBitW.W)) // [0.0, 1.0)
  })

  val out = WireDefault(0.U(bitW.outBitW.W))
  io.out := Pipe(io.in.valid, out, stage)

  if(order == 0) {
    assert(bitW.coefBits.length == 1)
    assert(bitW.outPrec == bitW.coefBits.reduce(_+_))

    out := io.in.bits.coeffs(0) // if order==0, then there is only one parameter.

  } else {
    // conversion to SInt (-1 ~ +1 range) must be done outside
    val x      = io.in.bits.x.get
    val coeffs = io.in.bits.coeffs

    def horner( coeff: SInt, curr: SInt, x: SInt ): SInt = {
      // minimize bit width in multiply
      val currW: Int = curr.getWidth
      val xPrec: Int = x.getWidth - 1 // [-1, 1) = [1'0000, 01111]
      val xminW: Int = min(currW, xPrec)

      val prod = x(xPrec, xPrec-xminW).asSInt * curr
      val sum  = (prod >> xminW) +& coeff
      sum
    }
    val res = coeffs.init.foldRight(coeffs.last)( (c, z) => horner(c, z, x) ).asUInt
    val resMax = Fill(bitW.outBitW, 1.U(1.W))

    out := Mux(res > resMax, resMax, res)
  }
}

// ------------------------------------------------------------
// software

object PolynomialEval {

  private def evalInterval(
    x: BigInt,
    bitW: PolynomialBitWidth,
    tableElem: Seq[Long],
  ): BigInt = {

    var calcW = bitW.calcBits(tableElem.length-1)
    var sum   = BigInt(tableElem (tableElem.length-1))

    // println(f"i=${tableElem.length-1}: inBitW=${bitW.inBitW}, coefW=${bitW.coefBits(tableElem.length-1)}, calcW=${calcW}, x=${x}, sum=${sum}, coef=${tableElem(tableElem.length-1)}")
    for(i <- (tableElem.length-2) to 0 by -1) {
      val xW   = min(calcW, bitW.inBitW) // minimum required bitW
      val prod = (x >> (bitW.inBitW - xW)) * sum

      //                                 -16                            15
      // [-1, 1] -> signed fixedpoint [10000,..11111, 00000, 00001,..01111]
      //
      // since x is signed and in [-1,1), the position of decimal point is
      // slightly different from its width.
      // if x is unsigned and in [0,1), decimal point is the same as its width.
      // but if x is signed and in [-1,1), decimal point locates width-1.
      sum   = (prod >> (xW-1)) + tableElem(i)
      calcW = bitW.calcBits(i)
      // println(f"i=${i}: inBitW=${bitW.inBitW}, coefW=${bitW.coefBits(i)}, calcW=${calcW}, x=${x}, prod=${prod>>xW}, sum=${sum}, coef=${tableElem(i)}")
    }
    sum
  }

  def eval(
    addr: Int,
    x: BigInt, // unsigned
    coeffs: PolynomialCoefficients,
  ): BigInt = {
    if(coeffs.bitW.nOrder == 0) {
      // just a table, no calc. no x.
      BigInt(coeffs.coeffs(addr)(0))
    } else {
      // convert unsigend [0, 1) -> signed [-1, 1)
      // 0000,..0111,1000,..1111 -> 1000,..1111,0000,..0111
      // then, (normally, on a hardware implementation) flip the MSB!
      // on a software impl, subtract `1000`.
      val sx = x - (BigInt(1) << coeffs.bitW.inPrec)
      PolynomialEval.evalInterval(sx, coeffs.bitW, coeffs.coeffs(addr))
    }
  }
}
