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
    // +1 because coefficient is treated as an SInt.
    assert(bitW.outBitW+1 == bitW.coefBits.reduce(_+_),
      f"outBitW=${bitW.outBitW}, coefBits=${bitW.coefBits.mkString("[",",","]")}, "+
      f"sum=${bitW.coefBits.reduce(_+_)}")
    assert(io.in.bits.coeffs(0).asUInt(bitW.outBitW) === 0.U,
      "if nOrder==0, coefficient must be in [0, 1) and the MSB must be zero")

    out := io.in.bits.coeffs(0).asUInt(bitW.outBitW-1, 0) // if order==0, then there is only one parameter.

  } else {
    // conversion to SInt (-1 ~ +1 range) must be done outside
    val x      = io.in.bits.x.get
    val coeffs = io.in.bits.coeffs
    val xPrec  = bitW.inPrec

    def mulByX(curr: SInt): SInt = {
      val xminW = min(curr.getWidth, xPrec) // min(calcBits, xPrec)
      val prod  = x(xPrec, xPrec - xminW).asSInt * curr
      (prod >> xminW).asSInt
    }

    val accVec = Wire(MixedVec((1 to order).map(i => SInt(bitW.calcBits(i).W))))

    // since 0-th (the last) order is `res`, we need to do -1
    accVec(order-1) := coeffs(order)

    for(i <- (order - 1) to 1 by -1) {
      val sumFull = mulByX(accVec(i)) +& coeffs(i)
      accVec(i-1) := sumFull

      assert(accVec(i-1) === sumFull, "PolynomialEval: horner intermediate (order=" +
        i + ") overflowed calcBits=" + bitW.calcBits(i))
    }

    // final stage (order 0): result may legitimately reach 1.0, saturate.
    val res    = (mulByX(accVec(0)) +& coeffs(0)).asUInt
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
      val xminW = min(calcW, bitW.inBitW - 1)
      val prod  = (x >> (bitW.inBitW - 1 - xminW)) * sum

      //                                 -16                            15
      // [-1, 1] -> signed fixedpoint [10000,..11111, 00000, 00001,..01111]
      //
      // since x is signed and in [-1,1), the position of decimal point is
      // slightly different from its width.
      // if x is unsigned and in [0,1), decimal point is the same as its width.
      // but if x is signed and in [-1,1), decimal point locates width-1.
      sum   = (prod >> xminW) + tableElem(i)
      calcW = bitW.calcBits(i)
      // println(f"i=${i}: inBitW=${bitW.inBitW}, coefW=${bitW.coefBits(i)}, calcW=${calcW}, x=${x}, prod=${prod>>xminW}, sum=${sum}, coef=${tableElem(i)}")
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
