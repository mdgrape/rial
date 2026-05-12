package rial.approx

import rial.util._

import chisel3._
import chisel3.util._

class PolynomialEval(
  cfg: ApproxFuncConfig,
  bitW: PolynomialBitWidth,
  stage: PipelineStageConfig,
) extends Module {
  assert(cfg.nOrder == bitW.nOrder)

  val order = cfg.nOrder
  val dxW   = cfg.inBitW

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
     val coeffs = UInt(bitW.coefBits.reduce(_+_).W)
     val x = if (order != 0) { Some(UInt(cfg.inBitW.W)) } else { None }
    }))
    val out = Valid(UInt(cfg.outPrec.W)) // [0.0, 1.0)
  })

  val out = WireDefault(0.U.asTypeOf(new Bundle {
    val valid = Bool()
    val bits  = UInt(cfg.outPrec.W)
  }))
  io.out := ShiftRegister(out, stage.total)

  if(order == 0) {
    assert(cfg.coefBits.length() == 1)
    assert(cfg.outPrec == cfg.coefBits.reduce(_+_))

    out.valid := io.in.valid
    out.bits  := io.in.bits.coeffs // if order==0, then there is only one parameter.

  } else {
    // conversion to SInt (-1 ~ +1 range) must be done outside
    val x = io.in.bits.dx.get.asSInt

    val ws     = bitW.coefBits.scanLeft(0)(_+_)
    val ranges = ws.zip(ws.tail)
    val coeffs = MixedVec(ranges.map { case (l,u) => io.in.bits.coeffs(l, u).asSInt })

    def horner( coeff: SInt, curr: SInt, x: SInt ): SInt = {
      // minimize bit width in multiply
      val currW: Int = curr.getWidth()
      val xPrec: Int = x.getWidth() - 1 // [-1, 1) = [1'0000, 01111]
      val xminW: Int = min(currW, xPrec)

      val prod = x(xPrec, xPrec-xminW).asSInt * curr
      val sum  = (prod >> xminW) +& coeff
      sum
    }
    val res = coeffs.init.foldRight(coeffs.last)( (c, z) => horner(c, z, x) ).asUInt
    val resMax = Fill(cfg.outPrec, 1.U(1.W))

    out.valid := io.in.valid
    out.bits  := Mux(res > resMax, resMax, res)
  }
}

// ------------------------------------------------------------
// software

object PolynomialEval {

  private def evalImpl(
    x: Long,
    bitWidth: PolynomialBitWidth,
    tableElem: Seq[Long],
  ): BigInt = {

    var calcW = bitWidth.calcBits(tableElem.length-1)
    var sum   = BigInt(tableElem (tableElem.length-1))

    // println(f"i=${tableElem.length-1}: inBitW=${cfg.inBitW}, coefW=${bitWidth.coefBits(tableElem.length-1)}, calcW=${calcW}, x=${x}, sum=${sum}, coef=${tableElem(tableElem.length-1)}")
    for(i <- (tableElem.length-2) to 0 by -1) {
      val xW   = min(calcW, cfg.inBitW) // minimum required bitW
      val prod = BigInt(x >> (cfg.inBitW - xW)) * sum

      //                                 -16                            15
      // [-1, 1] -> signed fixedpoint [10000,..11111, 00000, 00001,..01111]
      //
      // since x is signed and in [-1,1), the position of decimal point is
      // slightly different from its width.
      // if x is unsigned and in [0,1), decimal point is the same as its width.
      // but if x is signed and in [-1,1), decimal point locates width-1.
      sum   = (prod >> (xW-1)) + tableElem(i)
      calcW = bitWidth.calcBits(i)
      // println(f"i=${i}: inBitW=${cfg.inBitW}, coefW=${bitWidth.coefBits(i)}, calcW=${calcW}, x=${x}, prod=${prod>>xW}, sum=${sum}, coef=${tableElem(i)}")
    }
    sum
  }

  def eval(
    x: Double,
    coeffs: PolynomialCoefficients
  ): Double = {
    val idx = floor((x - cfg.xMin) * cfg.rIntervalW).toInt
    val xiv = (x - cfg.xMin - idx * cfg.intervalW) * cfg.rIntervalW * 2.0 - 1.0
    assert(-1.0 <= xiv && xiv <= 1.0)

    val xl0 = round(scalb(xiv, cfg.inPrec))
    val xlmin = -(1L << (cfg.inPrec))
    val xlmax =  (1L << (cfg.inPrec)) - 1L
    val xl = max(min(xl0, xlmax), xlmin)

    // println(f"xiv=${xiv}%8g, xl0=${xl0}, clamp(${xlmin}, ${xlmax})=${xl}, xiv-xl=${ abs(xiv - scalb(xl, -cfg.inBitW+1)) }%8g")

    val z = evalImpl(xl, coeffs.bitW, coeffs.coeffs(idx))
    scalb(z.toDouble, -cfg.outPrec)
  }
}
