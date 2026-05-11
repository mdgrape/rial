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
    val out = Valid(UInt(cfg.outBitW.W))
  })

  val out = WireDefault(0.U.asTypeOf(new Bundle {
    val valid = Bool()
    val bits  = UInt(cfg.outBitW.W)
  }))
  io.out := ShiftRegister(out, stage.total)

  if(order == 0) {
    assert(cfg.coefBits.length() == 1)
    assert(cfg.outBitW == cfg.coefBits.reduce(_+_))

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
      val sum  = (prod >> xminW) + coeff
      sum
    }
    val res = coeffs.init.foldRight(coeffs.last)( (c, z) => horner(c, z, x) )
    out.valid := io.in.valid
    out.bits  := res.asUInt
  }
}
