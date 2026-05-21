package rial.math

import chisel3._
import chisel3.util._

import rial.arith.{RealSpec, FloatChiselUtil}

private[rial] class DecomposedReal(
  val spec: RealSpec
) extends Bundle {
  val sgn  = UInt(1.W)
  val ex   = UInt(spec.exW.W)
  val man  = UInt(spec.manW.W)
  val zero = Bool()
  val inf  = Bool()
  val nan  = Bool()
}

private[rial] class DecomposeReal(
  val spec: RealSpec
) extends Module {

  val io = IO(new Bundle {
    val real   = Flipped(Valid(UInt(spec.W.W)))
    val decomp = Valid(new DecomposedReal(spec))
  })

  io.decomp.valid := io.real.valid

  val (sgn,  ex,  man) = FloatChiselUtil.decompose (spec, io.real.bits)
  val (zero, inf, nan) = FloatChiselUtil.checkValue(spec, io.real.bits)
  io.decomp.bits.sgn  := sgn
  io.decomp.bits.ex   := ex
  io.decomp.bits.man  := man
  io.decomp.bits.zero := zero
  io.decomp.bits.inf  := inf
  io.decomp.bits.nan  := nan
}


