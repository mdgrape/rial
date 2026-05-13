package rial.funcs

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
    val real   = Input(UInt(spec.W.W))
    val decomp = Output(new DecomposedReal(spec))
  })

  val (sgn,  ex,  man) = FloatChiselUtil.decompose (spec, io.real)
  val (zero, inf, nan) = FloatChiselUtil.checkValue(spec, io.real)
  io.decomp.sgn  := sgn
  io.decomp.ex   := ex
  io.decomp.man  := man
  io.decomp.zero := zero
  io.decomp.inf  := inf
  io.decomp.nan  := nan
}


