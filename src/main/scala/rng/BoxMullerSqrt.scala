package rial.rng

import scala.language.reflectiveCalls
import scala.math._

import chisel3._
import chisel3.util._

import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.arith.RealSpec
import rial.arith.FloatChiselUtil

import rial.math._

// -------------------------------------------------------------------------

class HTBoxMullerSqrt(
  val cfg: HTBoxMullerConfig
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val adrW  = polySpec.adrW
  val fracW = polySpec.fracW
  val dxW   = polySpec.dxW
  val order = polySpec.order

  val cbit  = SqrtTableCoeff.getCBits(realSpec, polySpec)

  val preStage  = cfg.polyPreStage
  val calcStage = cfg.polyCalcStage
  val postStage = cfg.polyPostStage
  val pcGap: Int = if(cfg.preCalcGap)   { 1 } else { 0 }
  val tcGap: Int = if(cfg.tableCalcGap) { 1 } else { 0 }
  val cpGap: Int = if(cfg.calcPostGap)  { 1 } else { 0 }

  def nStage(): Int = {
    preStage.total + pcGap + tcGap + calcStage.total + cpGap + postStage.total
  }

  // ---------------------------------------------------------------------------

  val io = IO(new Bundle {
    val en = Input  (Bool())
    val x  = Flipped(new DecomposedRealOutput(realSpec))
    val z  = Output (UInt(realSpec.W.W))
  })

  // ---------------------------------------------------------------------------

  val preproc = Module(new SqrtPreProcess(realSpec, polySpec, preStage))

  preproc.io.en     := io.en
  preproc.io.x.sgn  := io.x.sgn
  preproc.io.x.ex   := io.x.ex
  preproc.io.x.man  := io.x.man
  preproc.io.x.zero := io.x.zero
  preproc.io.x.inf  := io.x.inf
  preproc.io.x.nan  := io.x.nan

  // ---------------------------------------------------------------------------

  val table = Module(new SqrtTableCoeff(realSpec, polySpec, cbit))

  table.io.en  := ShiftRegister(io.en, preproc.nStage + pcGap, false.B, true.B)
  table.io.adr := ShiftRegister(preproc.io.adr, pcGap)

  // ---------------------------------------------------------------------------

  val eval = Module(new PolynomialEval(realSpec, polySpec, cbit, calcStage))

  val polynomialCoeff = ShiftRegister(table.io.cs.asUInt, tcGap)
  eval.io.coeffs := polynomialCoeff.asTypeOf(new TableCoeffInput(cbit))
  if(order != 0) {
    val dx = preproc.io.dx.get
    eval.io.dx.get := ShiftRegister(dx, pcGap + tcGap)
  }

  // ---------------------------------------------------------------------------

  val otherStage = PipelineStageConfig.atOut(preStage.total + pcGap + tcGap + calcStage.total)
  val other = Module(new SqrtOtherPath(realSpec, polySpec, otherStage))

  other.io.x.sgn  := io.x.sgn
  other.io.x.ex   := io.x.ex
  other.io.x.man  := io.x.man
  other.io.x.zero := io.x.zero
  other.io.x.inf  := io.x.inf
  other.io.x.nan  := io.x.nan

  // ---------------------------------------------------------------------------

  val postproc = Module(new RoundingPostProcess(realSpec, polySpec, postStage))

  postproc.io.en     := ShiftRegister(io.en, preStage.total + pcGap + tcGap + calcStage.total)
  postproc.io.zres   := ShiftRegister(eval.io.result,  cpGap)
  postproc.io.zother := ShiftRegister(other.io.zother, cpGap)

  io.z := postproc.io.z
}

