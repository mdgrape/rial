package rial.rng

import scala.language.reflectiveCalls
import scala.math._

import chisel3._
import chisel3.util._

import rial.arith._
import rial.math._
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._

//
// HTBoxMuller:
//   io.in.x -> [i2f] -> [log] -> [sqrt] -> [mul] -> io.out.z1
//   io.in.y -> [i2f] -> [sin] -----------> [mul] -> io.out.z2
//                    -> [cos] ----------->
//
case class HTBoxMullerConfig(
  val rndW:          Int,                 // input int width
  val realSpec:      RealSpec,            // output real width
  val polySpec:      PolynomialSpec,      // polynomial eval spec
  val int2floatStge: PipelineStageConfig, // Float01OpenCloseFromOneUInt stages
  val polyPreStage:  PipelineStageConfig, // polynomial eval preproc stages
  val polyCalcStage: PipelineStageConfig, // polynomial eval func calc stages
  val polyPostStage: PipelineStageConfig, // polynomial eval postproc stages
  val i2fPolyGap:    Boolean,             // register between i2f/polynomial
  val preCalcGap:    Boolean,             // register between preproc/calc
  val tableCalcGap:  Boolean,             // register between table/calc
  val calcPostGap:   Boolean,             // register between calc/postproc
  val mulStage:      PipelineStageConfig, // multiplier config
) {
  val i2fGap = if(i2fPolyGap)  {1} else {0}
  val pcGap = if(preCalcGap)   {1} else {0}
  val tcGap = if(tableCalcGap) {1} else {0}
  val cpGap = if(calcPostGap)  {1} else {0}

  val polyTotalStage = polyPreStage.total +
    pcGap + tcGap + polyCalcStage.total +
    cpGap + polyPostStage.total
}

// ----------------------------------------------------------------------------

class HTBoxMuller(
  val cfg: HTBoxMullerConfig
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  // x -> [i2f] -> [log] -> [sqrt] -> [mul] -> z1
  // y -> [i2f] -> [sin] -----------> [mul] -> z2
  //            -> [cos] ----------->

  // no stall, io.in is always ready.
  // And io.out must always be ready.
  val io = IO(new Bundle {
    val in = new Bundle {
      val ready = Output(Bool())
      val valid = Input(Bool())
      val x     = Input(UInt(cfg.rndW.W))
      val y     = Input(UInt(cfg.rndW.W))
    }
    val out = new Bundle {
      val ready = Input(Bool())
      val valid = Output(Bool())
      val z1    = Output(UInt(cfg.realSpec.W.W))
      val z2    = Output(UInt(cfg.realSpec.W.W))
    }
  })

  when(io.out.valid) {
    assert(io.out.ready === true.B,
      "HTBoxMuller never stalls. output must always be ready.")
  }
  io.in.ready := true.B

  // --------------------------------------------------------------------------
  // convert x and y into float.

  val convX = Module(new GenRandomFloat01Close(cfg.rndW, realSpec, cfg.int2floatStge))
  val convY = Module(new GenRandomFloat01Close(cfg.rndW, realSpec, cfg.int2floatStge))

  convX.io.rnd := io.in.x
  convY.io.rnd := io.in.y

  val decomposerX = Module(new DecomposeReal(realSpec))
  val decomposerY = Module(new DecomposeReal(realSpec))

  decomposerX.io.real := ShiftRegister(convX.io.z, cfg.i2fGap, 0.U, true.B)
  decomposerY.io.real := ShiftRegister(convY.io.z, cfg.i2fGap, 0.U, true.B)

  val xdcmp = decomposerX.io.decomp
  val ydcmp = decomposerY.io.decomp

  val convStage = convX.nStage + cfg.i2fGap

  // --------------------------------------------------------------------------
  // sqrt(-2log(x))

  val logModule  = Module(new HTBoxMullerLog(cfg))
  val sqrtModule = Module(new HTBoxMullerSqrt(cfg))
  val logStage  = logModule.nStage
  val sqrtStage = sqrtModule.nStage

  logModule.io.en := ShiftRegister(io.in.valid, convStage)
  logModule.io.x  := xdcmp

  val decomposerLog = Module(new DecomposeReal(realSpec))
  decomposerLog.io.real := logModule.io.z

  sqrtModule.io.en := ShiftRegister(io.in.valid, convStage + logStage)
  sqrtModule.io.x  := decomposerLog.io.decomp

  val sqrtlogx = sqrtModule.io.z

  // --------------------------------------------------------------------------
  // sin/cos(2piy)

  val sinModule = Module(new HTBoxMullerSinCos2Pi(cfg, isSin=true))
  val cosModule = Module(new HTBoxMullerSinCos2Pi(cfg, isSin=false))
  val sinStage  = sinModule.nStage
  val cosStage  = cosModule.nStage

  sinModule.io.en := ShiftRegister(io.in.valid, convStage)
  sinModule.io.x  := ydcmp
  cosModule.io.en := ShiftRegister(io.in.valid, convStage)
  cosModule.io.x  := ydcmp

  assert(logStage + sqrtStage >= sinStage)
  val sin2piy = ShiftRegister(sinModule.io.z, logStage + sqrtStage - sinStage)
  val cos2piy = ShiftRegister(cosModule.io.z, logStage + sqrtStage - cosStage)

  // --------------------------------------------------------------------------
  // z1 = sqrt(-2log(x)) * sin(2piy)
  // z2 = sqrt(-2log(x)) * cos(2piy)

  val mul1 = Module(new MultFPGeneric(realSpec, realSpec, realSpec, RoundSpec.roundToEven, cfg.mulStage))
  val mul2 = Module(new MultFPGeneric(realSpec, realSpec, realSpec, RoundSpec.roundToEven, cfg.mulStage))

  mul1.io.x := sqrtlogx
  mul1.io.y := sin2piy

  mul2.io.x := sqrtlogx
  mul2.io.y := cos2piy

  io.out.z1 := mul1.io.z
  io.out.z2 := mul2.io.z
  io.out.valid := ShiftRegister(io.in.valid, convStage + logStage + sqrtStage + cfg.mulStage.total, false.B, true.B)

  val nStage = convStage + logStage + sqrtStage + cfg.mulStage.total
}
