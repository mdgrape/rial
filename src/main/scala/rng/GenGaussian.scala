//% @file GenGaussian.scala
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
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
// generate Gaussian using Box-Muller.
//
// z1 = sqrt(-2log(x)) cos(2pi*y)
// z2 = sqrt(-2log(x)) sin(2pi*y)
//

// -----------------------------------------------------------------------------
// Fixed -> FP sin/cos(2pi*y)
//
// sin/cos are used to determine the argument (not an absolute value), so y can
// be a fixed point.
// By doing this, the resolution (the difference between the current and the
// next value) around 2piy ~ 0, pi, 2pi will be almost the same. It is not a
// problem.
//

class BoxMullerSinCos2PiPreProcessOutput(val rndW: Int, val spec: RealSpec) extends Bundle {
  val zone     = Output(Bool())
  val zsgn     = Output(UInt(1.W))
  val zzero    = Output(Bool())
  val xShift   = Output(UInt(log2Up(rndW-2).W))
  val xAligned = Output(UInt((rndW-2).W))
}

object BoxMullerSinCos2PiTableCoeff {
  def genTable(
    polySpec: PolynomialSpec
  ): FuncTableInt = {

    // --------------------------------------------------------------------------
    // in a range x in [0, 1/4), 4 < sin(2pix)/x <= 2pi.
    // So 1/2 < sin(2pix)/8x <= pi/4 = 0.78.. < 1.
    val order  = polySpec.order
    val fracW  = polySpec.fracW
    val adrW   = polySpec.adrW

    val tableD = new FuncTableDouble( x0 => {
      val x = x0 / 4.0 // convert [0, 1) to the input range, [0, 1/4)
      sin(x * Pi * 2.0) / (8.0*x)
    }, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  def getCBits(
    polySpec: PolynomialSpec
  ): Seq[Int] = {
    val tableI = genTable(polySpec)
    tableI.cbit
  }
}

class BoxMullerSinCos2PiPreProc(
    rndW: Int,      // width of input fixedpoint
    spec: RealSpec, // output width
    polySpec: PolynomialSpec,
    cbit: Seq[Int],
    stage: PipelineStageConfig
  ) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order  = polySpec.order
  val fracW  = polySpec.fracW
  val adrW   = polySpec.adrW
  val dxW    = polySpec.dxW

  assert(rndW > 2 + adrW + dxW) // manW + 2; satisfied by (i32, f32) or (i64, f64)

  // --------------------------------------------------------------------------

  val io = IO(new Bundle {
    val en    = Input(UInt(1.W))
    val isSin = Input(Bool())
    val rnd   = Input(UInt(rndW.W))

    val cs    = Flipped(new TableCoeffInput(cbit))
    val dx    = if (order != 0) { Some(Output(UInt(dxW.W))) } else { None }
    val out   = new BoxMullerSinCos2PiPreProcessOutput(rndW, spec)
  })

  val rnd = enable(io.en, io.rnd)

  // detect special values
  val xExactQuot = rnd(rndW-3, 0) === 0.U
  val xHighIs0   = rnd(rndW-1, rndW-2) === 0.U
  val xHighIs1   = rnd(rndW-1, rndW-2) === 1.U
  val xHighIs2   = rnd(rndW-1, rndW-2) === 2.U
  val xHighIs3   = rnd(rndW-1, rndW-2) === 3.U

  // convert x in [0, 1) to [0, 1/4).
  // Since we convert x range, in case of x is exactly N/4, the resulting x
  // will become indistinguishable from 0. but sin(0) and sin(2pi/4) completely
  // differ to each other. We need to check if x is exactly 1/4, 2/4, or 3/4.
  val zone = Wire(Bool())
  val zsgn = Wire(UInt(1.W))
  val x = Wire(UInt((rndW-2).W))
  when(io.isSin) {

    x := Mux(rnd(rndW-2) === 1.U, ~(rnd(rndW-3, 0))+1.U, rnd(rndW-3, 0))

    // x is exactly 1/4 or 3/4
    zone := (xHighIs1 || xHighIs3) && xExactQuot
    // if x is exactly 0 or 1/2, then the result is exactly zero, not -0.
    zsgn := rnd(rndW-1) & ~((xHighIs0 || xHighIs2) && xExactQuot).asUInt

  } otherwise {

    x := Mux(rnd(rndW-2) === 0.U, ~(rnd(rndW-3, 0))+1.U, rnd(rndW-3, 0))

    // x is exactly 0 or 1/2
    zone := (xHighIs0 || xHighIs2) && xExactQuot
    // if x is exactly 1/4 or 3/4, then the result is exactly zero, not -0.
    zsgn := (rnd(rndW-1) ^ rnd(rndW-2)) &
           ~((xHighIs1 || xHighIs3) && xExactQuot).asUInt
  }

  if(order != 0) {
    val dx = Cat(~x(x.getWidth-adrW-1), x(x.getWidth-adrW-2, x.getWidth-adrW-dxW))
    io.dx.get := ShiftRegister(enable(io.en, dx), nStage)
  }
  io.out.zsgn := ShiftRegister(zsgn, nStage)
  io.out.zone := ShiftRegister(zone, nStage)

  // --------------------------------------------------------------------------

  val adr = x(x.getWidth-1, x.getWidth-adrW)

  val tableI = BoxMullerSinCos2PiTableCoeff.genTable(polySpec)
  val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
  val coeff = getSlices(coeffTable(adr), coeffWidth)

  val coeffs = Wire(new TableCoeffInput(cbit))
  for (i <- 0 to order) {
    coeffs.cs(i) := coeff(i)
  }

  io.cs := ShiftRegister(enable(io.en, coeffs), nStage)

  // --------------------------------------------------------------------------

  val zzero  = x === 0.U
  val xShift = Mux(zzero, 0.U, PriorityEncoder(Reverse(x)))
  val xAligned = (x << xShift)(x.getWidth-1, 0)
  assert(xAligned(x.getWidth-1) === 1.U || zzero || io.en =/= 1.U)

  io.out.zzero    := ShiftRegister(zzero,    nStage)
  io.out.xShift   := ShiftRegister(xShift,   nStage)
  io.out.xAligned := ShiftRegister(xAligned, nStage)
}

class BoxMullerSinCos2PiPostProc(
    rndW: Int,      // width of input fixedpoint
    spec: RealSpec, // output width
    polySpec: PolynomialSpec,
    stage: PipelineStageConfig
  ) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order  = polySpec.order
  val fracW  = polySpec.fracW
  val adrW   = polySpec.adrW
  val dxW    = polySpec.dxW

  assert(rndW > 2 + adrW + dxW) // manW + 2; satisfied by (i32, f32) or (i64, f64)

  val io = IO(new Bundle {
    val en   = Input(Bool())
    val pre  = Flipped(new BoxMullerSinCos2PiPreProcessOutput(rndW, spec))
    val zres = Input(UInt(fracW.W))
    val z    = Output(UInt(spec.W.W))
  })

  val zone     = io.pre.zone
  val zsgn     = io.pre.zsgn
  val zzero    = io.pre.zzero
  val xShift   = io.pre.xShift
  val xAligned = enable(io.en, io.pre.xAligned)
  val zres     = enable(io.en, io.zres)

  val zProd      = xAligned * zres
  val zProdW     = xAligned.getWidth + fracW

  val zMoreThan2 = zProd(zProdW - 1)
  val zShifted   = Mux(zMoreThan2, zProd(zProdW-2, zProdW-manW-1),
                                   zProd(zProdW-3, zProdW-manW-2))

  val zRounded   = zShifted +& Mux(zMoreThan2, zProd(zProdW-manW-2),
                                               zProd(zProdW-manW-3))
  val zMoreThan2AfterRound = zRounded(manW)
  val zmanW1     = Mux(zMoreThan2AfterRound, Cat(1.U(1.W), 0.U(manW.W)),
                                             Cat(1.U(1.W), zRounded(manW-1, 0)))
  val zexInc     = zMoreThan2 + zMoreThan2AfterRound

  val zex  = Mux(zone, exBias.U(exW.W), Mux(zzero, 0.U, (exBias-1).U(exW.W) - xShift + zexInc))
  val zman = Mux(zone || zzero, 0.U, zmanW1(manW-1, 0))

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := ShiftRegister(enable(io.en, z), nStage)
}


// ===========================================================================
//
// z = -2log(1-x): Fixed -> FP
//
// x = [0, 1).
//

object BoxMullerLogTableCoeff {
  def genTable(
    polySpec: PolynomialSpec
  ): FuncTableInt = {

    val log2 = (x:Double) => {log(x) / log(2.0)}

    // --------------------------------------------------------------------------
    val order  = polySpec.order
    val fracW  = polySpec.fracW
    val adrW   = polySpec.adrW

    val tableD = new FuncTableDouble( x0 => {
      val x = x0 + 1.0
      log2(x)
    }, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  def getCBits(
    polySpec: PolynomialSpec
  ): Seq[Int] = {
    val tableI = genTable(polySpec)
    tableI.cbit
  }
}

class BoxMullerLogPreProcessOutput(val rndW: Int, val spec: RealSpec) extends Bundle {
  val zzero  = Output(Bool())
  val xShift = Output(UInt(log2Up(rndW-2).W))
}

//
// preproc takes FixedPoint x in 0x0000 ~ 0xFFFF and normalize (1 - x).
// let y = 1-x. y = 2^yex * 1.yman.
// table calculates log2(1.yman).
//
class BoxMullerLogPreProc(
  rndW: Int,      // width of input fixedpoint
  spec: RealSpec, // output width
  polySpec: PolynomialSpec,
  cbit: Seq[Int],
  stage: PipelineStageConfig
  ) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order  = polySpec.order
  val fracW  = polySpec.fracW
  val extraBits = fracW - manW
  val adrW   = polySpec.adrW
  val dxW    = polySpec.dxW

  val io = IO(new Bundle {
    val en = Input(UInt(1.W))
    val x  = Input(UInt(rndW.W))
    val cs = Flipped(new TableCoeffInput(cbit))
    val dx = if (order != 0) { Some(Output(UInt(dxW.W))) } else { None }
    val out = new BoxMullerLogPreProcessOutput(rndW, spec)
  })

  val negx = enable(io.en, ~io.x + 1.U)  // 1 - x
  val zzero = negx === 0.U // if 1 - x == 0, then x == 1. log(1) = 0.

  val nxShift   = Mux(zzero, 0.U, PriorityEncoder(Reverse(negx)))
  val nxShifted = (negx << nxShift)(rndW-1, 0)
  assert(zzero || nxShifted(rndW-1) === 1.U || io.en =/= 1.U)

  // ---------------------------------------------------------------------------
  // special values
  io.out.zzero  := ShiftRegister(zzero,   nStage)
  io.out.xShift := ShiftRegister(nxShift, nStage)

  // ---------------------------------------------------------------------------
  // table coeffs

  val xman = nxShifted(rndW-2, 0) // remove "hidden" bit to normalize
  val xmanW = xman.getWidth

  val adr = xman(xmanW-1, xmanW-adrW)

  val tableI = BoxMullerLogTableCoeff.genTable(polySpec)
  val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
  val coeff = getSlices(coeffTable(adr), coeffWidth)

  val coeffs = Wire(new TableCoeffInput(cbit))
  for (i <- 0 to order) {
    coeffs.cs(i) := coeff(i)
  }
  io.cs := ShiftRegister(enable(io.en, coeffs), nStage)

  if(order != 0) {
    val dx = Cat(~xman(xmanW-adrW-1), xman(xmanW-adrW-2, xmanW-adrW-dxW))
    io.dx.get := ShiftRegister(enable(io.en, dx), nStage)
  }
}

// BoxMullerLog calculates -2log(1-x). multiplying 2 is just adding 1 to ex.
//
// table calculates log2(1.yman).
// log2(y) = log2(2^yex * 1.yman) = yex + log2(1.yman). yex < 0.
// log_e(y) = log_e(2) * log2(y).
// yex is always negative and log2(1.man) is always positive.
// log2(1.man) < 1 and yex <= -1. So, it calculates (|yex| - log2(1.yman)) * log_e(2)
//
// yex is an integer and 0 <= log2(1.man) < 1. So subtraction requires only a
// negation and +1.
// `xShift` represents the shift required to make the MSB of 1-x to 1.
// `xShift===0` means the MSB was already 1, in other words, 0.5 < 1-x < 1.
// 0.5 < 1-x < 1 means (1-x).ex = -1. yex = -xShift-1.
//
// (|yex| - log2(1.yman)) * log_e(2) = (xShift + 1 - log2(1.yman)) * log_e(2).
//
class BoxMullerLogPostProc(
    rndW: Int,      // width of input fixedpoint
    spec: RealSpec, // output width
    polySpec: PolynomialSpec,
    stage: PipelineStageConfig
  ) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order  = polySpec.order
  val fracW  = polySpec.fracW
  val adrW   = polySpec.adrW
  val dxW    = polySpec.dxW

  assert(rndW > adrW + dxW)

  val io = IO(new Bundle {
    val en   = Input(Bool())
    val pre  = Flipped(new BoxMullerLogPreProcessOutput(rndW, spec))
    val zres = Input(UInt(fracW.W))
    val z    = Output(UInt(spec.W.W))
  })

  val zFrac = ~io.zres + 1.U

  val zzero = io.pre.zzero
  val xex   = io.pre.xShift +& (zFrac === 0.U).asUInt

  val log2x0     = Cat(xex, zFrac)
  val log2xShift = PriorityEncoder(Reverse(log2x0)) // xex might be zero
  val log2x      = enable(io.en, (log2x0 << log2xShift)(log2x0.getWidth-1, 0))
  assert(log2x(log2x0.getWidth-1) === 1.U || io.en =/= 1.U)

  val ln2 = new RealGeneric(spec, log(2.0))
  val ln2manW1 = ln2.manW1.toBigInt.U((manW+1).W)

  val zProd = log2x * ln2manW1
  val zProdW = log2x.getWidth + ln2manW1.getWidth

  val zMoreThan2 = zProd(zProdW - 1)
  val zShifted   = Mux(zMoreThan2, zProd(zProdW-2, zProdW-manW-1),
                                   zProd(zProdW-3, zProdW-manW-2))

  val zRounded   = zShifted +& Mux(zMoreThan2, zProd(zProdW-manW-2),
                                               zProd(zProdW-manW-3))
  val zMoreThan2AfterRound = zRounded(manW)
  val zmanW1     = Mux(zMoreThan2AfterRound, Cat(1.U(1.W), 0.U(manW.W)),
                                             Cat(1.U(1.W), zRounded(manW-1, 0)))
  val zexInc     = zMoreThan2 + zMoreThan2AfterRound

  val zex0 = (ln2.ex + xex.getWidth).U(exW.W) - log2xShift

  val zex  = Mux(zzero, 0.U, zex0 + zexInc)
  val zman = Mux(zzero, 0.U, zmanW1(manW-1, 0))
  val zsgn = 0.U(1.W)

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := ShiftRegister(enable(io.en, z), nStage)
}

// ============================================================================
//
// calc sqrt: FP -> FP. almost the same as math/sqrt.
//

class BoxMullerSqrtPreProcessOutput(val rndW: Int, val spec: RealSpec) extends Bundle {
  val zex = Output(UInt(spec.exW.W)) // if zex == 0, then the result is zero.
}

class BoxMullerSqrtPreProc(
    rndW: Int,      // width of input fixedpoint
    spec: RealSpec, // output width
    polySpec: PolynomialSpec,
    cbit: Seq[Int],
    stage: PipelineStageConfig
  ) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order  = polySpec.order
  val fracW  = polySpec.fracW
  val adrW   = polySpec.adrW
  val dxW    = Seq(polySpec.dxW, manW-adrW).min // now the input is FP.

  // --------------------------------------------------------------------------

  val io = IO(new Bundle {
    val en  = Input(UInt(1.W))
    val x   = Input(UInt(spec.W.W))
    val cs  = Flipped(new TableCoeffInput(cbit))
    val dx  = if (order != 0) { Some(Output(UInt(dxW.W))) } else { None }
    val out = new BoxMullerSqrtPreProcessOutput(rndW, spec)
  })

  // --------------------------------------------------------------------------

  val xexNobias = io.x(exW+manW-1, manW) - exBias.U
  val zexNobias = Cat(xexNobias(exW-1), xexNobias(exW-1, 1))
  val zex = zexNobias + exBias.U
  io.out.zex := ShiftRegister(zex, nStage)

  // --------------------------------------------------------------------------

  val adr = io.x(manW, manW-adrW)

  val tableI = SqrtSim.sqrtTableGeneration( order, adrW, manW, fracW )
  val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
  val coeff  = getSlices(coeffTable(adr), coeffWidth)

  val coeffs = Wire(new TableCoeffInput(cbit))
  for (i <- 0 to order) {
    coeffs.cs(i) := coeff(i)
  }
  io.cs := ShiftRegister(enable(io.en, coeffs), nStage)

  if(order != 0) {
    val dx = Cat(~io.x(manW-adrW-1), io.x(manW-adrW-2, manW-adrW-dxW))
    io.dx.get := ShiftRegister(enable(io.en, dx), nStage)
  }
}

class BoxMullerSqrtPostProc(
    rndW: Int,      // width of input fixedpoint
    spec: RealSpec, // output width
    polySpec: PolynomialSpec,
    stage: PipelineStageConfig
  ) extends Module {

  val nStage = stage.total

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val order  = polySpec.order
  val fracW  = polySpec.fracW
  val adrW   = polySpec.adrW
  val dxW    = polySpec.dxW

  val io = IO(new Bundle {
    val en   = Input(Bool())
    val pre  = Flipped(new BoxMullerSqrtPreProcessOutput(rndW, spec))
    val zres = Input(UInt(fracW.W))
    val z    = Output(UInt(spec.W.W))
  })

  val zsgn  = 0.U(1.W)
  val zex   = io.pre.zex
  val zzero = zex === 0.U

  val zman = Wire(UInt(manW.W))
  if(fracW == manW) {
    zman := io.zres
  } else {
    val extraBits = fracW - manW
    val zman0 = dropLSB(extraBits, io.zres) +& io.zres(extraBits-1)
    val polynomialOvf = zman0(manW)
    zman := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW-1,0))
  }

  val z = Mux(zzero, 0.U(spec.W.W), Cat(zsgn, zex, zman))
  io.z := ShiftRegister(enable(io.en, z), nStage)
}

// ============================================================================
//
// do Box-Muller.
//
// If the values are not used, it stops to reduce power consumption.
// It means that the random number queue should have a buffer having at least
// the same number of elements as the latency to keep the queue filled.
//

class BoxMuller(
    rndW: Int,      // width of input fixedpoint
    spec: RealSpec, // output width
    polySpec: PolynomialSpec,
    roundSpec: RoundSpec,
  ) extends Module {

  val order  = polySpec.order

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val io = IO(new Bundle {
    val consume = Output(Bool()) // it consumes the current random number input
    val valid   = Output(Bool()) // the z is valid gaussian random number output

    val x = Input(UInt(rndW.W))    // uniformly distributed random bits
    val z = Output(UInt(spec.W.W)) // gaussian-distribution in FP
  })

  val maxCbit = Seq(
    BoxMullerSinCos2PiTableCoeff.getCBits(polySpec),
    BoxMullerLogTableCoeff.getCBits(polySpec),
    SqrtTableCoeff.getCBits(spec, polySpec)
    ).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  // TODO:
  // - consider adding stages for each module

  //           pre   poly  post  multiply        | consume valid
  // --------------+ sqrt  cos                   | -------------
  // 0. rnd -> log +-----+ sqrt                  |    T      F
  // 1. rnd -> sin   log +-----+                 |    T      F
  // 2.     +> cos   sin   log | sqrt * sin -> z |    F      T
  // 3.        sqrt  cos   sin | sqrt * cos -> z |    F      T
  // --------------+ sqrt  cos +-----------------|---------------
  // 4. rnd -> log +-----+ sqrt                  |    T      F
  // 5. rnd -> sin   log +-----+                 |    T      F
  // 6.        cos   sin   log | sqrt * sin -> z |    F      T
  // 7.        sqrt  cos   sin | sqrt * cos -> z |    F      T
  // 8.              sqrt  cos +-----------------|---------------
  //

  val pc = RegInit(0.U(2.W)) // [log, sin, cos, sqrt]
  pc := pc + 1.U

  val rndReg = RegInit(0.U(rndW.W)) // to save the same rnd for sin and cos
  when(pc === 1.U) {
    rndReg := io.x
  }

  io.consume := pc(1) === 0.U
  io.valid   := pc(1) === 1.U

  // ---------------------------------------------------------------------------

  val sincosPre  = Module(new BoxMullerSinCos2PiPreProc (rndW, spec, polySpec, maxCbit, PipelineStageConfig.none))
  val sincosPost = Module(new BoxMullerSinCos2PiPostProc(rndW, spec, polySpec,          PipelineStageConfig.none))
  val logPre     = Module(new BoxMullerLogPreProc       (rndW, spec, polySpec, maxCbit, PipelineStageConfig.none))
  val logPost    = Module(new BoxMullerLogPostProc      (rndW, spec, polySpec,          PipelineStageConfig.none))
  val sqrtPre    = Module(new BoxMullerSqrtPreProc      (rndW, spec, polySpec, maxCbit, PipelineStageConfig.none))
  val sqrtPost   = Module(new BoxMullerSqrtPostProc     (rndW, spec, polySpec,          PipelineStageConfig.none))

  sincosPre.io.en    := (pc(0) ^ pc(1)) // pc == 1 or 2
  sincosPre.io.isSin := pc === 1.U
  sincosPre.io.rnd   := Mux(pc === 2.U, rndReg, io.x)

  logPre.io.en := pc === 0.U
  logPre.io.x  := io.x

  sqrtPre.io.en := pc === 3.U
  sqrtPre.io.x  := ShiftRegister(logPost.io.z, 1)

  // ---------------------------------------------------------------------------

  val polynomialEval = Module(new PolynomialEval(spec, polySpec, maxCbit, PipelineStageConfig.none))

  if(order != 0) {
    val polynomialDx = sincosPre.io.dx.get |
                          logPre.io.dx.get |
                         sqrtPre.io.dx.get
    polynomialEval.io.dx.get := ShiftRegister(polynomialDx, 1)
  }

  // table is accessed combinationally. There is no delay.
  val polynomialCoef = ShiftRegister(sincosPre.io.cs.asUInt |
                                       sqrtPre.io.cs.asUInt |
                                        logPre.io.cs.asUInt, 1)

  polynomialEval.io.coeffs := polynomialCoef.asTypeOf(new TableCoeffInput(maxCbit))

  val polynomialResult = ShiftRegister(polynomialEval.io.result, 1)

  // ---------------------------------------------------------------------------

  val sinResult  = RegInit(0.U(spec.W.W))
  val cosResult  = RegInit(0.U(spec.W.W))
  val logResult  = RegInit(0.U(spec.W.W))
  val sqrtResult = RegInit(0.U(spec.W.W))

  sincosPost.io.en   := ~(pc(0) ^ pc(1)) // pc == 3 or 0
  sincosPost.io.pre  := ShiftRegister(sincosPre.io.out, 2)
  sincosPost.io.zres := polynomialResult

  when(pc === 3.U) {
    sinResult := sincosPost.io.z
  }.elsewhen(pc === 0.U) {
    cosResult := sincosPost.io.z
  }

  logPost.io.en   := pc === 2.U
  logPost.io.pre  := ShiftRegister(logPre.io.out, 2)
  logPost.io.zres := polynomialResult

  when(pc === 2.U) {
    logResult := logPost.io.z
  }

  sqrtPost.io.en   := pc === 1.U
  sqrtPost.io.pre  := ShiftRegister(sqrtPre.io.out, 2)
  sqrtPost.io.zres := polynomialResult

  when(pc === 1.U) {
    sqrtResult := sqrtPost.io.z
  }

  // ---------------------------------------------------------------------------

  val a = Mux(pc === 2.U, sinResult,
          Mux(pc === 3.U, cosResult,
                          0.U))
  val b = enable(pc(1), sqrtResult)

  val amanW1 = Cat(1.U(1.W), a(manW-1, 0))
  val bmanW1 = Cat(1.U(1.W), b(manW-1, 0))
  val aex    = a(manW+exW-1, manW)
  val bex    = b(manW+exW-1, manW)
  val asgn   = a(spec.W-1)
  val bsgn   = b(spec.W-1)

  val azero = aex === 0.U
  val bzero = bex === 0.U

  val zProd = amanW1 * bmanW1
  val zProdW = manW+1 + manW+1
  assert(zProdW == zProd.getWidth)

  val zMoreThan2 = zProd(zProdW-1)
  val zShifted = Mux(zMoreThan2, zProd(zProdW-2, zProdW-manW-2),
                                 zProd(zProdW-3, zProdW-manW-3))
  val zLsb      = zShifted(1)
  val zRound    = zShifted(0)
  val zSticky   = (zProd(zProdW-manW-3) & zMoreThan2) | zProd(zProdW-manW-4, 0).orR
  val zRoundInc = FloatChiselUtil.roundIncBySpec(roundSpec, zLsb, zRound, zSticky)

  val zRounded = zShifted(manW, 1) +& zRoundInc
  val zMoreThan2AfterRound = zRounded(manW)
  val zman0 = Mux(zMoreThan2AfterRound, 0.U(manW.W), zRounded(manW-1, 0))

  val zexInc = zMoreThan2 + zMoreThan2AfterRound
  val zexSum = aex +& bex + zexInc - exBias.U
  val zex0   = zexSum(exW-1, 0)
  val zinf   = zexSum(exW)

  val zzero = azero || bzero
  val zsgn  = asgn ^ bsgn
  val zex   = Mux(zzero, 0.U, Mux(zinf, Fill(exW, 1.U(1.W)), zex0))
  val zman  = Mux(zzero, 0.U, Mux(zinf, 0.U,                 zman0))

  val z = Cat(zsgn, zex, zman)
  assert(z.getWidth == spec.W)

  io.z := z
}
