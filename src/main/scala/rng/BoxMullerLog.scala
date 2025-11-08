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
// calc -2ln(x) for x in (0.0, 1.0]
//
// ln(x) = ln2 * log2(x)
//       = ln2 * log2(2^ex * 1.man)
//       = ln2 * (ex + log2(1.man))
//
// case 0.0 <  x < 0.5:
//   calc -(ex + log2(1.man)) * 2ln2
//      = (-ex - log2(1.man)) * 2ln2
//      = ((exBias-xex) - log2(1.man)) * 2ln2
//      = (-ex - 1 + 1 - log2(1.man)) * 2ln2
//                   ^^^^^^^^^^^^^^^
//                   polynomial (0.0, 1.0]
//
// case 0.5 <= x < 1.0:
//   calc -ln(x)/(1-x) * 2(1-x)
//        ^^^^^^^^^^^^   ^^^^^^
//        polynomial     (0, 1]
//
// case x = 1.0:
//   return 0.0
//
private[rial] class HTBoxMullerLogPreProc(
  val cfg: HTBoxMullerConfig
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val manW = realSpec.manW

  val adrW  = polySpec.adrW
  val dxW   = polySpec.dxW
  val order = polySpec.order

  val nStage  = cfg.polyPreStage.total

  val io = IO(new Bundle {
    val en  = Input  (Bool())
    val x   = Flipped(new DecomposedRealOutput(realSpec))
    val adr = Output(UInt(adrW.W))
    val dx  = if(order != 0) { Some(Output(UInt(dxW.W))) } else { None }
  })

  val adr  = enableIf(io.en, io.x.man(manW-1, dxW))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx   = enableIf(io.en, Cat(~io.x.man(dxW-1), io.x.man(dxW-2, 0)))
    io.dx.get := ShiftRegister(dx, nStage)
  }
}

private[rial] object HTBoxMullerLog2TableCoeff{

  def genTable(
    polySpec: PolynomialSpec,
    calcW: Option[Seq[Int]] = None,
    cbit:  Option[Seq[Int]] = None
  ): FuncTableInt = {

    val f = (x0: Double) => {
      val x = 1.0 + x0 // [1, 2.0) -log2-> [0, 1)
      val log2x = log(x) / log(2.0)
      min(1.0 - log2x, 1.0 - pow(2.0, -53)) // avoid returning 1.0
    }

    val order   = polySpec.order
    val fracW   = polySpec.fracW
    val adrW    = polySpec.adrW

    val tableD = new FuncTableDouble( f, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)

    new FuncTableInt( tableD, fracW, calcW, cbit )
  }

  def getCBits(
    polySpec: PolynomialSpec
  ): Seq[Int] = {
    val tableI = genTable(polySpec)
    tableI.cbit
  }
  def getCalcW(
    polySpec: PolynomialSpec
  ): Seq[Int] = {
    val tableI = genTable(polySpec)
    tableI.calcWidth
  }
}

private[rial] object HTBoxMullerLnTableCoeff{

  def genTable(
    polySpec: PolynomialSpec,
    calcW: Option[Seq[Int]] = None,
    cbit:  Option[Seq[Int]] = None
  ): FuncTableInt = {

    val f = (x0: Double) => {
      val x = 0.5 * (1.0 + x0)
      if(x0 == 1.0) {
        0.0
      } else {
        -log(x) / (1.0 - x) - 1.0
      }
    }

    val order   = polySpec.order
    val fracW   = polySpec.fracW
    val adrW    = polySpec.adrW

    val tableD = new FuncTableDouble( f, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)

    new FuncTableInt( tableD, fracW, calcW, cbit )
  }

  def getCBits(
    polySpec: PolynomialSpec
  ): Seq[Int] = {
    val tableI = genTable(polySpec)
    tableI.cbit
  }
  def getCalcW(
    polySpec: PolynomialSpec
  ): Seq[Int] = {
    val tableI = genTable(polySpec)
    tableI.calcWidth
  }
}

private[rial] class HTBoxMullerLog2TableCoeff(
  val cfg: HTBoxMullerConfig
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val adrW  = polySpec.adrW
  val fracW = polySpec.fracW
  val dxW   = polySpec.dxW
  val order = polySpec.order

  val maxCbit = Seq(
    HTBoxMullerLog2TableCoeff.getCBits(polySpec),
    HTBoxMullerLnTableCoeff.getCBits(polySpec)
  ).reduce( (lhs, rhs) => {
    lhs.zip(rhs).map( c => max(c._1, c._2) )
  } )
  val maxCalcW = Seq(
    HTBoxMullerLog2TableCoeff.getCalcW(polySpec),
    HTBoxMullerLnTableCoeff.getCalcW(polySpec)
  ).reduce( (lhs, rhs) => {
    lhs.zip(rhs).map( c => max(c._1, c._2) )
  } )

  val io = IO(new Bundle {
    val en  = Input(Bool())
    val adr = Input(UInt(adrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {

    val tableI = VecInit((0L to 1L<<adrW).map(
      n => {
        val x = (n.toDouble / (1L<<adrW)) + 1.0
        val log2x = log(x) / log(2.0)
        val res = 1.0 - log2x

        assert(0.0 <= res && res <= 1.0, f"log2(x) = ${res}")
        val y = round(res * (1L<<fracW))
        if (y >= (1L<<fracW)) {
          maskL(fracW).U(fracW.W)
        } else if (y <= 0.0) {
          0.U(fracW.W)
        } else {
          y.U(fracW.W)
        }
      })
    )

    val coeff = tableI(io.adr)
    io.cs.cs(0) := enableIf(io.en, coeff)

  } else {

    val tableI = HTBoxMullerLog2TableCoeff.genTable(polySpec, Some(maxCalcW), Some(maxCbit))
    val cbit   = tableI.cbit
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(0 <= diffWidth, f"diffWidth = ${diffWidth} must be > 0")

      if(diffWidth != 0) {
        val ci  = coeff(i)
        val msb = ci(cbit(i)-1)
        coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffs.cs(i) := coeff(i)
      }
    }
    io.cs := enableIf(io.en, coeffs)
  }
}

private[rial] class HTBoxMullerLnTableCoeff(
  val cfg: HTBoxMullerConfig
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val adrW  = polySpec.adrW
  val fracW = polySpec.fracW
  val dxW   = polySpec.dxW
  val order = polySpec.order

  val maxCbit = Seq(
    HTBoxMullerLog2TableCoeff.getCBits(polySpec),
    HTBoxMullerLnTableCoeff.getCBits(polySpec)
  ).reduce( (lhs, rhs) => {
    lhs.zip(rhs).map( c => max(c._1, c._2) )
  } )
  val maxCalcW = Seq(
    HTBoxMullerLog2TableCoeff.getCalcW(polySpec),
    HTBoxMullerLnTableCoeff.getCalcW(polySpec)
  ).reduce( (lhs, rhs) => {
    lhs.zip(rhs).map( c => max(c._1, c._2) )
  } )

  val io = IO(new Bundle {
    val en  = Input(Bool())
    val adr = Input(UInt(adrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {

    val tableI = VecInit((0L to 1L<<adrW).map(
      n => {
        val x = 0.5 * ((n.toDouble / (1L<<adrW)) + 1.0) // 0~1 -> 0.5~1
        val res = if(n == 1L<<adrW) {0.0} else {-log(x) / (1.0 - x) - 1.0}
        assert(0.0 <= res && res <= 1.0, f"log(x)/(x-1) = ${res}")
        val y = round(res * (1L<<fracW))

        if (y >= (1L<<fracW)) {
          maskL(fracW).U(fracW.W)
        } else if (y <= 0.0) {
          0.U(fracW.W)
        } else {
          y.U(fracW.W)
        }
      })
    )

    val coeff = tableI(io.adr)
    io.cs.cs(0) := enableIf(io.en, coeff)

  } else {

    val tableI = HTBoxMullerLnTableCoeff.genTable(polySpec, Some(maxCalcW), Some(maxCbit))
    val cbit   = tableI.cbit
    val (coeffTable, coeffWidth) = tableI.getVectorUnified(/*sign mode =*/0)
    val coeff  = getSlices(coeffTable(io.adr), coeffWidth)

    val coeffs = Wire(new TableCoeffInput(maxCbit))
    for (i <- 0 to order) {
      val diffWidth = maxCbit(i) - cbit(i)
      assert(0 <= diffWidth)

      if(diffWidth != 0) {
        val ci  = coeff(i)
        val msb = ci(cbit(i)-1)
        coeffs.cs(i) := Cat(Fill(diffWidth, msb), ci) // sign extension
      } else {
        coeffs.cs(i) := coeff(i)
      }
    }
    io.cs := enableIf(io.en, coeffs)
  }
}

// case 0.0 <  x < 0.5:
//   calc -(ex + log2(1.man)) * 2ln2
//      = (-ex - 1 + 1 - log2(1.man)) * 2ln2
//   polynomial module calculates 1.0-log2(1.man) and result is in (0.0, 1.0].
//   if man === 0.U, 1-log2(1.0) == 1 and it overflows. we must check it.
//   otherwise, polynomial result is in (0.0, 1.0) and -ex-1 becomes the integer
//   part. then we normalize -ex-1+zres and multiply it with 2ln2.
//
private[rial] class HTBoxMullerLog2PostProcess(
  val cfg: HTBoxMullerConfig
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val exW  = realSpec.exW
  val manW = realSpec.manW
  val exBias = realSpec.exBias

  val fracW = polySpec.fracW

  val nStage = cfg.polyPostStage.total

  val io = IO(new Bundle {
    val en       = Input(Bool())
    val zres     = Input(UInt(fracW.W)) // 1 - log2(1.xman)
    val xManZero = Input(Bool())      // x.man === 0.U
    val xExM1    = Input(UInt(exW.W)) // -(x.ex-exBias)-1
    val z        = Output(UInt(realSpec.W.W))
  })

  // always positive.
  val zsgn = 0.U(1.W)

  // -2ln(x) = -log2(2^ex * 1.man) * 2ln2
  //         = -(ex + log2(1.man)) * 2ln2
  //         = (-ex - log2(1.man)) * 2ln2
  //         = (-ex - 1 + 1 - log2(1.man)) * 2ln2
  //            ^^^^^^^   ^^^^^^^^^^^^^^^
  //            zInt      zFrac (0, 1]. == 1 only if man == 0
  //
  // since x < 1/2, ex < -1. -ex > 1, so -ex - 1 > 0.

  val zInt  = Mux(io.xManZero, io.xExM1 + 1.U, io.xExM1)
  val zFrac = Mux(io.xManZero, 0.U(fracW.W),   io.zres)

  when(io.en) {assert(io.xExM1 =/= 0.U)}

  val z0Shift = PriorityEncoder(Reverse(zInt))

  val zex = (exBias + exW - 1).U(exW.W) - z0Shift

  val z0  = Cat(zInt, zFrac)
  val z0W = z0.getWidth
  val z0Aligned = (z0 << z0Shift)(z0W-1, 0)

  assert(z0W == exW + fracW)

  when(io.en) { assert(z0Aligned.head(1) === 1.U) } // XXX

  val zman0 = z0Aligned.tail(1) // remove hidden bit

  val extraBits = zman0.getWidth - manW
  val z0Rounded = dropLSB(extraBits, zman0) +& zman0(extraBits-1)
  val polynomialOvf = z0Rounded(manW)

  val zman = Mux(polynomialOvf, Fill(manW, 1.U(1.W)), z0Rounded(manW-1, 0))

  val z = enableIf(io.en, Cat(zsgn, zex, zman))
  assert(z.getWidth == realSpec.W)

//  printf("HTBoxMullerLog2PostProcess: zInt = %b, zFrac = %b (%d / %d)\n", zInt, io.zres, io.zres, (1<<fracW).U)
//  printf("HTBoxMullerLog2PostProcess: z0        = %b, z0Shift = %d\n", z0, z0Shift)
//  printf("HTBoxMullerLog2PostProcess: z0Aligned = %b\n", z0Aligned)

  io.z := ShiftRegister(z, nStage)
}

// case 0.5 <= x < 1.0:
//   calc -ln(x)/(1-x) * 2(1-x)
//        ^^^^^^^^^^^^
//        (1.0, 2.0)
//   polynomial calculates -ln(x)/(1-x) - 1 to fit it in [0.0, 1.0).
//   so we need to add hidden bit to zres, round it, and multiply with 2(1-x).
//
private[rial] class HTBoxMullerLnPostProcess(
  val cfg: HTBoxMullerConfig
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val exW  = realSpec.exW
  val manW = realSpec.manW
  val exBias = realSpec.exBias

  val fracW = polySpec.fracW
  val extraBits = polySpec.extraBits

  val nStage = cfg.polyPostStage.total

  val io = IO(new Bundle {
    val en       = Input(Bool())
    val zres     = Input(UInt(fracW.W))       // -ln(x)/(1-x) - 1
    val z        = Output(UInt(realSpec.W.W)) // -ln(x)/(1-x)
  })

  // always in [1.0, 2.0).
  val zsgn = 0.U(1.W)
  val zex  = exBias.U(exW.W)

  val zman = WireDefault(0.U(manW.W))
  if(extraBits == 0) {
    zman := io.zres
  } else {
    val zman0 = dropLSB(extraBits, io.zres) +& io.zres(extraBits-1)
    val polynomialOvf = zman0(manW)
    zman := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW-1,0))
  }

  val z = enableIf(io.en, Cat(zsgn, zex, zman))
  assert(z.getWidth == realSpec.W)

  io.z := ShiftRegister(z, nStage)
}

private[rial] class HTBoxMullerLogMultiplier(
  val cfg: HTBoxMullerConfig,
) extends Module {

  val realSpec = cfg.realSpec

  val exW  = realSpec.exW
  val manW = realSpec.manW

  val nStage = cfg.mulStage.total

  val io = IO(new Bundle {
    val en = Input(Bool())
    val x  = Input(UInt(realSpec.W.W))
    val y  = Input(UInt((1+manW).W))
    val z  = Output(UInt(realSpec.W.W)) // -ln(x)/(1-x)
  })

  val xex  = io.x(manW+exW-1, manW)
  val xman = io.x(manW-1, 0)

  val lhs = Cat(1.U(1.W), xman)
  val rhs = io.y
  val mul = lhs * rhs
  val mulW = mul.getWidth
  assert(mulW == manW+1 + manW+1)

  // 0001 -> 0
  // 001x -> 1
  // 01xx -> 2
  val mulShift = PriorityEncoder(Reverse(mul))

  val mulAligned = (mul << mulShift)(mulW-1, 0)
  when(io.en) {assert(mulAligned.head(1) === 1.U)}

  val lsb    = mulAligned(mulW-manW-1)
  val round  = mulAligned(mulW-manW-2)
  val sticky = mulAligned(mulW-manW-3, 0).orR
  val inc = FloatChiselUtil.roundIncBySpec(RoundSpec.roundToEven, lsb, round, sticky)

  val zmanRound = mulAligned(mulW-1, mulW-manW-1) +& inc
  val zmanMoreThan2 = zmanRound.head(1)

  val zman = Mux(zmanMoreThan2 === 1.U, 0.U(manW.W), zmanRound(manW-1, 0))

  val zsgn = 0.U(1.W)
  val zex = xex + 1.U + zmanMoreThan2 - mulShift

  val z0 = enableIf(io.en, Cat(zsgn, zex, zman))
  assert(z0.getWidth == realSpec.W)

  val z = Mux(io.x === 0.U, 0.U(realSpec.W.W), z0)

//  printf("HTBoxMullerLogMultiplier: xex=%d,xman=%b, y=%b\n", xex, xman, io.y)
//  printf("HTBoxMullerLogMultiplier: xexNoBias=%d, xfrac = %d / %d, y= %d / %d\n", xex - realSpec.exBias.U, lhs, 1.U << manW, rhs, 1.U << manW)
//  printf("HTBoxMullerLogMultiplier: mul        = %b, mulW     = %d, manW     = %d\n", mul, mulW.U, manW.U)
//  printf("HTBoxMullerLogMultiplier: mulAligned = %b, mulShift = %d\n", mulAligned, mulShift)
//  printf("HTBoxMullerLogMultiplier: zmanRound = %b, zman = %b, zex = %d\n", zmanRound, zman, zex)

  io.z := ShiftRegister(z, nStage)
}

// ============================================================================
//
//                     2ln2 -.
// x -+-> -(ex+log2(1.man)) ----.
//    +-> -ln(x)/(1-x)      ---mux--> mul -> -2ln(x)
//    '-> 2(1-x)            -'
//
private[rial] class HTBoxMullerLog(
  val cfg: HTBoxMullerConfig
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val manW   = realSpec.manW
  val exW    = realSpec.exW
  val exBias = realSpec.exBias

  val order = polySpec.order

  val maxCbit = Seq(
    HTBoxMullerLog2TableCoeff.getCBits(polySpec),
    HTBoxMullerLnTableCoeff.getCBits(polySpec)
  ).reduce( (lhs, rhs) => {
    lhs.zip(rhs).map( c => max(c._1, c._2) )
  } )

  val preStage  = cfg.polyPreStage
  val calcStage = cfg.polyCalcStage
  val postStage = cfg.polyPostStage
  val pcGap: Int = if(cfg.preCalcGap)   { 1 } else { 0 }
  val tcGap: Int = if(cfg.tableCalcGap) { 1 } else { 0 }
  val cpGap: Int = if(cfg.calcPostGap)  { 1 } else { 0 }
  val mulStage = cfg.mulStage

  def nStage(): Int = {
    preStage.total + pcGap + tcGap + calcStage.total + cpGap + postStage.total +
    mulStage.total
  }

  // ---------------------------------------------------------------------------

  val io = IO(new Bundle {
    val en = Input  (Bool())
    val x  = Flipped(new DecomposedRealOutput(realSpec))
    val z  = Output (UInt(realSpec.W.W))
  })

  val xIsZero     = io.x.zero
  val xIsOne      = io.x.ex === exBias.U     // return 0.0
  val xIsOverHalf = io.x.ex === (exBias-1).U // use ln(x)/(1-x) part.
                                             // otherwise, use log2 part.

  val useln   = io.en && !xIsOne &&  xIsOverHalf
  val uselog2 = io.en && !xIsOne && !xIsOverHalf

//  printf("HTBoxMullerLog: useln = %b, uselog2 = %b\n", useln, uselog2)
//  printf("HTBoxMullerLog: xex   = %d, xman = %b (%d / %d)\n", io.x.ex, io.x.man, Cat(1.U(1.W), io.x.man), (1<<manW).U)

  // ---------------------------------------------------------------------------

  val preproc = Module(new HTBoxMullerLogPreProc(cfg))

  preproc.io.en := io.en
  preproc.io.x  := io.x

  // ---------------------------------------------------------------------------

  val log2table = Module(new HTBoxMullerLog2TableCoeff(cfg))
  val lntable   = Module(new HTBoxMullerLnTableCoeff(cfg))

  log2table.io.en  := ShiftRegister(uselog2, preproc.nStage + pcGap, false.B, true.B)
  log2table.io.adr := ShiftRegister(preproc.io.adr, pcGap)

  lntable.io.en  := ShiftRegister(useln, preproc.nStage + pcGap, false.B, true.B)
  lntable.io.adr := ShiftRegister(preproc.io.adr, pcGap)

  val log2coeff = log2table.io.cs.asUInt
  val lncoeff   = lntable.io.cs.asUInt

  val polynomialCoeff = ShiftRegister(log2coeff | lncoeff, tcGap)

//  printf("HTBoxMullerLog: log2coeff = %b\n", log2coeff)
//  printf("HTBoxMullerLog: ln  coeff = %b\n", lncoeff)
//  printf("HTBoxMullerLog: coeff     = %b\n", polynomialCoeff)

  // ---------------------------------------------------------------------------

  val eval = Module(new PolynomialEval(realSpec, polySpec, maxCbit, calcStage))

  eval.io.coeffs := polynomialCoeff.asTypeOf(new TableCoeffInput(maxCbit))
  if(order != 0) {
    val dx = preproc.io.dx.get
    eval.io.dx.get := ShiftRegister(dx, pcGap + tcGap)
  }

  // ---------------------------------------------------------------------------

  val postProcInputDelay = preStage.total + pcGap + tcGap + calcStage.total + cpGap
  val log2postproc = Module(new HTBoxMullerLog2PostProcess(cfg))
  val lnpostproc   = Module(new HTBoxMullerLnPostProcess(cfg))

  log2postproc.io.en   := ShiftRegister(uselog2, postProcInputDelay)
  log2postproc.io.zres := ShiftRegister(eval.io.result, cpGap)
  log2postproc.io.xManZero := ShiftRegister(io.x.man === 0.U, postProcInputDelay)
  log2postproc.io.xExM1    := ShiftRegister((exBias-1).U - io.x.ex, postProcInputDelay)
  // -(x.ex-exBias) - 1 == -x.xex + exBias - 1

  lnpostproc.io.en   := ShiftRegister(useln, postProcInputDelay)
  lnpostproc.io.zres := ShiftRegister(eval.io.result, cpGap)

  val zIsZero = ShiftRegister(xIsOne, postProcInputDelay + postStage.total)

  val z0 = Mux(zIsZero, 0.U(realSpec.W.W), log2postproc.io.z | lnpostproc.io.z)

  // --------------------------------------------------------------------------

  val ln2x2    = new RealGeneric(realSpec, 2.0 * log(2.0))
  val log2Coef = Cat(1.U(1.W), (ln2x2.man.toBigInt).U(manW.W))
  val lnCoef0  = Cat(1.U(1.W), 0.U((manW+1).W)) - Cat(1.U(2.W), io.x.man) // 2(1-x)
  val lnCoef   = ShiftRegister(lnCoef0(manW+1-1, 0), postProcInputDelay + postStage.total)
  when(io.en) {assert(lnCoef0.head(1) === 0.U)}
  val uselog2coef = ShiftRegister(uselog2, postProcInputDelay + postStage.total)

  val coef = Mux(uselog2coef, log2Coef, lnCoef)

  val multiplier = Module(new HTBoxMullerLogMultiplier(cfg))

  multiplier.io.en := ShiftRegister(io.en, postProcInputDelay + postStage.total)
  multiplier.io.x  := z0
  multiplier.io.y  := coef

  val xWasZero = ShiftRegister(xIsZero, nStage)

  val fmax = Cat(0.U(1.W), realSpec.exMax.U(realSpec.exW.W), Fill(realSpec.manW, 1.U(1.W)))
  io.z := Mux(xWasZero, fmax, multiplier.io.z)
}
