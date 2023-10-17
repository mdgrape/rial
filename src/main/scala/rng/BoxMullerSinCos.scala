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

// sin/cos(2piy)
class HTBoxMullerSinCos2PiPreProc(
  val cfg: HTBoxMullerConfig,
  val isSin: Boolean
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val manW = realSpec.manW
  val exBias = realSpec.exBias

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

  // round 2pix to [0, pi/4]
  val x0piover2 = io.x.ex < (exBias-2).U
  val x1piover2 = (exBias-2).U === io.x.ex
  val x2piover2 = (exBias-1).U === io.x.ex && (io.x.man.head(1) === 0.U)
  val x3piover2 = (exBias-1).U === io.x.ex && (io.x.man.head(1) === 1.U)
  when(io.en) {assert(x0piover2 || x1piover2 || x2piover2 || x3piover2)}

  // ^   _           _ cos
  // |'.' '.       .'
  // |/ \   \     /    sin
  // +---\---\---/---/->
  //      \   \ /   /
  //       '._.'._.'
  // | 0 | 1 | 2 | 3 |
  // <--->
  // polynomial input [0,1)

  val xAligned = WireDefault(0.U((1+manW).W))
  if(isSin) {
    when(x0piover2) {
      xAligned := Cat(1.U(1.W), io.x.man) >> ((exBias-2).U - io.x.ex)
    }.elsewhen(x1piover2) {
      xAligned := ~io.x.man + 1.U
    }.elsewhen(x2piover2) {
      xAligned := Cat(io.x.man.tail(1), 0.U(1.W))
    }.otherwise {
      xAligned := ~Cat(io.x.man.tail(1), 0.U(1.W)) + 1.U
    }
  } else {
    when(x0piover2) {
      xAligned := ~(Cat(1.U(1.W), io.x.man) >> ((exBias-2).U - io.x.ex)) + 1.U
    }.elsewhen(x1piover2) {
      xAligned := io.x.man
    }.elsewhen(x2piover2) {
      xAligned := ~Cat(io.x.man.tail(1), 0.U(1.W)) + 1.U
    }.otherwise {
      xAligned := Cat(io.x.man.tail(1), 0.U(1.W))
    }
  }

  val adr  = enableIf(io.en, xAligned(manW-1, dxW))
  io.adr := ShiftRegister(adr, nStage)

  if(order != 0) {
    val dx   = enableIf(io.en, Cat(~xAligned(dxW-1), xAligned(dxW-2, 0)))
    io.dx.get := ShiftRegister(dx, nStage)
  }
}

object HTBoxMullerSinCos2PiTableCoeff {
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

class HTBoxMullerSinCos2PiTableCoeff(
  val cfg: HTBoxMullerConfig
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val adrW  = polySpec.adrW
  val fracW = polySpec.fracW
  val dxW   = polySpec.dxW
  val order = polySpec.order

  val maxCbit = HTBoxMullerSinCos2PiTableCoeff.getCBits(polySpec)

  val io = IO(new Bundle {
    val en  = Input(Bool())
    val adr = Input(UInt(adrW.W))
    val cs  = Flipped(new TableCoeffInput(maxCbit))
  })

  if(order == 0) {

    val tableI = VecInit((0L to 1L<<adrW).map(
      n => {
        val x = (n.toDouble / (1L<<adrW)) * 0.25
        val res = if(x == 0) {
          2 * Pi / 8
        } else {
          sin(2 * Pi * x) / (8 * x)
        }
        assert(0.0 <= res && res <= 1.0, f"sin(2pix)/8x = ${res}")

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

    val tableI = HTBoxMullerSinCos2PiTableCoeff.genTable(polySpec)
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


class HTBoxMullerSinCos2PiOtherPath(
  val cfg: HTBoxMullerConfig,
  val isSin: Boolean
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val exW  = realSpec.exW
  val manW = realSpec.manW
  val exBias = realSpec.exBias

  val fracW = polySpec.fracW
  val extraBits = polySpec.extraBits

  val preStage  = cfg.polyPreStage.total
  val calcStage = cfg.polyCalcStage.total
  val postStage = cfg.polyPostStage.total
  val pcGap: Int = if(cfg.preCalcGap)   { 1 } else { 0 }
  val tcGap: Int = if(cfg.tableCalcGap) { 1 } else { 0 }
  val cpGap: Int = if(cfg.calcPostGap)  { 1 } else { 0 }

  val nStage = preStage + pcGap + tcGap + calcStage + cpGap + postStage

  val io = IO(new Bundle {
    val en    = Input(Bool())
    val x     = Flipped(new DecomposedRealOutput(realSpec))
    val coef  = Output(UInt(realSpec.W.W))
    val zsgn  = Output(Bool())
    val zzero = Output(Bool())
    val zone  = Output(Bool())
  })

  val x0piover2 = io.x.ex < (exBias-2).U
  val x1piover2 = (exBias-2).U === io.x.ex
  val x2piover2 = (exBias-1).U === io.x.ex && (io.x.man.head(1) === 0.U)
  val x3piover2 = (exBias-1).U === io.x.ex && (io.x.man.head(1) === 1.U)
  when(io.en) {assert(x0piover2 || x1piover2 || x2piover2 || x3piover2)}

  // ^   _           _ cos
  // |'.' '.       .'
  // |/ \   \     /    sin
  // +---\---\---/---/->
  //      \   \ /   /
  //       '._.'._.'
  // | 0 | 1 | 2 | 3 |
  // <--->
  // polynomial input [0,1)

  val xman0    = io.x.man === 0.U
  val xmanhalf = io.x.man.head(1) === 1.U && io.x.man.tail(1) === 0.U

  if(isSin) {
    io.zsgn  := ShiftRegister(x2piover2 || x3piover2, nStage)
    io.zzero := ShiftRegister((io.x.ex === 0.U) || (io.x.ex === (exBias-1).U && xman0), nStage)
    io.zone  := ShiftRegister((io.x.ex === (exBias-2).U && xman0) ||
                              (io.x.ex === (exBias-1).U && xmanhalf), nStage)

    val yConv = WireDefault(0.U(manW.W))
    val yex   = WireDefault(0.U(exW.W))
    val yman  = WireDefault(0.U(manW.W))

    when(x0piover2) {
      yex  := io.x.ex
      yman := io.x.man
    }.otherwise {
      when(x1piover2) {
        yConv := ~io.x.man + 1.U
      }.elsewhen(x2piover2) {
        yConv := Cat(io.x.man, 0.U(1.W))(manW-1, 0)
      }.otherwise {
        yConv := ~Cat(io.x.man.tail(1), 0.U(1.W)) + 1.U
      }
      val yShift   = PriorityEncoder(Reverse(yConv))
      val yAligned = (yConv << yShift)(manW-1, 0)
      when(io.en) {assert(yAligned.head(1) === 1.U, "yConv = %b, yShift = %d, yAligned = %b\n", yConv, yShift, yAligned)}
      yex  := (exBias-3).U - yShift
      yman := Cat(yAligned, 0.U(1.W))(manW-1, 0)
    }
    io.coef := ShiftRegister(Cat(0.U(1.W), yex, yman), nStage)

  } else {

    io.zsgn  := ShiftRegister(x1piover2 || x2piover2, nStage)
    io.zzero := ShiftRegister((io.x.ex === (exBias-2).U && xman0) ||
                              (io.x.ex === (exBias-1).U && xmanhalf), nStage)
    io.zone  := ShiftRegister((io.x.ex === 0.U) || (io.x.ex === (exBias-1).U && xman0), nStage)

    val yConv = WireDefault(0.U(manW.W))
    val yex   = WireDefault(0.U(exW.W))
    val yman  = WireDefault(0.U(manW.W))
    when(x0piover2) {
      val xAligned = Cat(1.U(1.W), io.x.man) >> ((exBias-2).U - io.x.ex)
      when(io.en) {assert(xAligned.head(1) === 0.U)}
      yConv := ~xAligned(manW-1, 0) + 1.U
    }.elsewhen(x1piover2) {
      yConv := io.x.man
    }.elsewhen(x2piover2) {
      yConv := ~Cat(io.x.man.tail(1), 0.U(1.W)) + 1.U
    }.otherwise {
      yConv :=  Cat(io.x.man.tail(1), 0.U(1.W))
    }
    val yShift   = PriorityEncoder(Reverse(yConv))
    val yAligned = (Cat(yConv, 0.U(1.W)) << yShift)(manW-1, 0)
    yex  := (exBias-3).U - yShift
    yman := yAligned
    io.coef := ShiftRegister(Cat(0.U(1.W), yex, yman), nStage)
  }
}

class HTBoxMullerSinCos2PiPostProcess(
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
    val en   = Input(Bool())
    val zres = Input(UInt(fracW.W))       // [0.5,1), sin(2pix) / 8x
    val z    = Output(UInt(realSpec.W.W)) // sin(2pix)/x [4, 8)
  })
  when(io.en) {assert(io.zres.head(1) === 1.U)}

  val zsgn = 0.U(1.W)
  val zex  = (exBias+2).U(exW.W)

  val zman = WireDefault(0.U(manW.W))
  if(extraBits == 0) {
    zman := Cat(io.zres.tail(1), 0.U(1.W))
  } else if(extraBits == 1) {
    zman := io.zres.tail(1)
  } else {
    val z0    = io.zres.tail(1)
    val zman0 = dropLSB(extraBits-1, z0) +& z0(extraBits-2)
    val polynomialOvf = zman0(manW)
    zman := Mux(polynomialOvf, Fill(manW, 1.U(1.W)), zman0(manW-1,0))
  }

  val z = enableIf(io.en, Cat(zsgn, zex, zman))
  assert(z.getWidth == realSpec.W)

  io.z := ShiftRegister(z, nStage)
}

class HTBoxMullerSinCos2Pi(
  val cfg: HTBoxMullerConfig,
  val isSin: Boolean
) extends Module {

  val realSpec = cfg.realSpec
  val polySpec = cfg.polySpec

  val manW   = realSpec.manW
  val exW    = realSpec.exW
  val exBias = realSpec.exBias

  val order = polySpec.order

  val maxCbit = HTBoxMullerSinCos2PiTableCoeff.getCBits(polySpec)

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

  val preproc = Module(new HTBoxMullerSinCos2PiPreProc(cfg, isSin))

  preproc.io.en := io.en
  preproc.io.x  := io.x

  // ---------------------------------------------------------------------------

  val table = Module(new HTBoxMullerSinCos2PiTableCoeff(cfg))

  table.io.en  := ShiftRegister(io.en, preproc.nStage + pcGap, false.B, true.B)
  table.io.adr := ShiftRegister(preproc.io.adr, pcGap)
  val polynomialCoeff = ShiftRegister(table.io.cs.asUInt, tcGap)

  // ---------------------------------------------------------------------------

  val eval = Module(new PolynomialEval(realSpec, polySpec, maxCbit, calcStage))

  eval.io.coeffs := polynomialCoeff.asTypeOf(new TableCoeffInput(maxCbit))
  if(order != 0) {
    val dx = preproc.io.dx.get
    eval.io.dx.get := ShiftRegister(dx, pcGap + tcGap)
  }

  // ---------------------------------------------------------------------------

  val postProcInputDelay = preStage.total + pcGap + tcGap + calcStage.total + cpGap
  val postproc = Module(new HTBoxMullerSinCos2PiPostProcess(cfg))

  postproc.io.en   := ShiftRegister(io.en, postProcInputDelay)
  postproc.io.zres := ShiftRegister(eval.io.result, cpGap)

  val sinXoverX = postproc.io.z

  val sinXoverXex  = sinXoverX(manW+exW-1, manW)
  val sinXoverXman = sinXoverX(manW-1, 0)
//   printf("sinXoverX   = (%b|%d|%b (%d/%d))\n", sinXoverX.head(1), sinXoverXex, sinXoverXman, sinXoverXman, (1<<manW).U)

  // ---------------------------------------------------------------------------

  val other = Module(new HTBoxMullerSinCos2PiOtherPath(cfg, isSin))

  other.io.en := io.en
  other.io.x  := io.x

  val coef = other.io.coef
  val coefex  = coef(manW+exW-1, manW)
  val coefman = coef(manW-1, 0)
//   printf("coef        = (%b|%d|%b (%d/%d))\n", coef.head(1), coefex, coefman, coefman, (1<<manW).U)

  // ---------------------------------------------------------------------------

  val multiplier = Module(new MultFPGeneric(realSpec, realSpec, realSpec, RoundSpec.roundToEven, mulStage))

  multiplier.io.x := sinXoverX
  multiplier.io.y := coef

  val z0 = multiplier.io.z

  val zsgn  = ShiftRegister(other.io.zsgn,  mulStage.total, false.B, true.B)
  val zzero = ShiftRegister(other.io.zzero, mulStage.total, false.B, true.B)
  val zone  = ShiftRegister(other.io.zone,  mulStage.total, false.B, true.B)

  val zex  = Mux(zzero, 0.U(exW.W),
             Mux(zone, exBias.U(exW.W), z0(manW+exW-1, manW)))
  val zman = Mux(zone || zzero, 0.U(manW.W), z0(manW-1, 0))

  io.z := Cat(zsgn, zex, zman)
}
