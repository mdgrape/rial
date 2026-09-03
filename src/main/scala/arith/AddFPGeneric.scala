//% @file AddFPGeneric.scala
//
// Generic Floating-Point Adder Unit
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.arith

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._
import chisel3.util.experimental.BoringUtils
import rial.table._
import rial.util._
import rial.util.RialChiselUtil._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig._
import rial.util.DebugControl

/** Generic Floating-Point Adder.
 *
 * It takes `x` and `y` and returns `z` as `z = x + y`.
 *
 * All the input/output floating-point spec are configurable through [[rial.arith.RealSpec]].
 *
 * The latency of this module can be configured via [[rial.util.PipelineStageConfig]].
 * Unlike the modules that simply put all their registers at their output, this
 * adder declares the boundaries of its datapath, so the registers can be placed
 * inside it - either automatically or by name:
 *
 *  - `"in"`     the module input
 *  - `"unpack"` after decompose / classify / exponent difference / near-far
 *               decision / mantissa padding
 *  - `"swap"`   after the operand swap muxes
 *  - `"align"`  after the far-path alignment shift and sticky, and the
 *               near-path pre-shift
 *  - `"add"`    after the far-path significand add and the near-path subtract
 *  - `"norm"`   after the leading-zero count and normalisation shift, the
 *               mantissa / round selection and the far-near merge
 *  - `"round"`  after the rounding increment and the exponent adjustment
 *  - `"out"`    the module output
 *
 * {{{
 * // three registers, placed by the module
 * new AddFPGeneric(spec, spec, spec, RoundSpec.roundToEven, 3)
 * // hand-picked positions
 * new AddFPGeneric(spec, spec, spec, RoundSpec.roundToEven,
 *   PipelineStageConfig.specify("add" -> 1, "norm" -> 1, "round" -> 1))
 * // everything at the output, as before
 * new AddFPGeneric(spec, spec, spec, RoundSpec.roundToEven, PipelineStageConfig.atOut(3))
 * }}}
 *
 * {{{
 * class AddFPGeneric(...) extends Module {
 *   val io = IO(new Bundle{
 *     val x = Input (UInt(xSpec.W.W))
 *     val y = Input (UInt(ySpec.W.W))
 *     val z = Output(UInt(zSpec.W.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * @constructor create a chisel Module AddFPGeneric
 * @param xSpec          floating point spec of the left hand side input
 * @param ySpec          floating point spec of the right hand side input
 * @param zSpec          floating point spec of the output
 * @param roundSpec      rounding algorithm
 * @param stage          number of pipeline stages in this module.
 * @param enableDebug
 *   add some additional printf and assertions. by default, false.
 * @param useParallelLZA
 *   do LZA for faster latency. by default, false.
 *
 * For parameters, see also:
 * @see [[rial.arith.RealSpec]]
 * @see [[rial.arith.RoundSpec]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 *
 */
class AddFPGeneric(
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig,
  val enableDebug : Boolean = false,
  val useParallelLZA : Boolean = false
) extends Module with DebugControl {

  /** Boundaries of the datapath at which pipeline registers can be placed.
   * The weights are rough, target-independent estimates of the delay of the
   * logic that ends at each boundary. */
  val stagePoints = AddFPGeneric.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getParam = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x   = Input(UInt(xSpec.W.W))
    val y   = Input(UInt(ySpec.W.W))
    val z   = Output(UInt(zSpec.W.W))
  })

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  val xIn = ShiftRegister(io.x, stages.at("in"))
  val yIn = ShiftRegister(io.y, stages.at("in"))

  //----------------------------------------------------------------------
  // S1 : decompose / classify / exponent difference / mantissa padding
  val (xsgn, xex, xman) = FloatChiselUtil.decompose(xSpec, xIn)
  val (ysgn, yex, yman) = FloatChiselUtil.decompose(ySpec, yIn)
  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, xIn)
  val (yzero, yinf, ynan) = FloatChiselUtil.checkValue(ySpec, yIn)
  val diffSgn = (xsgn ^ ysgn).asBool

  // inf + inf = inf, -inf - inf = -inf, inf - inf = nan
  val xyNaN  = xnan || ynan || (xinf && yinf && xsgn =/= ysgn)
  val xyInf  = (xinf || yinf) && !xyNaN
  val xyInfSgn = ! ((xinf && !xsgn) || (yinf && !ysgn))
  val xyBothZero = xzero && yzero
  dbgPrintf("%x %x\n",  xIn, yIn)
  dbgPrintf("%b %b %b\n",  xyInf, xyNaN, xyBothZero)

  //----------------------------------------------------------------------
  // Exponent comparison
  val diffMaxXY = xSpec.exMax - ySpec.exMin // x:max, y:min
  val diffMaxYX = ySpec.exMax - xSpec.exMin // y:max, x:min
  val diffExW = log2Up(max(diffMaxXY,diffMaxYX)+1)+1
  val diffExXY = ( xex.pad(diffExW).asSInt-yex.pad(diffExW).asSInt
    +(-xSpec.exBias+ySpec.exBias).S(diffExW.W) )
  // The following is just negative of diffExXY but in some case the separate calculation
  // decreases the latency
  val diffExYX = ( yex.pad(diffExW).asSInt-xex.pad(diffExW).asSInt
    +(xSpec.exBias-ySpec.exBias).S(diffExW.W) )
  val near = (!xzero) && (!yzero) && diffSgn && ((diffExXY(diffExW-1,1) === 0.U)|| (diffExXY=== -1.S)) // abs(xEx-yEx) <= 1
  val nearShift = diffExXY(0)

  // exponent of x >= y
  // - x=FP32 -127(==zero) + y=FP64 -300 => zex should be -300
  //   diffExXY > 0: true
  //   xzero: true
  //   yzero: false
  // - x=FP64 -300 + y=FP32 -127(==zero) => zex should be -300
  //   diffExXY > 0: false
  //   xzero: false
  //   yzero: true
  val exXlargerThanY = (diffExXY(diffExW-1) === 0.U && (!xzero)) || yzero

  val zexBaseMax = max( xSpec.exMax + zSpec.exBias, ySpec.exMax + zSpec.exBias)
  val zexBaseMin = max( xSpec.exMin + zSpec.exBias, ySpec.exMin + zSpec.exBias)
  val zexBaseW = getWidthToRepresentNumbersAlwaysSigned( Seq(zexBaseMin,  zexBaseMax) )
  //println(f"diffExW=$diffExW%d zexBaseW = $zexBaseW%d")
  val zexNorm = Wire(SInt(zexBaseW.W))
  if ((zSpec.exBias == xSpec.exBias) && (zSpec.exBias == ySpec.exBias)) {
    zexNorm := (Mux (exXlargerThanY, xex, yex)).zext
  } else {
    val xexPad = xex.pad(zexBaseW).asSInt
    val yexPad = yex.pad(zexBaseW).asSInt

    val zexDiff = Mux (exXlargerThanY,
      (zSpec.exBias-xSpec.exBias).S(zexBaseW.W),
      (zSpec.exBias-ySpec.exBias).S(zexBaseW.W) )

    zexNorm := Mux (exXlargerThanY, xexPad, yexPad) + zexDiff
  }
  dbgPrintf("near=%b zexNorm=%x diffExXY/XY=%x / %x\n", near, zexNorm, diffExXY, diffExYX)

  //------------------------------------------------------------------------
  // Mantissa padding
  val minManW    = xSpec.manW.min(ySpec.manW)
  val maxManWxy  = xSpec.manW.max(ySpec.manW)
  val maxManWxyz = maxManWxy.max(zSpec.manW)
  //val manXlargerThanY = xman(xSpec.manW-1, xSpec.manW-minManW) > yman(ySpec.manW-1, ySpec.manW-minManW)

  val xmanPadXyz = Mux(xzero, 0.U((maxManWxyz+1).W), Cat(1.U(1.W), padLSB(maxManWxyz, xman)))
  val ymanPadXyz = Mux(yzero, 0.U((maxManWxyz+1).W), Cat(1.U(1.W), padLSB(maxManWxyz, yman)))

  val yzeroFar = Mux(exXlargerThanY, xzero, yzero)

  //---------------------------- stage "unpack" --------------------------
  private def cutUnpack[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("unpack"))
  }
  // Signals read by `dbgPrintf` only must not create registers in a normal
  // (non-debug) build, so they are piped only when the printfs exist.
  private def cutUnpackDbg[T <: Data](x: T): T = {
    if (enableDebug) {
      cutUnpack(x)
    } else {
      x
    }
  }

  val xmanPadXyzUnpack     = cutUnpack(xmanPadXyz)
  val ymanPadXyzUnpack     = cutUnpack(ymanPadXyz)
  val exXlargerThanYUnpack = cutUnpack(exXlargerThanY)
  val diffExXYUnpack       = cutUnpack(diffExXY)
  val diffExYXUnpack       = cutUnpack(diffExYX)
  val nearUnpack           = cutUnpack(near)
  val nearShiftUnpack      = cutUnpack(nearShift)
  val zexNormUnpack        = cutUnpack(zexNorm)
  val diffSgnUnpack        = cutUnpack(diffSgn)
  val xmanUnpack           = cutUnpackDbg(xman) // dbgPrintf only
  val ymanUnpack           = cutUnpackDbg(yman) // dbgPrintf only
  val xsgnUnpack           = cutUnpack(xsgn)
  val ysgnUnpack           = cutUnpack(ysgn)
  val xyNaNUnpack          = cutUnpack(xyNaN)
  val xyInfUnpack          = cutUnpack(xyInf)
  val xyInfSgnUnpack       = cutUnpack(xyInfSgn)
  val xyBothZeroUnpack     = cutUnpack(xyBothZero)

  //----------------------------------------------------------------------
  // S2 : operand swap (larger exponent -> xFar)
  //   Far path
  //   when exponent difference > 1 or sign is the same
  //   in this case, assume x is normalized as 1<=x<2 and abs(y)<x,
  //     -1/2 < y < 2, and then 1/2 < x+y < 4.
  //
  //   currently, separate barrel shifters/adders are used for far/near paths.
  val xFar = Cat(0.U(1.W), Mux(exXlargerThanYUnpack, xmanPadXyzUnpack, ymanPadXyzUnpack))
  val yFar = Mux(exXlargerThanYUnpack, ymanPadXyzUnpack, xmanPadXyzUnpack)
  val xFarSign = Mux(exXlargerThanYUnpack, xsgnUnpack, ysgnUnpack).asBool
  // Drop sign bit
  val diffFar = Mux(exXlargerThanYUnpack, diffExXYUnpack(diffExW-2,0), diffExYXUnpack(diffExW-2,0))

  //---------------------------- stage "swap" ----------------------------
  private def cutSwap[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("swap"))
  }
  private def cutSwapDbg[T <: Data](x: T): T = {
    if (enableDebug) {
      cutSwap(x)
    } else {
      x
    }
  }

  val xFarSwap       = cutSwap(xFar)
  val yFarSwap       = cutSwap(yFar)
  val xFarSignSwap   = cutSwap(xFarSign)
  val diffFarSwap    = cutSwap(diffFar)
  val nearSwap       = cutSwap(nearUnpack)
  val nearShiftSwap  = cutSwap(nearShiftUnpack)
  val zexNormSwap    = cutSwap(zexNormUnpack)
  val diffSgnSwap    = cutSwap(diffSgnUnpack)
  val xmanSwap       = cutSwapDbg(xmanUnpack)       // dbgPrintf only
  val ymanSwap       = cutSwapDbg(ymanUnpack)       // dbgPrintf only
  val xmanPadXyzSwap = cutSwapDbg(xmanPadXyzUnpack) // dbgPrintf only
  val ymanPadXyzSwap = cutSwapDbg(ymanPadXyzUnpack) // dbgPrintf only
  val xsgnSwap       = cutSwap(xsgnUnpack)
  val ysgnSwap       = cutSwap(ysgnUnpack)
  val xyNaNSwap      = cutSwap(xyNaNUnpack)
  val xyInfSwap      = cutSwap(xyInfUnpack)
  val xyInfSgnSwap   = cutSwap(xyInfSgnUnpack)
  val xyBothZeroSwap = cutSwap(xyBothZeroUnpack)

  //----------------------------------------------------------------------
  // S3 : far path alignment shift + sticky, near path pre-shift
  val yFarR = Cat(yFarSwap, 0.U(2.W))

  // Shift y
  //   Actually, shift amount here is always larger than or equal to 2.
  val shiftW0 = log2Up(maxManWxyz+1+2+1) // +2 for round bits, +1 for Leading 1
  // Note diffExW includes the width for sign bit
  val (shiftW, shiftOut) =
    if   (shiftW0<(diffExW-2)) {(shiftW0, dropLSB(shiftW0,diffFarSwap).orR)}
    else {(diffExW-1, false.B)} // Usually this never happens
  val yFarShift = Mux(shiftOut, 0.U, yFarR >> diffFarSwap(shiftW-1, 0))
  println(f"shiftW = $shiftW%d diffExW = $diffExW%d ${diffFarSwap.getWidth}%d")
  when (!nearSwap) {
    dbgPrintf("xman=%x yman=%x\n",  xmanSwap, ymanSwap)
    dbgPrintf("xmanPadXyz=%x ymanPadXyz=%x\n", xmanPadXyzSwap, ymanPadXyzSwap)
    dbgPrintf("xFar=%x yFarShift=%x\n", xFarSwap, yFarShift)
    //dbgPrintf("%b\n%b\n", yFar, yFarShift)
  }
  val diffFarMinus2    = diffFarSwap(shiftW-1,0) - 2.U
  val diffFarLessThan2 = !diffFarSwap(diffFarSwap.getWidth-1,shiftW-1).orR && diffFarMinus2(shiftW-1)

  // Y shifted = (main value)-round0-round1-stickey
  //    these round0/round1 can be combined to stickey
  //    when zSpec.manW is smaller than x/ySpec.manW.
  //    However, we expect such situations are rare.
  //
  // Get stickey bit
  //   if stickey bit = Or ( all bits [shift-3:0] )
  //   PriorityEncoder in Chisel returns bit position of 1 from lsb
  //     0 returns all 1; unfortunately if bit width of input is 2^w,
  //     2^(w-1) will return the same result as 0.
  val onePosFromLSB = PriorityEncoder(yFarSwap)
  val yFarW = yFarSwap.getWidth
  // Since always a leading 1 bit exists, stickey == 1 when shiftOut, but round should be 0
  val stickeyFar = ( onePosFromLSB < diffFarMinus2) && !diffFarLessThan2
  // Increment for negative value occurs only rounded bits are all 0
  val guardFarNeg = (~Cat(yFarShift(1,0), stickeyFar)) +& 1.U
  val guardFar    = Mux(diffSgnSwap, guardFarNeg(2,1), yFarShift(1,0))
  val negIncFar   = diffSgnSwap && guardFarNeg(3)
  val yAddFar = Cat(diffSgnSwap,
    Mux(diffSgnSwap,~yFarShift(maxManWxyz+2,2),yFarShift(maxManWxyz+2,2)))

  //----------------------------------------------------------------------
  // Near path
  //   when exponent difference <= 1 and sign is different
  //   in this case, assume x is normalized as 1<=x<2 and abs(y)<x,
  //     -1/2 >= y > -2, and 0 < x+y < 1+1/2.
  //   if compare only exponent, the result can be negative.
  //     in this case, -1 < x+y < 1+1/2.
  //
  // Implementation without mantissa comparison - use XFar/YFar
  val xNear = Cat(getMSB(2+maxManWxy, xFarSwap), 0.U(1.W))
  val yNear = getMSB(2+maxManWxy, yFarSwap)

  val yNearShift = Mux( nearShiftSwap, Cat(0.U(2.W), yNear), Cat(0.U(1.W), yNear, 0.U(1.W)) )

  //--------------------------- stage "align" ----------------------------
  private def cutAlign[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("align"))
  }
  private def cutAlignDbg[T <: Data](x: T): T = {
    if (enableDebug) {
      cutAlign(x)
    } else {
      x
    }
  }

  val xFarAlign          = cutAlign(xFarSwap)
  val yAddFarAlign       = cutAlign(yAddFar)
  val negIncFarAlign     = cutAlign(negIncFar)
  val guardFarAlign      = cutAlign(guardFar)
  val stickeyFarAlign    = cutAlign(stickeyFar)
  val xNearAlign         = cutAlign(xNear)
  val yNearShiftAlign    = cutAlign(yNearShift)
  val nearAlign          = cutAlign(nearSwap)
  val xFarSignAlign      = cutAlign(xFarSignSwap)
  val zexNormAlign       = cutAlign(zexNormSwap)
  val xsgnAlign          = cutAlign(xsgnSwap)
  val ysgnAlign          = cutAlign(ysgnSwap)
  val xyNaNAlign         = cutAlign(xyNaNSwap)
  val xyInfAlign         = cutAlign(xyInfSwap)
  val xyInfSgnAlign      = cutAlign(xyInfSgnSwap)
  val xyBothZeroAlign    = cutAlign(xyBothZeroSwap)
  val yNearAlign         = cutAlignDbg(yNear)         // dbgPrintf only
  val nearShiftAlign     = cutAlignDbg(nearShiftSwap) // dbgPrintf only
  val yFarAlign          = cutAlignDbg(yFarSwap)      // dbgPrintf only
  val onePosFromLSBAlign = cutAlignDbg(onePosFromLSB) // dbgPrintf only
  val diffFarAlign       = cutAlignDbg(diffFarSwap)   // dbgPrintf only
  val shiftOutAlign      = cutAlignDbg(shiftOut)      // dbgPrintf only
  val diffSgnAlign       = cutAlignDbg(diffSgnSwap)   // dbgPrintf only

  //----------------------------------------------------------------------
  // S4 : significand add (far) / subtract (near)
  // Width = maxManWxyz+2
  val sumFar = xFarAlign + yAddFarAlign + negIncFarAlign
  when (!nearAlign) {
    dbgPrintf("yFar=%x onePosFromLSB=%d diffFar=%d(%x) shiftOut=%d\n",
      yFarAlign, onePosFromLSBAlign, diffFarAlign, diffFarAlign, shiftOutAlign)
    dbgPrintf("guardFar=%x stickeyFar=%d\n", guardFarAlign, stickeyFarAlign)
    dbgPrintf("diffSgn=%d sumFar=%x yAddFar=%x negIncFar=%b\n",
      diffSgnAlign, sumFar, yAddFarAlign, negIncFarAlign)
  }

  val sumNear = xNearAlign - yNearShiftAlign
  // The below includes 1.xxxx+round bit
  //  width = maxManWxy+2
  val sumNearSign = getMSB1(0,sumNear)
  val sumNearAbs = Mux(sumNearSign, -sumNear(maxManWxy+1,0), sumNear(maxManWxy+1,0))
  val sumNearZero = !sumNear.orR

  // leading zero anticipation, in parallel with the subtraction above
  // (TODO: carry correction does not run in parallel)
  val shiftApprox = if(useParallelLZA) {
    val sumNearG  =   xNearAlign  & (~yNearShiftAlign) // x > y  (Greater)
    val sumNearS  = (~xNearAlign) &   yNearShiftAlign  // x < y  (Smaller)
    val sumNearE  = ~(sumNearG | sumNearS)       // x == y (otherwise)
    val sumNearF  = VecInit(Seq.tabulate(sumNearE.getWidth)( i => {
        val eNext = if(i+1 < sumNearE.getWidth) {sumNearE(i+1)} else {1.U(1.W)}
        val sPrev = if(0 <= i-1) {sumNearS(i-1)} else {0.U(1.W)}
        val gPrev = if(0 <= i-1) {sumNearG(i-1)} else {0.U(1.W)}
        val sCurr = sumNearS(i)
        val gCurr = sumNearG(i)
        val term1 =   eNext  & ((gCurr & (~sPrev)) | (sCurr & (~gPrev)))
        val term2 = (~eNext) & ((sCurr & (~sPrev)) | (gCurr & (~gPrev)))
        term1 | term2
      })).asUInt

    assert(sumNearF.getWidth == sumNearE.getWidth)

    Some(PriorityEncoder(Reverse(sumNearF)))
  } else {
    None
  }

  //---------------------------- stage "add" -----------------------------
  private def cutAdd[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("add"))
  }
  private def cutAddDbg[T <: Data](x: T): T = {
    if (enableDebug) {
      cutAdd(x)
    } else {
      x
    }
  }

  val sumFarAdd      = cutAdd(sumFar)
  val guardFarAdd    = cutAdd(guardFarAlign)
  val stickeyFarAdd  = cutAdd(stickeyFarAlign)
  val sumNearAdd     = cutAddDbg(sumNear)         // dbgPrintf only
  val sumNearAbsAdd  = cutAdd(sumNearAbs)
  val sumNearSignAdd = cutAdd(sumNearSign)
  val sumNearZeroAdd = cutAdd(sumNearZero)
  val shiftApproxAdd = shiftApprox.map(cutAdd(_))
  val xNearAdd       = cutAddDbg(xNearAlign)      // dbgPrintf only
  val yNearAdd       = cutAddDbg(yNearAlign)      // dbgPrintf only
  val nearShiftAdd   = cutAddDbg(nearShiftAlign)  // dbgPrintf only
  val yNearShiftAdd  = cutAddDbg(yNearShiftAlign) // dbgPrintf only
  val nearAdd        = cutAdd(nearAlign)
  val xFarSignAdd    = cutAdd(xFarSignAlign)
  val zexNormAdd     = cutAdd(zexNormAlign)
  val xsgnAdd        = cutAdd(xsgnAlign)
  val ysgnAdd        = cutAdd(ysgnAlign)
  val xyNaNAdd       = cutAdd(xyNaNAlign)
  val xyInfAdd       = cutAdd(xyInfAlign)
  val xyInfSgnAdd    = cutAdd(xyInfSgnAlign)
  val xyBothZeroAdd  = cutAdd(xyBothZeroAlign)

  //----------------------------------------------------------------------
  // S5 : normalisation (near), mantissa / round selection (far), merge
  val sumFarWithGuard = Cat(sumFarAdd, guardFarAdd)

  // Rounding
  val zmanFar     = Wire(UInt(zSpec.manW.W)) // This omit leading zero
  val roundFarSel = Wire(Bool())
  val additionalStickeyFar = Wire(Bool())
  val exSubFar    = Wire(SInt(2.W))
  val diffWz = maxManWxyz-zSpec.manW // >=0
  when (sumFarAdd(maxManWxyz+1).asBool) { // >=2
    zmanFar      := sumFarAdd(maxManWxyz, diffWz+1)
    roundFarSel  := sumFarAdd(diffWz)
    additionalStickeyFar := orRLSB(diffWz+2, sumFarWithGuard)
    exSubFar     := -1.S
  }.elsewhen (sumFarAdd(maxManWxyz)===1.U) { // >=1
    zmanFar      := sumFarAdd(maxManWxyz-1, diffWz)
    roundFarSel  := sumFarWithGuard(diffWz+1)
    additionalStickeyFar := orRLSB(diffWz+1, sumFarWithGuard)
    exSubFar     := 0.S
  }.otherwise { // >=1/2
    zmanFar      := sumFarWithGuard(maxManWxyz, diffWz+1)
    roundFarSel  := sumFarWithGuard(diffWz)
    additionalStickeyFar := orRLSB(diffWz, sumFarWithGuard)
    exSubFar     := 1.S
  }
  when (!nearAdd) { dbgPrintf("zmanFar=%x\n", zmanFar) }

  val incFar = FloatChiselUtil.roundIncBySpec(roundSpec, zmanFar(0), roundFarSel,
    stickeyFarAdd || additionalStickeyFar)

  // Near path normalisation
  val (shiftSum, sumNearNorm) = if(useParallelLZA) {
    val shiftApproxR = shiftApproxAdd.get
    val sumNearNorm0 = (sumNearAbsAdd << shiftApproxR)(maxManWxy+2,0)
    val sumNearNorm  = Mux(sumNearNorm0(maxManWxy+2) === 1.U, // carry correction
      sumNearNorm0(maxManWxy+1,1), sumNearNorm0(maxManWxy,0))
    val shiftSum = shiftApproxR - sumNearNorm0(maxManWxy+2)

    if(enableDebug) {
      val shiftRef = PriorityEncoder(Reverse(sumNearAbsAdd))
      assert(shiftRef === shiftSum)
    }

    (shiftSum, sumNearNorm)

  } else {
    // 0 .. >=1
    // 1 .. >=1/2 ...
    val shiftSum    = PriorityEncoder(Reverse(sumNearAbsAdd))
    val sumNearNorm = (sumNearAbsAdd << shiftSum)(maxManWxy,0) // include (possible) round bit
    (shiftSum, sumNearNorm)
  }

  val nearSign = (!sumNearZeroAdd) && (sumNearSignAdd ^ xFarSignAdd)

  // Rounding
  val zmanNear = Wire(UInt((zSpec.manW).W))
  val incNear  = Wire(Bool())
  if (zSpec.manW < maxManWxy) { // Output precision is lower
    val roundNear   = getMSB1(zSpec.manW+1, sumNearNorm)
    val stickeyNear = orRLSB(maxManWxy-zSpec.manW, sumNearNorm) // This is usually zero
    incNear  := FloatChiselUtil.roundIncBySpec(roundSpec, getMSB1(zSpec.manW,sumNearNorm), roundNear,
      stickeyNear)
    zmanNear := sumNearNorm(maxManWxy, maxManWxy-zSpec.manW+1)
  } else if (zSpec.manW == maxManWxy) { // Output precision is the same as input
    incNear  := FloatChiselUtil.roundIncBySpec(roundSpec, sumNearNorm(1), sumNearNorm(0), false.B)
    zmanNear := sumNearNorm(maxManWxy, 1)
  } else { // Output precision is higher, no rounding
    incNear  := 0.U(1.W)
    zmanNear := padLSB( zSpec.manW, sumNearNorm(maxManWxy, 0) )
  }
  when (nearAdd) {
    dbgPrintf("xNear=%b yNear=%b nearShift=%d yNearShift=%b\n",
      xNearAdd, yNearAdd, nearShiftAdd, yNearShiftAdd)
    dbgPrintf("sumNear=%b sumNearAbs=%b\n", sumNearAdd, sumNearAbsAdd)
    dbgPrintf("sumNearAbs=%b\n", sumNearAbsAdd)
    dbgPrintf("shiftSum=%d(%b) sumNearNorm=%b incNear=%b\n", shiftSum, shiftSum, sumNearNorm, incNear)
  }

  //----------------------------------------------------------------------
  // Merge Far / Near
  val zmanSel = Mux(nearAdd, zmanNear, zmanFar)
  val incSel  = Mux(nearAdd, incNear, incFar)
  val zeroSel = nearAdd && sumNearZeroAdd
  val exSub   = Mux(nearAdd, shiftSum.zext, exSubFar) // SInt
  val zsgn0   = Mux(nearAdd, nearSign, xFarSignAdd)

  //---------------------------- stage "norm" ----------------------------
  private def cutNorm[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("norm"))
  }

  val zmanSelNorm    = cutNorm(zmanSel)
  val incSelNorm     = cutNorm(incSel)
  val zeroSelNorm    = cutNorm(zeroSel)
  val exSubNorm      = cutNorm(exSub)
  val zsgn0Norm      = cutNorm(zsgn0)
  val zexNormNorm    = cutNorm(zexNormAdd)
  val xsgnNorm       = cutNorm(xsgnAdd)
  val ysgnNorm       = cutNorm(ysgnAdd)
  val xyNaNNorm      = cutNorm(xyNaNAdd)
  val xyInfNorm      = cutNorm(xyInfAdd)
  val xyInfSgnNorm   = cutNorm(xyInfSgnAdd)
  val xyBothZeroNorm = cutNorm(xyBothZeroAdd)

  //----------------------------------------------------------------------
  // S6 : rounding increment and exponent adjustment
  val zman0   = zmanSelNorm +& incSelNorm

  // Calculate zex Range
  val zexMax = zexBaseMax+2
  // In cases when zSpec.manW is smaller, it will increase 2 (result become 4)
  val zexMin = zexBaseMin - sumNearAbs.getWidth
  val zexW   = getWidthToRepresentNumbersAlwaysSigned( Seq(zexMin,  zexMax) )
  val zex0   = zexNormNorm.pad(zexW) - exSubNorm + zman0(zSpec.manW).zext
  val zeroAfterAdd = if(zSpec.disableSubnormal){
    (zex0 <= 0.S) || zeroSelNorm // if subnormal is disabled, include ex==0 as zero.
  } else {
    (zex0 <  0.S) || zeroSelNorm
  }
  val infAfterAdd = zex0 >= maskI(zSpec.exW).S

  dbgPrintf("zex0 = %d, zeroAfterAdd = %d\n", zex0, zeroAfterAdd)

  //--------------------------- stage "round" ----------------------------
  private def cutRound[T <: Data](x: T): T = {
    ShiftRegister(x, stages.at("round"))
  }

  val zmanBitsRound     = cutRound(zman0(zSpec.manW-1, 0))
  val zexBitsRound      = cutRound(zex0(zSpec.exW-1, 0))
  val zeroAfterAddRound = cutRound(zeroAfterAdd)
  val infAfterAddRound  = cutRound(infAfterAdd)
  val zsgn0Round        = cutRound(zsgn0Norm)
  val xsgnRound         = cutRound(xsgnNorm)
  val ysgnRound         = cutRound(ysgnNorm)
  val xyNaNRound        = cutRound(xyNaNNorm)
  val xyInfRound        = cutRound(xyInfNorm)
  val xyInfSgnRound     = cutRound(xyInfSgnNorm)
  val xyBothZeroRound   = cutRound(xyBothZeroNorm)

  //----------------------------------------------------------------------
  // S7 : Final (special case muxes)
  val infOrNaN = xyInfRound || xyNaNRound || infAfterAddRound
  val zero     = zeroAfterAddRound||xyBothZeroRound

  val zex = Mux (infOrNaN, maskU(zSpec.exW),
    Mux (zero, 0.U(zSpec.exW.W), zexBitsRound))
  val zman = Mux(infOrNaN || zero, Cat(xyNaNRound, 0.U((zSpec.manW-1).W)), zmanBitsRound)
  val zsgn = if(roundSpec == RoundSpec.truncate) {
    Mux(xyNaNRound, 0.U(1.W), Mux(xyBothZeroRound, xsgnRound | ysgnRound, Mux(xyInfRound, xyInfSgnRound, zsgn0Round)))
  } else {
    Mux(xyNaNRound, 0.U(1.W), Mux(xyBothZeroRound, xsgnRound & ysgnRound, Mux(xyInfRound, xyInfSgnRound, zsgn0Round)))
  }
  // currently it does not distinguish positive NaN and negative NaN

  val z0 = Cat(zsgn, zex, zman)

  io.z   := ShiftRegister(z0, stages.at("out"))
  dbgPrintf("x=%x y=%x z=%x\n", io.x, io.y, io.z)
}

/** Pipeline stage points of [[rial.arith.AddFPGeneric]].
 *
 * @see [[rial.util.StagePoint]]
 */
object AddFPGeneric {
  /** Ordered boundaries of the adder datapath, with a rough,
   *  target-independent estimate of the delay of the logic ending at each. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("unpack", 3), StagePoint("swap",  1), StagePoint("align", 3),
    StagePoint("add",    3), StagePoint("norm",  2), StagePoint("round", 4),
    StagePoint("out",    1))
}

/** Specialization of AddFPGeneric for FP64.
 *
 * @see [[rial.arith.AddFPGeneric]]
 */
class AddFP64( stage : PipelineStageConfig )
    extends AddFPGeneric( RealSpec.Float64Spec,
      RealSpec.Float64Spec,  RealSpec.Float64Spec,
      RoundSpec.roundToEven, stage) {
}
