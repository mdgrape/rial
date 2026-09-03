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

  val xmanPadXyz = Mux(xzero, 0.U((maxManWxyz+1).W), 1.U(1.W) ## padLSB(maxManWxyz, xman))
  val ymanPadXyz = Mux(yzero, 0.U((maxManWxyz+1).W), 1.U(1.W) ## padLSB(maxManWxyz, yman))

  val yzeroFar = Mux(exXlargerThanY, xzero, yzero)

  //---------------------------- stage "unpack" --------------------------
  private def cutUnpack[T <: Data](x: T): T = ShiftRegister(x, stages.at("unpack"))
  // Signals read by `dbgPrintf` only must not create registers in a normal
  // (non-debug) build, so they are piped only when the printfs exist.
  private def cutUnpackDbg[T <: Data](x: T): T = if (enableDebug) cutUnpack(x) else x

  val xmanPadXyz_1   = cutUnpack(xmanPadXyz)
  val ymanPadXyz_1   = cutUnpack(ymanPadXyz)
  val exXlargerThanY_1 = cutUnpack(exXlargerThanY)
  val diffExXY_1     = cutUnpack(diffExXY)
  val diffExYX_1     = cutUnpack(diffExYX)
  val near_1         = cutUnpack(near)
  val nearShift_1    = cutUnpack(nearShift)
  val zexNorm_1      = cutUnpack(zexNorm)
  val diffSgn_1      = cutUnpack(diffSgn)
  val xman_1         = cutUnpackDbg(xman)   // dbgPrintf only
  val yman_1         = cutUnpackDbg(yman)   // dbgPrintf only
  val xsgn_1         = cutUnpack(xsgn)
  val ysgn_1         = cutUnpack(ysgn)
  val xyNaN_1        = cutUnpack(xyNaN)
  val xyInf_1        = cutUnpack(xyInf)
  val xyInfSgn_1     = cutUnpack(xyInfSgn)
  val xyBothZero_1   = cutUnpack(xyBothZero)

  //----------------------------------------------------------------------
  // S2 : operand swap (larger exponent -> xFar)
  //   Far path
  //   when exponent difference > 1 or sign is the same
  //   in this case, assume x is normalized as 1<=x<2 and abs(y)<x,
  //     -1/2 < y < 2, and then 1/2 < x+y < 4.
  //
  //   currently, separate barrel shifters/adders are used for far/near paths.
  val xFar = 0.U(1.W) ## Mux(exXlargerThanY_1, xmanPadXyz_1, ymanPadXyz_1)
  val yFar = Mux(exXlargerThanY_1, ymanPadXyz_1, xmanPadXyz_1)
  val xFarSign = Mux(exXlargerThanY_1, xsgn_1, ysgn_1).asBool
  // Drop sign bit
  val diffFar = Mux(exXlargerThanY_1, diffExXY_1(diffExW-2,0), diffExYX_1(diffExW-2,0))

  //---------------------------- stage "swap" ----------------------------
  private def cutSwap[T <: Data](x: T): T = ShiftRegister(x, stages.at("swap"))
  private def cutSwapDbg[T <: Data](x: T): T = if (enableDebug) cutSwap(x) else x

  val xFar_2       = cutSwap(xFar)
  val yFar_2       = cutSwap(yFar)
  val xFarSign_2   = cutSwap(xFarSign)
  val diffFar_2    = cutSwap(diffFar)
  val near_2       = cutSwap(near_1)
  val nearShift_2  = cutSwap(nearShift_1)
  val zexNorm_2    = cutSwap(zexNorm_1)
  val diffSgn_2    = cutSwap(diffSgn_1)
  val xman_2       = cutSwapDbg(xman_1)     // dbgPrintf only
  val yman_2       = cutSwapDbg(yman_1)     // dbgPrintf only
  val xmanPadXyz_2 = cutSwapDbg(xmanPadXyz_1) // dbgPrintf only
  val ymanPadXyz_2 = cutSwapDbg(ymanPadXyz_1) // dbgPrintf only
  val xsgn_2       = cutSwap(xsgn_1)
  val ysgn_2       = cutSwap(ysgn_1)
  val xyNaN_2      = cutSwap(xyNaN_1)
  val xyInf_2      = cutSwap(xyInf_1)
  val xyInfSgn_2   = cutSwap(xyInfSgn_1)
  val xyBothZero_2 = cutSwap(xyBothZero_1)

  //----------------------------------------------------------------------
  // S3 : far path alignment shift + sticky, near path pre-shift
  val yFarR = yFar_2 ## 0.U(2.W)

  // Shift y
  //   Actually, shift amount here is always larger than or equal to 2.
  val shiftW0 = log2Up(maxManWxyz+1+2+1) // +2 for round bits, +1 for Leading 1
  // Note diffExW includes the width for sign bit
  val (shiftW, shiftOut) =
    if   (shiftW0<(diffExW-2)) {(shiftW0, dropLSB(shiftW0,diffFar_2).orR)}
    else {(diffExW-1, false.B)} // Usually this never happens
  val yFarShift = Mux(shiftOut, 0.U, yFarR >> diffFar_2(shiftW-1, 0))
  println(f"shiftW = $shiftW%d diffExW = $diffExW%d ${diffFar_2.getWidth}%d")
  when (!near_2) {
    dbgPrintf("xman=%x yman=%x\n",  xman_2, yman_2)
    dbgPrintf("xmanPadXyz=%x ymanPadXyz=%x\n", xmanPadXyz_2, ymanPadXyz_2)
    dbgPrintf("xFar=%x yFarShift=%x\n", xFar_2, yFarShift)
    //dbgPrintf("%b\n%b\n", yFar, yFarShift)
  }
  val diffFarMinus2    = diffFar_2(shiftW-1,0) - 2.U
  val diffFarLessThan2 = !diffFar_2(diffFar_2.getWidth-1,shiftW-1).orR && diffFarMinus2(shiftW-1)

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
  val onePosFromLSB = PriorityEncoder(yFar_2)
  val yFarW = yFar_2.getWidth
  // Since always a leading 1 bit exists, stickey == 1 when shiftOut, but round should be 0
  val stickeyFar = ( onePosFromLSB < diffFarMinus2) && !diffFarLessThan2
  // Increment for negative value occurs only rounded bits are all 0
  val guardFarNeg = (~(yFarShift(1,0) ## stickeyFar)) +& 1.U
  val guardFar    = Mux(diffSgn_2, guardFarNeg(2,1), yFarShift(1,0))
  val negIncFar   = diffSgn_2 && guardFarNeg(3)
  val yAddFar = diffSgn_2 ## Mux(diffSgn_2,~yFarShift(maxManWxyz+2,2),yFarShift(maxManWxyz+2,2))

  //----------------------------------------------------------------------
  // Near path
  //   when exponent difference <= 1 and sign is different
  //   in this case, assume x is normalized as 1<=x<2 and abs(y)<x,
  //     -1/2 >= y > -2, and 0 < x+y < 1+1/2.
  //   if compare only exponent, the result can be negative.
  //     in this case, -1 < x+y < 1+1/2.
  //
  // Implementation without mantissa comparison - use XFar/YFar
  val xNear = getMSB(2+maxManWxy, xFar_2) ## 0.U(1.W)
  val yNear = getMSB(2+maxManWxy, yFar_2)

  val yNearShift = Mux( nearShift_2, 0.U(2.W) ## yNear, 0.U(1.W) ## yNear ## 0.U(1.W) )

  //--------------------------- stage "align" ----------------------------
  private def cutAlign[T <: Data](x: T): T = ShiftRegister(x, stages.at("align"))
  private def cutAlignDbg[T <: Data](x: T): T = if (enableDebug) cutAlign(x) else x

  val xFar_3       = cutAlign(xFar_2)
  val yAddFar_3    = cutAlign(yAddFar)
  val negIncFar_3  = cutAlign(negIncFar)
  val guardFar_3   = cutAlign(guardFar)
  val stickeyFar_3 = cutAlign(stickeyFar)
  val xNear_3      = cutAlign(xNear)
  val yNearShift_3 = cutAlign(yNearShift)
  val near_3       = cutAlign(near_2)
  val xFarSign_3   = cutAlign(xFarSign_2)
  val zexNorm_3    = cutAlign(zexNorm_2)
  val xsgn_3       = cutAlign(xsgn_2)
  val ysgn_3       = cutAlign(ysgn_2)
  val xyNaN_3      = cutAlign(xyNaN_2)
  val xyInf_3      = cutAlign(xyInf_2)
  val xyInfSgn_3   = cutAlign(xyInfSgn_2)
  val xyBothZero_3 = cutAlign(xyBothZero_2)
  val yNear_3         = cutAlignDbg(yNear)           // dbgPrintf only
  val nearShift_3     = cutAlignDbg(nearShift_2)     // dbgPrintf only
  val yFar_3          = cutAlignDbg(yFar_2)          // dbgPrintf only
  val onePosFromLSB_3 = cutAlignDbg(onePosFromLSB)   // dbgPrintf only
  val diffFar_3       = cutAlignDbg(diffFar_2)       // dbgPrintf only
  val shiftOut_3      = cutAlignDbg(shiftOut)        // dbgPrintf only
  val diffSgn_3       = cutAlignDbg(diffSgn_2)       // dbgPrintf only

  //----------------------------------------------------------------------
  // S4 : significand add (far) / subtract (near)
  // Width = maxManWxyz+2
  val sumFar = xFar_3 + yAddFar_3 + negIncFar_3
  when (!near_3) {
    dbgPrintf("yFar=%x onePosFromLSB=%d diffFar=%d(%x) shiftOut=%d\n",
      yFar_3, onePosFromLSB_3, diffFar_3, diffFar_3, shiftOut_3)
    dbgPrintf("guardFar=%x stickeyFar=%d\n", guardFar_3, stickeyFar_3)
    dbgPrintf("diffSgn=%d sumFar=%x yAddFar=%x negIncFar=%b\n",
      diffSgn_3, sumFar, yAddFar_3, negIncFar_3)
  }

  val sumNear = xNear_3 - yNearShift_3
  // The below includes 1.xxxx+round bit
  //  width = maxManWxy+2
  val sumNearSign = getMSB1(0,sumNear)
  val sumNearAbs = Mux(sumNearSign, -sumNear(maxManWxy+1,0), sumNear(maxManWxy+1,0))
  val sumNearZero = !sumNear.orR

  // leading zero anticipation, in parallel with the subtraction above
  // (TODO: carry correction does not run in parallel)
  val shiftApprox = if(useParallelLZA) {
    val sumNearG  =   xNear_3  & (~yNearShift_3) // x > y  (Greater)
    val sumNearS  = (~xNear_3) &   yNearShift_3  // x < y  (Smaller)
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
  private def cutAdd[T <: Data](x: T): T = ShiftRegister(x, stages.at("add"))
  private def cutAddDbg[T <: Data](x: T): T = if (enableDebug) cutAdd(x) else x

  val sumFar_4       = cutAdd(sumFar)
  val guardFar_4     = cutAdd(guardFar_3)
  val stickeyFar_4   = cutAdd(stickeyFar_3)
  val sumNear_4      = cutAddDbg(sumNear)      // dbgPrintf only
  val sumNearAbs_4   = cutAdd(sumNearAbs)
  val sumNearSign_4  = cutAdd(sumNearSign)
  val sumNearZero_4  = cutAdd(sumNearZero)
  val shiftApprox_4  = shiftApprox.map(cutAdd(_))
  val xNear_4        = cutAddDbg(xNear_3)      // dbgPrintf only
  val yNear_4        = cutAddDbg(yNear_3)      // dbgPrintf only
  val nearShift_4    = cutAddDbg(nearShift_3)  // dbgPrintf only
  val yNearShift_4   = cutAddDbg(yNearShift_3) // dbgPrintf only
  val near_4         = cutAdd(near_3)
  val xFarSign_4     = cutAdd(xFarSign_3)
  val zexNorm_4      = cutAdd(zexNorm_3)
  val xsgn_4         = cutAdd(xsgn_3)
  val ysgn_4         = cutAdd(ysgn_3)
  val xyNaN_4        = cutAdd(xyNaN_3)
  val xyInf_4        = cutAdd(xyInf_3)
  val xyInfSgn_4     = cutAdd(xyInfSgn_3)
  val xyBothZero_4   = cutAdd(xyBothZero_3)

  //----------------------------------------------------------------------
  // S5 : normalisation (near), mantissa / round selection (far), merge
  val sumFarWithGuard = sumFar_4 ## guardFar_4

  // Rounding
  val zmanFar     = Wire(UInt(zSpec.manW.W)) // This omit leading zero
  val roundFarSel = Wire(Bool())
  val additionalStickeyFar = Wire(Bool())
  val exSubFar    = Wire(SInt(2.W))
  val diffWz = maxManWxyz-zSpec.manW // >=0
  when (sumFar_4(maxManWxyz+1).asBool) { // >=2
    zmanFar      := sumFar_4(maxManWxyz, diffWz+1)
    roundFarSel  := sumFar_4(diffWz)
    additionalStickeyFar := orRLSB(diffWz+2, sumFarWithGuard)
    exSubFar     := -1.S
  }.elsewhen (sumFar_4(maxManWxyz)===1.U) { // >=1
    zmanFar      := sumFar_4(maxManWxyz-1, diffWz)
    roundFarSel  := sumFarWithGuard(diffWz+1)
    additionalStickeyFar := orRLSB(diffWz+1, sumFarWithGuard)
    exSubFar     := 0.S
  }.otherwise { // >=1/2
    zmanFar      := sumFarWithGuard(maxManWxyz, diffWz+1)
    roundFarSel  := sumFarWithGuard(diffWz)
    additionalStickeyFar := orRLSB(diffWz, sumFarWithGuard)
    exSubFar     := 1.S
  }
  when (!near_4) { dbgPrintf("zmanFar=%x\n", zmanFar) }

  val incFar = FloatChiselUtil.roundIncBySpec(roundSpec, zmanFar(0), roundFarSel,
    stickeyFar_4 || additionalStickeyFar)

  // Near path normalisation
  val (shiftSum, sumNearNorm) = if(useParallelLZA) {
    val shiftApproxR = shiftApprox_4.get
    val sumNearNorm0 = (sumNearAbs_4 << shiftApproxR)(maxManWxy+2,0)
    val sumNearNorm  = Mux(sumNearNorm0(maxManWxy+2) === 1.U, // carry correction
      sumNearNorm0(maxManWxy+1,1), sumNearNorm0(maxManWxy,0))
    val shiftSum = shiftApproxR - sumNearNorm0(maxManWxy+2)

    if(enableDebug) {
      val shiftRef = PriorityEncoder(Reverse(sumNearAbs_4))
      assert(shiftRef === shiftSum)
    }

    (shiftSum, sumNearNorm)

  } else {
    // 0 .. >=1
    // 1 .. >=1/2 ...
    val shiftSum    = PriorityEncoder(Reverse(sumNearAbs_4))
    val sumNearNorm = (sumNearAbs_4 << shiftSum)(maxManWxy,0) // include (possible) round bit
    (shiftSum, sumNearNorm)
  }

  val nearSign = (!sumNearZero_4) && (sumNearSign_4 ^ xFarSign_4)

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
  when (near_4) {
    dbgPrintf("xNear=%b yNear=%b nearShift=%d yNearShift=%b\n",
      xNear_4, yNear_4, nearShift_4, yNearShift_4)
    dbgPrintf("sumNear=%b sumNearAbs=%b\n", sumNear_4, sumNearAbs_4)
    dbgPrintf("sumNearAbs=%b\n", sumNearAbs_4)
    dbgPrintf("shiftSum=%d(%b) sumNearNorm=%b incNear=%b\n", shiftSum, shiftSum, sumNearNorm, incNear)
  }

  //----------------------------------------------------------------------
  // Merge Far / Near
  val zmanSel = Mux(near_4, zmanNear, zmanFar)
  val incSel  = Mux(near_4, incNear, incFar)
  val zeroSel = near_4 && sumNearZero_4
  val exSub   = Mux(near_4, shiftSum.zext, exSubFar) // SInt
  val zsgn0   = Mux(near_4, nearSign, xFarSign_4)

  //---------------------------- stage "norm" ----------------------------
  private def cutNorm[T <: Data](x: T): T = ShiftRegister(x, stages.at("norm"))

  val zmanSel_5    = cutNorm(zmanSel)
  val incSel_5     = cutNorm(incSel)
  val zeroSel_5    = cutNorm(zeroSel)
  val exSub_5      = cutNorm(exSub)
  val zsgn0_5      = cutNorm(zsgn0)
  val zexNorm_5    = cutNorm(zexNorm_4)
  val xsgn_5       = cutNorm(xsgn_4)
  val ysgn_5       = cutNorm(ysgn_4)
  val xyNaN_5      = cutNorm(xyNaN_4)
  val xyInf_5      = cutNorm(xyInf_4)
  val xyInfSgn_5   = cutNorm(xyInfSgn_4)
  val xyBothZero_5 = cutNorm(xyBothZero_4)

  //----------------------------------------------------------------------
  // S6 : rounding increment and exponent adjustment
  val zman0   = zmanSel_5 +& incSel_5

  // Calculate zex Range
  val zexMax = zexBaseMax+2
  // In cases when zSpec.manW is smaller, it will increase 2 (result become 4)
  val zexMin = zexBaseMin - sumNearAbs.getWidth
  val zexW   = getWidthToRepresentNumbersAlwaysSigned( Seq(zexMin,  zexMax) )
  val zex0   = zexNorm_5.pad(zexW) - exSub_5 + zman0(zSpec.manW).zext
  val zeroAfterAdd = if(zSpec.disableSubnormal){
    (zex0 <= 0.S) || zeroSel_5 // if subnormal is disabled, include ex==0 as zero.
  } else {
    (zex0 <  0.S) || zeroSel_5
  }
  val infAfterAdd = zex0 >= maskI(zSpec.exW).S

  dbgPrintf("zex0 = %d, zeroAfterAdd = %d\n", zex0, zeroAfterAdd)

  //--------------------------- stage "round" ----------------------------
  private def cutRound[T <: Data](x: T): T = ShiftRegister(x, stages.at("round"))

  val zmanBits_6     = cutRound(zman0(zSpec.manW-1, 0))
  val zexBits_6      = cutRound(zex0(zSpec.exW-1, 0))
  val zeroAfterAdd_6 = cutRound(zeroAfterAdd)
  val infAfterAdd_6  = cutRound(infAfterAdd)
  val zsgn0_6        = cutRound(zsgn0_5)
  val xsgn_6         = cutRound(xsgn_5)
  val ysgn_6         = cutRound(ysgn_5)
  val xyNaN_6        = cutRound(xyNaN_5)
  val xyInf_6        = cutRound(xyInf_5)
  val xyInfSgn_6     = cutRound(xyInfSgn_5)
  val xyBothZero_6   = cutRound(xyBothZero_5)

  //----------------------------------------------------------------------
  // S7 : Final (special case muxes)
  val infOrNaN = xyInf_6 || xyNaN_6 || infAfterAdd_6
  val zero     = zeroAfterAdd_6||xyBothZero_6

  val zex = Mux (infOrNaN, maskU(zSpec.exW),
    Mux (zero, 0.U(zSpec.exW.W), zexBits_6))
  val zman = Mux(infOrNaN || zero, xyNaN_6 ## 0.U((zSpec.manW-1).W), zmanBits_6)
  val zsgn = if(roundSpec == RoundSpec.truncate) {
    Mux(xyNaN_6, 0.U(1.W), Mux(xyBothZero_6, xsgn_6 | ysgn_6, Mux(xyInf_6, xyInfSgn_6, zsgn0_6)))
  } else {
    Mux(xyNaN_6, 0.U(1.W), Mux(xyBothZero_6, xsgn_6 & ysgn_6, Mux(xyInf_6, xyInfSgn_6, zsgn0_6)))
  }
  // currently it does not distinguish positive NaN and negative NaN

  val z0 = zsgn ## zex ## zman

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
