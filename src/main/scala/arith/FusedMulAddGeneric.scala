//% @file FusedMulAddFPGeneric.scala
//
// Generic Floating-Point FMA Unit
// Copyright (C) Toru Niina RIKEN BDR 2021
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

/** Generic Floating-Point Fused Multiply-Adder.
 *
 * It takes `x`, `y`, `z` and returns `w` as `w = x * y + z`.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`     the module input
 *  - `"unpack"` after the operand unpacking
 *  - `"mul"`    after the mantissa product
 *  - `"align"`  after the alignment shift
 *  - `"add"`    after the significand add
 *  - `"norm"`   after the normalisation
 *  - `"round"`  after the rounding
 *  - `"out"`    the module output
 *
 * [[getStage]] gives the resulting latency, and `out.valid` is `in.valid`
 * delayed by that many cycles.
 *
 * {{{
 * // three registers, placed by the module
 * new FusedMulAddFPGeneric(spec, spec, spec, spec, RoundSpec.roundToEven, 3)
 * // hand-picked positions
 * new FusedMulAddFPGeneric(spec, spec, spec, spec, RoundSpec.roundToEven,
 *   PipelineStageConfig.specify("mul" -> 1, "add" -> 1, "round" -> 1))
 * // all the registers at the module output
 * new FusedMulAddFPGeneric(spec, spec, spec, spec, RoundSpec.roundToEven,
 *   PipelineStageConfig.atOut(3))
 * // combinational
 * new FusedMulAddFPGeneric(spec, spec, spec, spec, RoundSpec.roundToEven,
 *   PipelineStageConfig.none)
 * }}}
 *
 * {{{
 * class FusedMulAddFPGeneric(...) extends Module {
 *   val io = IO(new Bundle{
 *     val in  = Flipped(Valid(new FusedMulAddFPGeneric.Input(xSpec, ySpec, zSpec)))
 *     val out = Valid(UInt(wSpec.W.W))
 *   })
 *   //...
 * }
 * }}}
 *
 * @constructor create a chisel Module FusedMulAddFPGeneric
 * @param xSpec          floating point spec of the first input
 * @param ySpec          floating point spec of the second input
 * @param zSpec          floating point spec of the third input
 * @param wSpec          floating point spec of the output
 * @param roundSpec      rounding algorithm
 * @param stage          pipeline registers of this module.
 * @param enableDebug
 *   add some additional printf and assertions. by default, false.
 *
 * For parameters, see also:
 * @see [[rial.arith.RealSpec]]
 * @see [[rial.arith.RoundSpec]]
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 *
 */
class FusedMulAddFPGeneric(
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, wSpec : RealSpec,
  roundSpec : RoundSpec,
  stage : PipelineStageConfig,
  val enableDebug : Boolean = false
) extends Module with DebugControl {

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = FusedMulAddFPGeneric.stagePoints

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getParam = { (xSpec, ySpec, zSpec, wSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new FusedMulAddFPGeneric.Input(xSpec, ySpec, zSpec)))
    val out = Valid(UInt(wSpec.W.W))
  })

  // The far-product / far-addend / near paths run in parallel, so every
  // boundary but "in" and "out" is crossed by each of them: the `cut<Point>`
  // helper of a boundary is used once per path, plus once for the control
  // signals the three paths share. A signal that is only read by a
  // `dbgPrintf` must not create a register in a normal (non-debug) build, so
  // it is piped with the `cut<Point>Dbg` variant, which is the identity when
  // `enableDebug` is false.

  //----------------------------------------------------------------------
  // S0 : input registers ("in")
  private val cutIn = new SharedPipe("in", io.in.valid, stages.at("in"))

  val xIn = cutIn(io.in.bits.x)
  val yIn = cutIn(io.in.bits.y)
  val zIn = cutIn(io.in.bits.z)

  //----------------------------------------------------------------------
  // S1 : decompose / classify / exponent comparison / path decision
  val (xsgn, xex, xman0) = FloatChiselUtil.decomposeWithLeading1(xSpec, xIn)
  val (ysgn, yex, yman0) = FloatChiselUtil.decomposeWithLeading1(ySpec, yIn)
  val (zsgn, zex, zman0) = FloatChiselUtil.decomposeWithLeading1(zSpec, zIn)

  dbgPrintf("x   = %d|%d(%d)|%d\n", xsgn, xex, xex.asSInt-xSpec.exBias.S, xman0)
  dbgPrintf("y   = %d|%d(%d)|%d\n", ysgn, yex, yex.asSInt-ySpec.exBias.S, yman0)
  dbgPrintf("z   = %d|%d(%d)|%d\n", zsgn, zex, zex.asSInt-zSpec.exBias.S, zman0)

  val (xzero, xinf, xnan) = FloatChiselUtil.checkValue(xSpec, xIn)
  val (yzero, yinf, ynan) = FloatChiselUtil.checkValue(ySpec, yIn)
  val (zzero, zinf, znan) = FloatChiselUtil.checkValue(zSpec, zIn)

  // if subnormal numbers are disabled, the mantissa of "zero(exponent == 0)"
  // may contain non-zero value and it can affect to the result. to avoid it,
  // we replace the whole mantissa of z with zeroes.
  val xman = Mux(xzero, 0.U, xman0)
  val yman = Mux(yzero, 0.U, yman0)
  val zman = Mux(zzero, 0.U, zman0)

  val xysgn  = xsgn ^ ysgn
  val xyzop  = (xysgn ^ zsgn).asBool // 0->(+), 1->(-)
  val xyzero = xzero || yzero
  val xyinf  = xinf  || yinf
  val xynan  = xnan  || ynan || (xinf && yzero) || (xzero && yinf)
  val xyznan = xynan || znan || (xyinf && zinf && xyzop)
  val xyzinf = (xyinf && !znan && !zinf) || (xyinf && zinf && !xyzop)

  // -------------------------------------------------------------------------
  // exponent check

  val xyMaxEx  = xSpec.exMax + ySpec.exMax
  val xyMinEx  = xSpec.exMin + ySpec.exMin
  val xyExW0   = max(log2Up(xyMaxEx+1), log2Up(abs(xyMinEx) + 1))
  val xyExW    = if (xyMinEx<0) {xyExW0+1} else {xyExW0}
  val xyExBias = xSpec.exBias + ySpec.exBias

  val diffMaxXYMinZ = xyMaxEx - zSpec.exMin
  val diffMaxZMinXY = zSpec.exMax - xyMinEx
  val diffExW       = log2Up(max(diffMaxXYMinZ, diffMaxZMinXY) + 1) + 1

  dbgPrintf(f"xyMaxEx  = $xyMaxEx%d = xMaxEx(${xSpec.exMax}%d) + yMaxEx(${ySpec.exMax}%d)\n")
  dbgPrintf(f"xyMinEx  = $xyMinEx%d = xMaxEx(${xSpec.exMin}%d) + yMaxEx(${ySpec.exMin}%d)\n")
  dbgPrintf(f"xyExW0   = $xyExW0%d\n")
  dbgPrintf(f"xyExW    = $xyExW%d\n")
  dbgPrintf(f"xyExBias = $xyExBias%d\n")

  val xyExNobias = xex.pad(xyExW).asSInt + yex.pad(xyExW).asSInt - xyExBias.S(xyExW.W)
  val zExNobias  = zex.pad(diffExW).asSInt - zSpec.exBias.S(diffExW.W)

  val zEx0  = Mux(zzero,  xyExNobias, zExNobias)
  val xyEx0 = Mux(xyzero, zExNobias, xyExNobias)

  dbgPrintf("xex   = %d\n", xex.pad(xyExW))
  dbgPrintf("yex   = %d\n", yex.pad(xyExW))
  dbgPrintf("zex   = %d\n", zex)
  dbgPrintf("xyEx0 = %d\n", xyEx0)
  dbgPrintf("zexnb = %d\n", zEx0)

  dbgPrintf(f"diffMaxXYMinZ = $diffMaxXYMinZ%d\n")
  dbgPrintf(f"diffMaxZMinXY = $diffMaxZMinXY%d\n")
  dbgPrintf(f"diffExW       = $diffExW%d\n")

  val diffExXYZ = xyEx0.pad(diffExW).asSInt - zEx0.pad(diffExW).asSInt
  val diffExZXY = zEx0.pad(diffExW).asSInt - xyEx0.pad(diffExW).asSInt

  // later we may take abs of diffEx, so we calculate + and - beforehand.
  dbgPrintf("diffExXYZ = %d\n", diffExXYZ)
  dbgPrintf("diffExZXY = %d\n", diffExZXY)

  // zEx - xyEx == +2, +1, 0, -1 and signs are different
  val isNear = (!xyzero) && (!zzero) && (xysgn =/= zsgn) &&
               ((diffExZXY === 2.S) || (diffExZXY ===  1.S) ||
                (diffExZXY === 0.S) || (diffExZXY === -1.S))
  val exXYLargerThanZ = diffExXYZ(diffExW-1) === 0.U // exponent of xy >= z

  // xyEx - zEx >= 0 (if 0 or 1, signs are the same)
  val isFarProd   = (!isNear) && exXYLargerThanZ
  val isFarAddend = (!isNear) && (!isFarProd)    // zEx - xyEx >= 1

  //--------------------------- stage "unpack" ---------------------------
  private val cutUnpack = new SharedPipe("unpack", cutIn.valid, stages.at("unpack"))

  private def cutUnpackDbg[T <: Data](x: T): T = {
    if (enableDebug) {
      cutUnpack(x)
    } else {
      x
    }
  }

  val xmanUnpack        = cutUnpack(xman)
  val ymanUnpack        = cutUnpack(yman)
  val zmanUnpack        = cutUnpack(zman)
  val diffExXYZUnpack   = cutUnpack(diffExXYZ)
  val diffExZXYUnpack   = cutUnpack(diffExZXY)
  val xyEx0Unpack       = cutUnpack(xyEx0)
  val zExNobiasUnpack   = cutUnpack(zExNobias)
  val xyzopUnpack       = cutUnpack(xyzop)
  val xysgnUnpack       = cutUnpack(xysgn)
  val isNearUnpack      = cutUnpack(isNear)
  val isFarProdUnpack   = cutUnpack(isFarProd)
  val xyznanUnpack      = cutUnpack(xyznan)
  val xyzinfUnpack      = cutUnpack(xyzinf)
  val xyzeroUnpack      = cutUnpack(xyzero)
  val zzeroUnpack       = cutUnpack(zzero)
  val zsgnUnpack        = cutUnpack(zsgn)
  val isFarAddendUnpack = cutUnpackDbg(isFarAddend) // dbgPrintf only
  val zexUnpack         = cutUnpackDbg(zex)         // dbgPrintf only

  // -------------------------------------------------------------------------
  // S2 : multiply

  val xyprod = xmanUnpack * ymanUnpack;

  //---------------------------- stage "mul" -----------------------------
  private val cutMul = new SharedPipe("mul", cutUnpack.valid, stages.at("mul"))

  private def cutMulDbg[T <: Data](x: T): T = {
    if (enableDebug) {
      cutMul(x)
    } else {
      x
    }
  }

  val xyprodMul      = cutMul(xyprod)
  val zmanMul        = cutMul(zmanUnpack)
  val diffExXYZMul   = cutMul(diffExXYZUnpack)
  val diffExZXYMul   = cutMul(diffExZXYUnpack)
  val xyEx0Mul       = cutMul(xyEx0Unpack)
  val zExNobiasMul   = cutMul(zExNobiasUnpack)
  val xyzopMul       = cutMul(xyzopUnpack)
  val xysgnMul       = cutMul(xysgnUnpack)
  val isNearMul      = cutMul(isNearUnpack)
  val isFarProdMul   = cutMul(isFarProdUnpack)
  val xyznanMul      = cutMul(xyznanUnpack)
  val xyzinfMul      = cutMul(xyzinfUnpack)
  val xyzeroMul      = cutMul(xyzeroUnpack)
  val zzeroMul       = cutMul(zzeroUnpack)
  val zsgnMul        = cutMul(zsgnUnpack)
  val isFarAddendMul = cutMulDbg(isFarAddendUnpack) // dbgPrintf only
  val zexMul         = cutMulDbg(zexUnpack)         // dbgPrintf only

  // -------------------------------------------------------------------------
  // S3 : far path (product is larger)
  //
  // xyEx > zEx, so |deltaEx| = xyEx - zEx = diffExXYZ >= 2.
  // if diffExXYZ == 1 || 0 and signs are different, use near path.
  // Note that, if diffExXYZ == 1 || 0 and signs are the same, use this path.
  //
  // +--carry
  // |      +---- 1+xSpec.manW + 1+ySpec.manW = 48
  // |      |
  // v.-----+-----.
  // 0**.**********00         |   001.0000000
  //     '---+---'^           |  -  0.00000zzzz
  //         |    |           |     0.11111wwww
  //         |    + round     |     ^
  //         + wSpec.manW+1   |     ex -= 1

  val diffExFarProd = diffExXYZMul.asUInt
  val farProdWidth = max(2 + (1+xSpec.manW) + (1+ySpec.manW)+1,
                         3 + (1+wSpec.manW) + 1)
  val farProdFracWidth = farProdWidth - 4

  when(isFarProdMul) {
    dbgPrintf("diffExFarProd     = %d\n", diffExFarProd)
    dbgPrintf(f"farProdWidth      = $farProdWidth%d\n")
    dbgPrintf(f"farProdFracWidth  = $farProdFracWidth%d\n")
  }

  //      +-- 1+xSpec.manW(=3) + 1+ySpec.manW(=3)
  //  .---+---.|
  // 0**.******|
  // 00z.zzzz  |(no shift)
  //     zzzzz |(>>1)
  //           |zzzzz (>>7 = xyFracWidth + 1)
  //
  val zShiftOutLimit = farProdFracWidth + 1
  val zShiftOutW     = log2Up(zShiftOutLimit)
  val zShiftOut      = (diffExFarProd >= zShiftOutLimit.U)

  when(isFarProdMul) {
    dbgPrintf(f"zShiftOutLimit  = $zShiftOutLimit%d\n")
    dbgPrintf(f"zShiftOutW      = $zShiftOutW%d\n")
    dbgPrintf("zShiftOut       = %d\n", zShiftOut)
  }

  val zManPad       = Cat(0.U(3.W), zmanMul, 0.U((farProdWidth - (1+zSpec.manW) - 3).W))
  val zManInverted  = Mux(xyzopMul, -(zManPad.asSInt), zManPad.asSInt)
  val zManAligned0  = zManInverted.asSInt >> diffExXYZMul(zShiftOutW-1, 0)
  // Fully shifted out: an effective subtraction (xyzop) still owes a borrow
  // (see sumFarProd below); all-ones is that borrow folded into the operand.
  val zManAligned   = Mux(zShiftOut, Fill(farProdWidth, xyzopMul),
    zManAligned0(zManAligned0.getWidth - 1, zManAligned0.getWidth - 1 - (farProdWidth-1)))
  val prodPad       = Cat(0.U(2.W), xyprodMul, 0.U((farProdWidth - 2 - xyprodMul.getWidth).W)) // XXX FIXME

  //---------------- stage "align" (far product path + shared) --------------
  private val cutAlign = new SharedPipe("align", cutMul.valid, stages.at("align"))

  private def cutAlignDbg[T <: Data](x: T): T = {
    if (enableDebug) {
      cutAlign(x)
    } else {
      x
    }
  }

  val prodPadAlign      = cutAlign(prodPad)
  val zManAlignedAlign  = cutAlign(zManAligned)
  val zShiftOutAlign    = cutAlign(zShiftOut)
  // shared by the three paths
  val xyprodAlign       = cutAlign(xyprodMul)
  val xyzopAlign        = cutAlign(xyzopMul)
  val xysgnAlign        = cutAlign(xysgnMul)
  val xyEx0Align        = cutAlign(xyEx0Mul)
  val zExNobiasAlign    = cutAlign(zExNobiasMul)
  val isNearAlign       = cutAlign(isNearMul)
  val isFarProdAlign    = cutAlign(isFarProdMul)
  val xyznanAlign       = cutAlign(xyznanMul)
  val xyzinfAlign       = cutAlign(xyzinfMul)
  val xyzeroAlign       = cutAlign(xyzeroMul)
  val zzeroAlign        = cutAlign(zzeroMul)
  val zsgnAlign         = cutAlign(zsgnMul)
  val isFarAddendAlign  = cutAlignDbg(isFarAddendMul) // dbgPrintf only
  val zexAlign          = cutAlignDbg(zexMul)         // dbgPrintf only
  val zManPadAlign      = cutAlignDbg(zManPad)        // dbgPrintf only
  val zManInvertedAlign = cutAlignDbg(zManInverted)   // dbgPrintf only
  val zManAligned0Align = cutAlignDbg(zManAligned0)   // dbgPrintf only

  // -------------------------------------------------------------------------
  // S4 : far path (product is larger), significand add
  //
  // Effective subtraction (xyzop) with the subtrahend fully shifted out of the
  // window (zShiftOut) loses the borrow: the true value is product - epsilon,
  // just BELOW the windowed product, so a round-to-even tie must round DOWN.
  // zManAligned carries that borrow as all-ones, so this is a plain add.
  val sumFarProd    = prodPadAlign + zManAlignedAlign
  // shiftSticky only reads signals of the segment that ends at "align"; the
  // source computes it after sumFarProd, so its register follows here.
  val shiftSticky   = PriorityEncoder(zManInverted) < diffExXYZMul(zShiftOutW-1, 0) || zShiftOut
  val shiftStickyAlign = cutAlign(shiftSticky)

  when(isFarProdAlign) {
    dbgPrintf("zManPad       = 0b%b(%d), width= %d\n", zManPadAlign     , zManPadAlign     , zManPadAlign     .getWidth.U)
    dbgPrintf("zManInverted  = 0b%b(%d), width= %d\n", zManInvertedAlign, zManInvertedAlign, zManInvertedAlign.getWidth.U)
    dbgPrintf("zManAligned0  = 0b%b(%d), width= %d\n", zManAligned0Align, zManAligned0Align, zManAligned0Align.getWidth.U)
    dbgPrintf("xyProd        = 0b  %b(%d), width=  %d\n", xyprodAlign     , xyprodAlign     , xyprodAlign     .getWidth.U)
    dbgPrintf("zManAligned   = 0b%b(%d), width= %d\n", zManAlignedAlign , zManAlignedAlign , zManAlignedAlign .getWidth.U)
    dbgPrintf("ProdPad       = 0b%b(%d), width= %d\n", prodPadAlign     , prodPadAlign     , prodPadAlign     .getWidth.U)
    dbgPrintf("sumFarProd    = 0b%b(%d), width= %d\n", sumFarProd  , sumFarProd  , sumFarProd  .getWidth.U)
  }

  val farProdSgn    = Mux(xysgnAlign.asBool, ~sumFarProd(sumFarProd.getWidth-1), sumFarProd(sumFarProd.getWidth-1))
  val farProdAbs    = Mux(sumFarProd(sumFarProd.getWidth-1), -sumFarProd, sumFarProd)

  //----------------- stage "add" (far product path + shared) ---------------
  private val cutAdd = new SharedPipe("add", cutAlign.valid, stages.at("add"))

  private def cutAddDbg[T <: Data](x: T): T = {
    if (enableDebug) {
      cutAdd(x)
    } else {
      x
    }
  }

  // shared by the three paths
  val xyEx0Add       = cutAdd(xyEx0Align)
  val zExNobiasAdd   = cutAdd(zExNobiasAlign)
  val xysgnAdd       = cutAdd(xysgnAlign)
  val isNearAdd      = cutAdd(isNearAlign)
  val isFarProdAdd   = cutAdd(isFarProdAlign)
  val xyznanAdd      = cutAdd(xyznanAlign)
  val xyzinfAdd      = cutAdd(xyzinfAlign)
  val xyzeroAdd      = cutAdd(xyzeroAlign)
  val zzeroAdd       = cutAdd(zzeroAlign)
  val zsgnAdd        = cutAdd(zsgnAlign)
  val isFarAddendAdd = cutAddDbg(isFarAddendAlign) // dbgPrintf only
  val zexAdd         = cutAddDbg(zexAlign)         // dbgPrintf only
  val xyprodAdd      = cutAddDbg(xyprodAlign)      // dbgPrintf only

  //--------------------- stage "norm" (shared) -----------------------------
  private val cutNorm = new SharedPipe("norm", cutAdd.valid, stages.at("norm"))

  private def cutNormDbg[T <: Data](x: T): T = {
    if (enableDebug) {
      cutNorm(x)
    } else {
      x
    }
  }

  // shared by the three paths
  val xyEx0Norm        = cutNorm(xyEx0Add)
  val zExNobiasNorm    = cutNorm(zExNobiasAdd)
  val xysgnNorm        = cutNorm(xysgnAdd)
  val isNearNorm       = cutNorm(isNearAdd)
  val isFarProdNorm    = cutNorm(isFarProdAdd)
  val xyznanNorm       = cutNorm(xyznanAdd)
  val xyzinfNorm       = cutNorm(xyzinfAdd)
  val xyzeroNorm       = cutNorm(xyzeroAdd)
  val zzeroNorm        = cutNorm(zzeroAdd)
  val zsgnNorm         = cutNorm(zsgnAdd)
  val isFarAddendNorm  = cutNormDbg(isFarAddendAdd) // dbgPrintf only
  val zexNorm          = cutNormDbg(zexAdd)         // dbgPrintf only

  // -------------------------------------------------------------------------
  // far path (addend is larger)
  //
  // xyEx < zEx, so deltaEx = zEx - xyEx = diffExZXY >= 3.
  // if deltaEx == 2 || 1 || 0 || -1, use near path.
  //
  // 00*.******
  // 000.000xxxxxxxxxxxxxxxx

  val diffExFarAdd = diffExZXYMul.asUInt // >= 3
  val farAddWidth = max(1 + (1+xSpec.manW) + (1+ySpec.manW), 2 + (1+wSpec.manW) + 1)
  val farAddFracWidth = farAddWidth - 2

  when(isFarAddendMul) {
    dbgPrintf("diffExFarAdd   = %d\n", diffExFarAdd)
  }

  val xyShiftOutLimit = farAddFracWidth + 2
  val xyShiftOutW     = log2Up(xyShiftOutLimit)
  val xyShiftOut      = (diffExFarAdd >= xyShiftOutLimit.U)

  val zManFarAddPad      = Cat(0.U(2.W), zmanMul, 0.U((farAddWidth - (1+zSpec.manW) - 2).W))
  val zManFarAddInverted = Mux(xyzopMul, -(zManFarAddPad.asSInt), zManFarAddPad.asSInt)
  val xyManFarAddAligned = Cat(0.U(1.W), (xyprodMul >> diffExFarAdd(xyShiftOutW-1, 0)))

  //--------------------- stage "align" (far addend path) -------------------
  val zManFarAddInvertedAlign = cutAlign(zManFarAddInverted)
  val xyManFarAddAlignedAlign = cutAlign(xyManFarAddAligned)
  val xyShiftOutAlign         = cutAlign(xyShiftOut)
  val zManFarAddPadAlign      = cutAlignDbg(zManFarAddPad) // dbgPrintf only

  val sumFarAdd          = Mux(xyShiftOutAlign, zManFarAddInvertedAlign,
    xyManFarAddAlignedAlign.asSInt + zManFarAddInvertedAlign)
  // as for shiftSticky above, this only reads signals of the segment that ends
  // at "align", so its register follows here.
  val shiftStickyFarAdd  = PriorityEncoder(xyprodMul) < diffExFarAdd(xyShiftOutW-1, 0) || xyShiftOut
  val shiftStickyFarAddAlign = cutAlign(shiftStickyFarAdd)

  when(isFarAddendAlign) {
    dbgPrintf("xyProd              = 0b%b(%d), width= %d\n", xyprodAlign, xyprodAlign, xyprodAlign.getWidth.U)
    dbgPrintf("xyManFarAddAligned  = 0b%b(%d), width= %d\n", xyManFarAddAlignedAlign, xyManFarAddAlignedAlign, xyManFarAddAlignedAlign.getWidth.U)
    dbgPrintf("zManFarAddPad       = 0b%b(%d), width= %d\n", zManFarAddPadAlign     , zManFarAddPadAlign     , zManFarAddPadAlign     .getWidth.U)
    dbgPrintf("zManFarAddInverted  = 0b%b(%d), width= %d\n", zManFarAddInvertedAlign, zManFarAddInvertedAlign, zManFarAddInvertedAlign.getWidth.U)
    dbgPrintf("sumFarAdd           = 0b%b(%d), width= %d\n", sumFarAdd         , sumFarAdd         , sumFarAdd         .getWidth.U)
    dbgPrintf("shiftStickyFarAdd   = 0b%b(%d), width= %d\n", shiftStickyFarAddAlign, shiftStickyFarAddAlign, shiftStickyFarAddAlign.getWidth.U)
  }

  val farAddSgn   = Mux(xysgnAlign.asBool, ~sumFarAdd(sumFarAdd.getWidth-1), sumFarAdd(sumFarAdd.getWidth-1))
  // On this path xyEx <= zEx - 3, so an effective subtraction (xyzop) always
  // leaves sumFarAdd strictly negative: the borrow rides the negate's carry-in
  // (-s - 1 == ~s) and the "not zero" guard cannot fire.
  val farAddDec   = xyzopAlign && shiftStickyFarAddAlign
  val farAddAbs   = Mux(sumFarAdd(sumFarAdd.getWidth-1),
    (~sumFarAdd).asUInt + (!farAddDec).asUInt, sumFarAdd.asUInt)

  // -------------------------------------------------------------------------
  // near path (signs are different)
  //
  // since x*y is a product, the integer part could be larger than 2.
  //
  //  dEx === 2  |  dEx === 1 |  dEx === 0 | dEx === -1  (dEx = zEx - xyEx)
  //   11.xxxxxx |  11.xxxxxx |  01.xxxxxx |  01.xxxxxx
  // -zzz.zzz    | -zz.zzzz   |  -z.zzzzz  |  -0.zzzzzz
  //
  //  2+xmanW+ymanW
  // 0011.xxxxxxxxx
  // 000z.zzz000000
  //  3+zmanW+Pad = 2+xmanW+ymanW
  //          Pad = xmanW+ymanW - zmanW - 2
  //
  // if 2+x.manW+y.manW < w.manW
  //
  // 0011.xxxxxxxxx000
  // 000z.zzz000000000
  //      '--w.manW--'
  //
  // TODO: the case when wSpec.manW > xSpec.manW + ySpec.manW
  //            and when zSpec.manW + 1 > xSpec.manW + ySpec.manW
  //       the intermediate value width might not be correct because FMA
  //       requires intermediate values to have infinite precision.
  //       e.g. if zSpec.manW > xSpec.manW + ySpec.manW, then we need to add
  //            more bits to the least significant part to avoid truncation of
  //            LSBs of z while shifting

  val zNearWidth     = 2+max((1+xSpec.manW)+(1+ySpec.manW), 1+wSpec.manW)
  val zNearFracWidth = xSpec.manW + ySpec.manW
  val zNearPad       = Cat(0.U(3.W), zmanMul, 0.U((zNearWidth - (1+zSpec.manW) - 3).W))
  val zNearShift     = Mux(diffExZXYMul(diffExW-1) === 0.U, (zNearPad << diffExZXYMul(1,0))(zNearWidth-1, 0), (Cat(0.U(1.W), zNearPad) >> 1));

  //------------------------ stage "align" (near path) ----------------------
  val zNearShiftAlign = cutAlign(zNearShift)
  val zNearPadAlign   = cutAlignDbg(zNearPad) // dbgPrintf only

  val nearSum        = Cat(0.U(2.W), xyprodAlign) - zNearShiftAlign
  val nearSgn        = Mux(xysgnAlign.asBool, ~nearSum(zNearWidth-1), nearSum(zNearWidth-1))
  val nearAbs        = Mux(nearSum(zNearWidth-1), -nearSum, nearSum)
  // exact cancellation flag (see nearExactZero below); it is formed next to
  // nearAbs so that only the flag has to cross the remaining boundaries
  val nearZero       = (nearAbs === 0.U)

  // -------------------------------------------------------------------------
  // S5/S6 : one normalisation and one rounder for the three paths
  //
  // The paths are mutually exclusive and their windows have the same shape --
  // the leading one at the MSB, the mantissa below it, then the round bit and
  // the sticky bits -- so the three sums are left justified into one window
  // and share the leading zero count, the normalising shifter and the rounder.
  // The zeros appended below each sum join the sticky, which they do not
  // change, and they do not move the leading one, so the shared window rounds
  // exactly as the three separate ones did. Only the exponent base and the
  // sticky of the bits shifted out at "align" stay per path.
  //
  // xxxx.xxxx   : shift 0, ex += 3
  // 0xxx.xxxxx  : shift 1, ex += 2
  // 00xx.xxxxxx : shift 2, ex += 1
  // 000x.xxxxxxx: shift 3, ex += 0

  val sumWidth = max(farProdWidth, max(farAddWidth + 2, zNearWidth + 1))

  /** `x`, the magnitude of the sum of a path `w` bits wide, in the window. */
  private def justifySum(x: UInt, w: Int): UInt = {
    require(x.getWidth == w, s"expected a ${w} bit sum, got ${x.getWidth}")
    if (w == sumWidth) {
      x
    } else {
      Cat(x, 0.U((sumWidth - w).W))
    }
  }

  val sumAbs = Mux(isNearAlign, justifySum(nearAbs, zNearWidth),
    Mux(isFarProdAlign, justifySum(farProdAbs, farProdWidth),
                        justifySum(farAddAbs, farAddWidth)))
  val sumSgn = Mux(isNearAlign, nearSgn,
    Mux(isFarProdAlign, farProdSgn, farAddSgn))
  val sumShiftSticky = Mux(isNearAlign, false.B,
    Mux(isFarProdAlign, shiftStickyAlign, shiftStickyFarAddAlign))

  //---------------------------- stage "add" --------------------------------
  val sumAbsAdd          = cutAdd(sumAbs)
  val sumSgnAdd          = cutAdd(sumSgn)
  val sumShiftStickyAdd  = cutAdd(sumShiftSticky)
  val nearZeroAdd        = cutAdd(nearZero)
  val nearSumAdd         = cutAddDbg(nearSum)         // dbgPrintf only
  val zNearPadAdd        = cutAddDbg(zNearPadAlign)   // dbgPrintf only
  val zNearShiftAdd      = cutAddDbg(zNearShiftAlign) // dbgPrintf only

  val sumShift = PriorityEncoder(Reverse(sumAbsAdd))
  val sumNorm  = (sumAbsAdd << sumShift)(sumWidth-1, 0)

  when(isNearAdd) {
    dbgPrintf(f"zNearWidth     = $zNearWidth%d\n")
    dbgPrintf(f"zNearFracWidth = $zNearFracWidth%d\n")
    dbgPrintf("zNearPad       = 0b%b(%d), width= %d\n",   zNearPadAdd  , zNearPadAdd  , zNearPadAdd  .getWidth.U)
    dbgPrintf("zNearShift     = 0b%b(%d), width= %d\n",   zNearShiftAdd, zNearShiftAdd, zNearShiftAdd.getWidth.U)
    dbgPrintf("nearSum        = 0b%b(%d), width= %d\n",   nearSumAdd   , nearSumAdd   , nearSumAdd   .getWidth.U)
  }

  //---------------------------- stage "norm" -------------------------------
  val sumNormNorm        = cutNorm(sumNorm)
  val sumShiftNorm       = cutNorm(sumShift)
  val sumSgnNorm         = cutNorm(sumSgnAdd)
  val sumShiftStickyNorm = cutNorm(sumShiftStickyAdd)
  val nearZeroNorm       = cutNorm(nearZeroAdd)

  // the window of the far-addend path is the addend's, the other two are the
  // product's, so only the base of the exponent is selected per path
  val wexBase = Mux(isFarProdNorm || isNearNorm,
    xyEx0Norm + (3.S - sumShiftNorm.asSInt),
    zExNobiasNorm.asSInt - (sumShiftNorm - 2.U).asSInt)

  val wsgn0 = sumSgnNorm
  val wman0 = Wire(UInt(wSpec.manW.W))
  val wex0  = Wire(SInt((xyExW+1).W))
  if (wSpec.manW < min(farProdFracWidth, min(farAddFracWidth, zNearFracWidth))) {
    //                                         + leading 1
    //                            + MSB        |   + mantissa width
    //                          .-+----------. v .-+--------.
    val lsbSum    = sumNormNorm(sumWidth-1-1-wSpec.manW+1)
    val roundSum  = sumNormNorm(sumWidth-1-1-wSpec.manW)
    val stickySum = sumNormNorm(sumWidth-1-1-wSpec.manW-1, 0).orR.asBool || sumShiftStickyNorm

    val wIncSum = FloatChiselUtil.roundIncBySpec(roundSpec, lsbSum, roundSum, stickySum)
    val wmanSum = sumNormNorm(sumWidth-1, sumWidth-1-wSpec.manW)

    // 1.11111111 + Round = 10.00000000 -> ex+1 1.000000
    val carryRoundSum = wIncSum && (wmanSum === maskL(wSpec.manW+1).U)

    wman0 := Mux(carryRoundSum, 0.U, (wmanSum + wIncSum.asUInt)(wSpec.manW-1, 0))
    wex0  := wexBase + wSpec.exBias.S((xyExW+1).W) + Mux(carryRoundSum, 1.S, 0.S)

  } else {
    // w is wider than the exact result, so there is nothing to round
    wman0 := padLSB(wSpec.manW+1, sumNormNorm)
    wex0  := wexBase + wSpec.exBias.S((xyExW+1).W)
  }

  //----------------------------------------------------------------------
  // Special values

  val wZeroAfterAdd = (wex0.asSInt < 0.S)
  val wInfAfterAdd  = wex0 >= maskI(wSpec.exW).S

  // Exact cancellation on the near path: x*y + z == 0 exactly. nearAbs is the
  // magnitude of the near-path significand difference, so nearAbs === 0 (i.e.
  // nearZero) means the exact result is zero. When that happens the normalising
  // shift count is meaningless and the normalized exponent
  // is bogus, so wZeroAfterAdd (which only catches wex0 < 0) does not fire and a
  // spurious tiny normal is emitted. Force +0, matching the ideal
  // RealGeneric.fmadd.
  val nearExactZero = isNearNorm && nearZeroNorm

  //--------------------------- stage "round" ----------------------------
  private val cutRound = new SharedPipe("round", cutNorm.valid, stages.at("round"))

  val wman0Round         = cutRound(wman0)
  val wexBitsRound       = cutRound(wex0(wSpec.exW-1,0))
  val wsgn0Round         = cutRound(wsgn0)
  val wZeroAfterAddRound = cutRound(wZeroAfterAdd)
  val wInfAfterAddRound  = cutRound(wInfAfterAdd)
  val nearExactZeroRound = cutRound(nearExactZero)
  val xyznanRound        = cutRound(xyznanNorm)
  val xyzinfRound        = cutRound(xyzinfNorm)
  val xysgnRound         = cutRound(xysgnNorm)
  val xyzeroRound        = cutRound(xyzeroNorm)
  val zzeroRound         = cutRound(zzeroNorm)
  val zsgnRound          = cutRound(zsgnNorm)

  //----------------------------------------------------------------------
  // Final
  val infOrNaN = wInfAfterAddRound  || xyzinfRound || xyznanRound
  val wzero    = wZeroAfterAddRound || (xyzeroRound && zzeroRound) || nearExactZeroRound

  val wex  = Mux(infOrNaN, maskU(wSpec.exW), Mux(wzero, 0.U(wSpec.exW.W), wexBitsRound))
  val wman = Mux(infOrNaN || wzero, Cat(xyznanRound, 0.U((wSpec.manW-1).W)), wman0Round(wSpec.manW-1, 0))
  // IEEE 754 (2019) 6.3 : if both the product and z are zero, the zero result
  // is negative only when both of the two signs are negative. `truncate`
  // stands for roundTowardNegative here (the same convention as
  // AddFPGeneric), where it is negative unless both are positive. This
  // matches RealGeneric.fmadd.
  val wzeroBothSgn = if(roundSpec == RoundSpec.truncate) {
    xysgnRound | zsgnRound
  } else {
    xysgnRound & zsgnRound
  }
  val wsgn = Mux(xyznanRound, 0.U,
    Mux(xyzinfRound, xysgnRound,
      Mux(xyzeroRound && zzeroRound, wzeroBothSgn,
        Mux(nearExactZeroRound, 0.U, wsgn0Round)))) // exact cancel -> +0
  // +inf + inf = +inf, -inf - inf = inf, but inf - inf = nan.
  // so if xyzinf == true, then xysgn and zsgn should be the same.
  // if infAfterAdd == true, then xyzinf != true.

  dbgPrintf("wsgn  = %b(%d), width = %d\n", wsgn, wsgn, wsgn.getWidth.U)
  dbgPrintf("wex   = %b(%d), width = %d\n", wex,  wex,  wex .getWidth.U)
  dbgPrintf("wman  = %b(%d), width = %d\n", wman, wman, wman.getWidth.U)
  dbgPrintf("w0    = %d|%d|%d\n", wsgn, wex, wman)

  val w0 = Cat(wsgn, wex, wman)

  //---------------------------- stage "out" -----------------------------
  private val cutOut = new SharedPipe("out", cutRound.valid, stages.at("out"))

  io.out.bits  := cutOut(w0)
  io.out.valid := cutOut.valid
  dbgPrintf("x=%x y=%x z=%x, w=%x\n",
    io.in.bits.x, io.in.bits.y, io.in.bits.z, io.out.bits)
}

/** Input bundle and pipeline stage points of [[rial.arith.FusedMulAddFPGeneric]].
 *
 * @see [[rial.util.StagePoint]]
 */
object FusedMulAddFPGeneric {

  /** Input of [[rial.arith.FusedMulAddFPGeneric]]: the two factors `x`, `y`
   *  and the addend `z`.
   *
   * @param xSpec floating point spec of the first input
   * @param ySpec floating point spec of the second input
   * @param zSpec floating point spec of the third input
   */
  class Input(xSpec: RealSpec, ySpec: RealSpec, zSpec: RealSpec) extends Bundle {
    val x = UInt(xSpec.W.W)
    val y = UInt(ySpec.W.W)
    val z = UInt(zSpec.W.W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   *  [[rial.arith.FusedMulAddFPGeneric]]. */
  val stagePoints: Seq[StagePoint] = Seq(
    StagePoint("unpack", 3), StagePoint("mul",   4), StagePoint("align", 3),
    StagePoint("add",    3), StagePoint("norm",  3), StagePoint("round", 3),
    StagePoint("out",    1))
}

/** specialization of FusedMulAddFPGeneric for FP64.
 *
 * @see [[rial.arith.FusedMulAddFPGeneric]]
 */
class FusedMulAddFP64( stage : PipelineStageConfig )
    extends FusedMulAddFPGeneric( RealSpec.Float64Spec,
      RealSpec.Float64Spec,  RealSpec.Float64Spec, RealSpec.Float64Spec,
      RoundSpec.roundToEven, stage) {
}
