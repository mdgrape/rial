package rial.table

import scala.math._

import rial.util._

import chisel3._
import chisel3.util._

/** Evaluates one interval polynomial with Horner's scheme.
 *
 * The coefficients of the interval and, for a polynomial of order greater than
 * zero, the position `x` inside it are given on `io.in`; the result is
 * `io.out`, an unsigned value in `[0.0, 1.0)` scaled by `1 << bitW.outBitW`.
 *
 * The pipeline registers are placed by [[rial.util.PipelineStageConfig]]. The
 * declared boundaries, in datapath order, are
 *
 *  - `"step1"` ... `"stepN"` (`N` = `bitW.nOrder`): one boundary after each
 *    multiply-add of Horner's scheme. The accumulator, `x`, the coefficients
 *    still to be used and the `valid` flag cross these boundaries;
 *  - `"out"`: the module output. The clamping of the result belongs to the
 *    logic just before it.
 *
 * A polynomial of order zero is a plain table lookup and declares `"out"`
 * only. `"in"` is implicit, as always.
 *
 * {{{
 * // let the module place two registers
 * new PolynomialEval(bitW, PipelineStageConfig.default(2))
 * // one register after the first multiply-add
 * new PolynomialEval(bitW, PipelineStageConfig.specify("step1" -> 1))
 * // all the registers at the output, which is what a plain Int means
 * new PolynomialEval(bitW, 2)
 * }}}
 *
 * @constructor create a chisel Module PolynomialEval
 * @param bitW           widths of the input, the coefficients and the
 *                       intermediate results
 * @param stage          pipeline configuration of this module
 * @param clampNegToZero clamp a negative result to zero instead of letting it
 *                       wrap around
 *
 * @see [[rial.util.PipelineStageConfig]]
 * @see [[rial.util.StagePoint]]
 */
class PolynomialEval(
  bitW: PolynomialBitWidth,
  stage: PipelineStageConfig,
  clampNegToZero: Boolean = false,
) extends Module {

  /** `n` registers at the output, as a plain `Int` has always meant here. */
  def this(bitW: PolynomialBitWidth, stage: Int, clampNegToZero: Boolean) = {
    this(bitW, PipelineStageConfig.atOut(stage), clampNegToZero)
  }

  /** `n` registers at the output, as a plain `Int` has always meant here. */
  def this(bitW: PolynomialBitWidth, stage: Int) = {
    this(bitW, PipelineStageConfig.atOut(stage), false)
  }

  val order = bitW.nOrder

  /** Boundaries of the datapath at which pipeline registers can be placed. */
  val stagePoints = PolynomialEval.stagePoints(order)

  val stages = stage.resolve(stagePoints)

  val nStage = stages.total

  def getStage = nStage

  val io = IO(new Bundle {
    val in = Flipped(Valid(new Bundle {
     val coeffs = Input(MixedVec(bitW.coefBits.map{w => SInt(w.W)}))
     val x = if (order != 0) { Some(SInt(bitW.inBitW.W)) } else { None }
    }))
    val out = Valid(UInt(bitW.outBitW.W)) // [0.0, 1.0)
  })

  // ----------------------------------------------------------------------
  // input registers ("in")

  val nIn      = stages.at("in")
  val validIn  = ShiftRegister(io.in.valid, nIn, false.B, true.B)
  val coeffsIn = ShiftRegister(io.in.bits.coeffs, nIn)
  val xInOpt   = if (order != 0) {
    Some(ShiftRegister(io.in.bits.x.get, nIn))
  } else {
    None
  }

  val out      = WireDefault(0.U(bitW.outBitW.W))
  val validOut = WireDefault(false.B)
  io.out := Pipe(validOut, out, stages.at("out"))

  if(order == 0) {
    assert(bitW.coefBits.length == 1)
    // +1 because coefficient is treated as an SInt.
    assert(bitW.outBitW+1 == bitW.coefBits.reduce(_+_),
      f"outBitW=${bitW.outBitW}, coefBits=${bitW.coefBits.mkString("[",",","]")}, "+
      f"sum=${bitW.coefBits.reduce(_+_)}")
    assert(!validIn || coeffsIn(0).asUInt(bitW.outBitW) === 0.U,
      "if nOrder==0, coefficient must be in [0, 1) and the MSB must be zero")

    out := coeffsIn(0).asUInt(bitW.outBitW-1, 0) // if order==0, then there is only one parameter.

    validOut := validIn

  } else {

    // conversion to SInt (-1 ~ +1 range) must be done outside
    val xPrec  = bitW.inPrec

    //           ((c2) * x + c1) * x + c0
    // accVec(2)= ^^^^
    // accVec(1)= ^^^^^^^^^^^^^^
    // accVec(0)= ^^^^^^^^^^^^^^^^^^^^^^^ -> out
    //
    val accVec = Wire(MixedVec((0 to order).map(k => SInt(bitW.calcBits(k).W))))
    accVec(order) := coeffsIn(order)

    // the signals that cross the "step" boundaries: the accumulator, x, the
    // coefficients that are still to be used, and valid
    var accCur   = accVec(order)
    var xCur     = xInOpt.get
    var coefCur  = coeffsIn
    var validCur = validIn

    for(k <- (order-1) to 0 by -1) {
      val xminW   = min(accCur.getWidth, xPrec) // min(calcBits, xPrec)
      val prod    = xCur(xPrec, xPrec - xminW).asSInt * accCur
      val sumFull = (prod >> xminW).asSInt +& coefCur(k)
      accVec(k) := sumFull // drop the MSB (carry)

      assert(!validCur || (k == 0).B || accVec(k) === sumFull, "PolynomialEval: horner "+
        f"intermediate (k=${k}) overflowed calcBits=${bitW.calcBits(k)}")

      val nStep = stages.at(s"step${order - k}")
      accCur   = ShiftRegister(accVec(k), nStep)
      xCur     = ShiftRegister(xCur,      nStep)
      coefCur  = ShiftRegister(coefCur,   nStep)
      validCur = ShiftRegister(validCur,  nStep, false.B, true.B)
    }

    val resMax = Fill(bitW.outBitW, 1.U(1.W))
    if (clampNegToZero) {
      val res = accCur
      out := Mux(res < 0.S, 0.U,
               Mux(res > resMax.zext, resMax, res(bitW.outBitW-1, 0)))
    } else {
      val res = accCur.asUInt
      out := Mux(res > resMax, resMax, res)
    }

    validOut := validCur
  }
}

// ------------------------------------------------------------
// software

object PolynomialEval {

  /** Boundaries of the datapath of a polynomial of order `nOrder`.
   *
   * One boundary after each multiply-add of Horner's scheme, named `"step1"`
   * to `"stepN"` in evaluation order, then the module output. A polynomial of
   * order zero is a plain table lookup and has no `"step"` boundary. The
   * weights are rough, target-independent estimates of the delay of the logic
   * that ends at each boundary: a multiply-add against the final clamp and
   * select.
   */
  def stagePoints(nOrder: Int): Seq[StagePoint] = {
    (1 to nOrder).map(k => StagePoint(s"step$k", 3.0)) :+ StagePoint("out", 1.0)
  }

  private def evalInterval(
    x: BigInt,
    bitW: PolynomialBitWidth,
    tableElem: Seq[Long],
  ): BigInt = {

    var calcW = bitW.calcBits(tableElem.length-1)
    var sum   = BigInt(tableElem (tableElem.length-1))

    // println(f"i=${tableElem.length-1}: inBitW=${bitW.inBitW}, coefW=${bitW.coefBits(tableElem.length-1)}, calcW=${calcW}, x=${x}, sum=${sum}, coef=${tableElem(tableElem.length-1)}")
    for(i <- (tableElem.length-2) to 0 by -1) {
      val xminW = min(calcW, bitW.inBitW - 1)
      val prod  = (x >> (bitW.inBitW - 1 - xminW)) * sum

      //                                 -16                            15
      // [-1, 1] -> signed fixedpoint [10000,..11111, 00000, 00001,..01111]
      //
      // since x is signed and in [-1,1), the position of decimal point is
      // slightly different from its width.
      // if x is unsigned and in [0,1), decimal point is the same as its width.
      // but if x is signed and in [-1,1), decimal point locates width-1.
      sum   = (prod >> xminW) + tableElem(i)
      calcW = bitW.calcBits(i)
      // println(f"i=${i}: inBitW=${bitW.inBitW}, coefW=${bitW.coefBits(i)}, calcW=${calcW}, x=${x}, prod=${prod>>xminW}, sum=${sum}, coef=${tableElem(i)}")
    }
    sum
  }

  def eval(
    addr: Int,
    x: BigInt, // unsigned
    coeffs: PolynomialCoefficients,
  ): BigInt = {
    if(coeffs.bitW.nOrder == 0) {
      // just a table, no calc. no x.
      BigInt(coeffs.coeffs(addr)(0))
    } else {
      // convert unsigend [0, 1) -> signed [-1, 1)
      // 0000,..0111,1000,..1111 -> 1000,..1111,0000,..0111
      // then, (normally, on a hardware implementation) flip the MSB!
      // on a software impl, subtract `1000`.
      val sx = x - (BigInt(1) << coeffs.bitW.inPrec)
      PolynomialEval.evalInterval(sx, coeffs.bitW, coeffs.coeffs(addr))
    }
  }
}
