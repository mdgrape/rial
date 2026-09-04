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
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`, `"out"`: the module input and output. The clamping of the result
 *    belongs to the logic just before `"out"`;
 *  - `"prod1"`, `"step1"` ... `"prodN"`, `"stepN"` (`N` = `bitW.nOrder`): two
 *    boundaries per multiply-add of Horner's scheme, in evaluation order.
 *    `"prodK"` is directly at the output of the multiplier, before the
 *    alignment shift and the addition of the coefficient; `"stepK"` is at the
 *    end of the multiply-add. A polynomial of order zero is a plain table
 *    lookup and has neither.
 *
 * [[getStage]] gives the resulting latency, and `io.out.valid` is
 * `io.in.valid` delayed by that many cycles.
 *
 * {{{
 * // two registers, placed by the module
 * new PolynomialEval(bitW, PipelineStageConfig.default(2))
 * // one register at the multiplier of the first multiply-add, one after it
 * new PolynomialEval(bitW, PipelineStageConfig.specify("prod1" -> 1, "step1" -> 1))
 * // two registers at the module output, which is what a plain Int means
 * new PolynomialEval(bitW, 2)
 * // combinational
 * new PolynomialEval(bitW, PipelineStageConfig.none)
 * }}}
 *
 * @constructor create a chisel Module PolynomialEval
 * @param bitW           widths of the input, the coefficients and the
 *                       intermediate results
 * @param stage          pipeline registers of this module
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

  /** `n` registers at the module output. */
  def this(bitW: PolynomialBitWidth, stage: Int, clampNegToZero: Boolean) = {
    this(bitW, PipelineStageConfig.atOut(stage), clampNegToZero)
  }

  /** `n` registers at the module output. */
  def this(bitW: PolynomialBitWidth, stage: Int) = {
    this(bitW, PipelineStageConfig.atOut(stage), false)
  }

  val order = bitW.nOrder

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
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

  val nIn = stages.at("in")

  private val cutIn = new SharedPipe("in", io.in.valid, nIn)

  val validIn  = cutIn.valid
  val coeffsIn = cutIn(io.in.bits.coeffs)
  val xInOpt   = if (order != 0) {
    Some(cutIn(io.in.bits.x.get))
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
      val idx     = order - k // position of the multiply-add in evaluation order
      val xminW   = min(accCur.getWidth, xPrec) // min(calcBits, xPrec)
      val prod    = xCur(xPrec, xPrec - xminW).asSInt * accCur

      // the raw product, right at the output of the multiplier, so that the
      // register can be absorbed into the DSP block
      val cutProd   = new SharedPipe(s"prod${idx}", validCur, stages.at(s"prod${idx}"))
      val prodProd  = cutProd(prod)
      val xProd     = cutProd(xCur)
      val coefProd  = cutProd(coefCur)
      val validProd = cutProd.valid

      val sumFull = (prodProd >> xminW).asSInt +& coefProd(k)
      accVec(k) := sumFull // drop the MSB (carry)

      assert(!validProd || (k == 0).B || accVec(k) === sumFull, "PolynomialEval: horner "+
        f"intermediate (k=${k}) overflowed calcBits=${bitW.calcBits(k)}")

      val cutStep = new SharedPipe(s"step${idx}", validProd, stages.at(s"step${idx}"))
      accCur   = cutStep(accVec(k))
      xCur     = cutStep(xProd)
      coefCur  = cutStep(coefProd)
      validCur = cutStep.valid
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

  /** Boundaries at which pipeline registers can be placed for a polynomial of
   * order `nOrder`; see [[rial.table.PolynomialEval]].
   */
  def stagePoints(nOrder: Int): Seq[StagePoint] = {
    (1 to nOrder).flatMap(k => {
      Seq(StagePoint(s"prod$k", 2.0), StagePoint(s"step$k", 1.0))
    }) :+ StagePoint("out", 1.0)
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
