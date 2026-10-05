//% @file MathFunctions.scala
//
// Math functions module
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import chisel3._
import chisel3.util._

import rial.table._
import rial.arith.{MulImpl, RealSpec, RoundSpec, SplitMultiplier}
import rial.util._

import scala.collection.mutable.ArrayBuffer
import scala.reflect.ClassTag

/** Evaluates one of a configured set of math functions on a shared datapath.
 *
 * The function to evaluate is selected every cycle by `io.in.bits.sel`, a
 * signal built with [[rial.math.MathFuncConfig.signal]]; `io.in.bits.x` is its
 * argument and `io.in.bits.y`, present only when
 * [[rial.math.FuncKind.ATan2Phase1]] is configured, its second one. The result
 * of the selected function appears on `io.out`.
 *
 * The arguments and the result are floating point numbers of the format the
 * tables of `cfg` are generated for, `cfg.spec`.
 *
 * Pipeline registers can be enabled at the following points of the datapath
 * (names for `PipelineStageConfig.specify`):
 *
 *  - `"in"`        the module input
 *  - `"preMul.*"`  the boundaries of the shared pre-process multiplier: with
 *                  `MulImpl.split` those of [[rial.arith.SplitMultiplier]],
 *                  and `"preMul.out"` alone for a plain product
 *  - `"pre.*"`     the argument reduction of each function
 *  - `"table"`     after the coefficient table lookup
 *  - `"coef"`      after the coefficient selection
 *  - `"poly.*"`    the boundaries of the [[rial.table.PolynomialEval]]
 *                  submodule, see its doc for the names; `"poly"` stands for
 *                  `"poly.out"`
 *  - `"mulArgs.*"` the operands of the post-process multiplier: `"form"`
 *                  after each function has formed the operands it offers and
 *                  `"out"` after the selected ones; `"mulArgs"` stands for
 *                  `"mulArgs.out"`
 *  - `"postMul.*"` the boundaries of the shared post-process multiplier, see
 *                  [[rial.math.PostProcMultiplier]] for the names
 *  - `"post.*"`    the result assembly of each function
 *  - `"out"`       the module output
 *
 * `"pre"` and `"post"` are sets of parallel paths, one per configured
 * function, named with [[rial.math.FuncKind.pointName]]: the boundaries of
 * the argument reduction of `Sqrt` are `"pre.sqrt.*"` and those of its result
 * assembly `"post.sqrt.*"`, with the names below them in the doc of the
 * corresponding submodule. `"pre" -> n` (or `"post" -> n`) adds `n` registers
 * to every path of the set; unless each function keeps its own latency (see
 * below), the paths of a set are padded to the longest of them.
 *
 * `"preMul.*"`, `"mulArgs.*"` and `"postMul.*"` exist only when one of the
 * configured functions needs the corresponding multiplier;
 * [[MathFunctions.stagePoints]] gives the points of a given configuration.
 * The modules that compute the exponent and the special values of a function
 * place their own registers, so they have no point of their own.
 *
 * [[getStage]] gives the resulting latency, and `io.out.valid` is
 * `io.in.valid` with a `sel` that names a configured function, delayed by that
 * many cycles.
 *
 * With `perFunctionLatency` each function keeps the latency of its own path
 * instead of waiting for the longest one, and skips a shared multiplier it
 * does not use: `getStage(fn)` is the latency of `fn` and `getStage` the
 * largest of them. The coefficient table and the output port stay shared, so
 * two operations must not reach them in the same cycle; `io.issueOk`, indexed
 * by the `sel` value itself, says which functions may be issued this cycle
 * and has to be honoured. The results of one function come out in the order
 * it was issued in; the results of different functions come out in the order
 * their latencies give them.
 *
 * With `tagW` an opaque tag rides along: `io.in.bits.tag` enters with the
 * operation and `io.out.bits.tag` comes back with its result.
 *
 * {{{
 * // six registers, placed by the module
 * new MathFunctions(cfg, 6)
 * // hand-picked positions, one of them inside the polynomial
 * new MathFunctions(cfg,
 *   PipelineStageConfig.specify("pre" -> 1, "poly.prod1" -> 1, "post" -> 1))
 * // combinational
 * new MathFunctions(cfg, PipelineStageConfig.none)
 * }}}
 *
 * {{{
 * // one latency per function, with an 8 bit tag
 * val math = Module(new MathFunctions(cfg, 6, perFunctionLatency=true, tagW=8))
 *
 * val sqrtOk = math.io.issueOk.get(cfg.signal(FuncKind.Sqrt))
 * val issue  = wantSqrt && sqrtOk
 *
 * math.io.in.valid        := issue
 * math.io.in.bits.sel     := Mux(issue, cfg.signal(FuncKind.Sqrt), cfg.signalNone())
 * math.io.in.bits.x       := x
 * math.io.in.bits.tag.get := tag
 *
 * when (math.io.out.valid) {
 *   result(math.io.out.bits.tag.get) := math.io.out.bits.z
 * }
 * }}}
 *
 * @param cfg         the functions to support and the floating point format
 *                    of the arguments and the result.
 * @param stage       pipeline registers of this module.
 * @param preMulImpl  structure of the product of the pre-process multiplier.
 * @param postMulImpl structure of the product of the post-process multiplier.
 * @param polyMulImpl structure of the products of the polynomial.
 * @param perFunctionLatency whether each function keeps the latency of its own
 *                    path, see above; the caller then has to honour `io.issueOk`.
 * @param tagW        width of the tag that rides along with an operation; 0 for no tag.
 *
 * @see [[rial.arith.MulImpl]]
 * @see [[rial.util.PipelineStageConfig]]
 */
class MathFunctions(
  val cfg:         MathFuncConfig,
  val stage:       PipelineStageConfig,
  val preMulImpl:  MulImpl = MulImpl.plain,
  val postMulImpl: MulImpl = MulImpl.plain,
  val polyMulImpl: MulImpl = MulImpl.plain,
  val perFunctionLatency: Boolean = false,
  val tagW:               Int     = 0
) extends Module {

  import FuncKind._

  val spec = cfg.spec

  assert(!spec.disableSign, "assuming the float type allows negative values")

  val exW    = spec.exW
  val manW   = spec.manW
  val exBias = spec.exBias

  val bitW   = cfg.bitW
  val coeffs = cfg.coeffs

  private def bitWOf(fn: FuncKind): PolynomialBitWidth = { cfg.bitWOf(fn) }

  /** Boundaries at which pipeline registers can be placed; see the class doc. */
  val stagePoints = MathFunctions.stagePoints(cfg, preMulImpl, postMulImpl,
    polyMulImpl, perFunctionLatency)

  val stages = stage.resolve(stagePoints)

  val nIn  = stages.at("in")
  val nOut = stages.at("out")

  val io = IO(new Bundle {
    val in  = Flipped(Valid(new MathFunctions.Input(cfg, tagW)))
    val out = Valid(new MathFunctions.Output(spec, tagW))
    val issueOk = if (perFunctionLatency) {
      Some(Output(Vec(1 + cfg.funcs.length, Bool())))
    } else {
      None
    }
  })

  /** pipes a `Valid` through one cut of the datapath. */
  private def cutValid[T <: Data](cut: SharedPipe, x: Valid[T]): Valid[T] = {
    val w = Wire(Valid(chiselTypeOf(x.bits)))
    w.valid := cut.valid
    w.bits  := cut(x.bits)
    w
  }

  /** `x` at every age of one cut, `taps(0)` being `x` entering it. */
  private def cutTaps[T <: Data](cut: SharedPipe, x: Valid[T]): Seq[Valid[T]] = {
    cut.taps(x.bits).zip(cut.valids).map { case (b, v) => {
      val w = Wire(Valid(chiselTypeOf(x.bits)))
      w.valid := v
      w.bits  := b
      w
    }}
  }

  // ==========================================================================
  // stage "in": the input, decomposed into its floating point fields
  // ==========================================================================

  private val cutIn = new SharedPipe("in",
    io.in.valid && cfg.isValid(io.in.bits.sel), nIn)

  private val selInTaps = cutIn.taps(io.in.bits.sel)
  private val tagInTaps = io.in.bits.tag.map(t => cutIn.taps(t))

  val selIn = selInTaps.last

  val xdecomp = Module(new DecomposeReal(spec))
  xdecomp.io.real.valid := cutIn.valid
  xdecomp.io.real.bits  := cutIn(io.in.bits.x)

  val xIn: Valid[DecomposedReal] = xdecomp.io.decomp

  // ATan2 also needs the y operand decomposed.
  val yIn: Option[Valid[DecomposedReal]] = if (cfg.has(ATan2Phase1)) {
    val m = Module(new DecomposeReal(spec))
    m.io.real.valid := cutIn.valid
    m.io.real.bits  := cutIn(io.in.bits.y.get)
    Some(m.io.decomp)
  } else {
    None
  }

  // ==========================================================================
  // Shared pre/post-process multipliers (only built if some configured
  // function needs them; sized to the minimum width that satisfies all
  // such functions).
  // ==========================================================================

  val needPreMul = PreProcMultiplier.isNeeded(cfg.funcs)

  val nPreMul = if (needPreMul) { stages.latency("preMul") } else { 0 }

  val preMulLhsW = PreProcMultiplier.lhsW(cfg.funcs, spec)

  private def preMulOperands = new MathFunctions.MulOperands(preMulLhsW, 1 + manW)

  require(!needPreMul || !preMulImpl.isMacc,
    "the pre-process multiplier has no addend to fold into a multiply-add; " +
    "ask for MulImpl.split to slice its product")

  /** the operands of the selected function, driven below once every function
   *  has offered its own set. */
  val preMulIn: Option[Valid[MathFunctions.MulOperands]] = if (needPreMul) {
    Some(Wire(Valid(preMulOperands)))
  } else {
    None
  }

  /** the product of the pre-process multiplier and the valid that leaves it. */
  private val preMulOut: Option[(Bool, UInt)] = preMulIn.map(in => {
    if (preMulImpl.isSplit) {
      val mult = Module(new SplitMultiplier(preMulLhsW, 1 + manW,
        preMulImpl.sliceW.get, stages.sub("preMul")))

      mult.io.in.valid  := in.valid
      mult.io.in.bits.x := in.bits.lhs
      mult.io.in.bits.y := in.bits.rhs
      (mult.io.out.valid, mult.io.out.bits)
    } else {
      // the product as it is written by hand: `MulImpl.plain` asks nothing of
      // the operands, so there is nothing to shape
      val cutMulIn  = new SharedPipe("preMulIn",  in.valid,       stages.at("preMul.in" ))
      val cutMulOut = new SharedPipe("preMulOut", cutMulIn.valid, stages.at("preMul.out"))
      (cutMulOut.valid, cutMulOut(cutMulIn(in.bits.lhs) * cutMulIn(in.bits.rhs)))
    }
  })

  val preMulProd:  Option[UInt] = preMulOut.map({ case (_,  prod) => prod  })
  val preMulValid: Option[Bool] = preMulOut.map({ case (valid, _) => valid })

  // Pre-mult operands, one set per function: each function that wants the
  // multiplier drives its own set in its block below and offers 0 otherwise.
  val preMulArgs: Option[Map[FuncKind, Valid[MathFunctions.MulOperands]]] =
    if (needPreMul) {
      Some(cfg.funcs.map(fn => {
        fn -> WireDefault(0.U.asTypeOf(Valid(preMulOperands)))
      }).toMap)
    } else {
      None
    }

  if (preMulIn.isDefined) {
    preMulIn.get := MuxLookup(selIn, 0.U.asTypeOf(Valid(preMulOperands)))(
      cfg.funcs.map(fn => cfg.signal(fn) -> preMulArgs.get(fn)))
  }

  // ==========================================================================
  // stage "preMul": everything that bypasses the pre-process multiplier
  // ==========================================================================

  private val cutPreMul = new SharedPipe("preMul", cutIn.valid, nPreMul)

  private val selPreTaps = cutPreMul.taps(selIn)
  private val xPreTaps   = cutTaps(cutPreMul, xIn)
  private val yPreTaps   = yIn.map(y => cutTaps(cutPreMul, y))
  private val tagPreTaps = tagInTaps.map(t => cutPreMul.taps(t.last))

  val selPre = selPreTaps.last
  val xPre   = xPreTaps.last
  val yPre   = yPreTaps.map(y => y.last)

  // --------------------------------------------------------------------------

  val needPostMul = PostProcMultiplier.isNeeded(cfg.funcs)

  val nPostMul = if (needPostMul) { stages.latency("postMul") } else { 0 }
  val nMulArgs = if (needPostMul) { stages.latency("mulArgs") } else { 0 }

  /** registers between the operands a function forms and their selection; the
   *  `"in"` of the group sits where the polynomial output already is, so its
   *  registers stand before the operands are formed as well. */
  val nMulArgsForm = if (needPostMul) {
    stages.at("mulArgs.in") + stages.at("mulArgs.form")
  } else {
    0
  }

  /** registers after the selection of the operands. */
  val nMulArgsOut = if (needPostMul) { stages.at("mulArgs.out") } else { 0 }

  /** latency from the polynomial output to the post-process input. */
  val nMulToPost = nMulArgs + nPostMul

  val postMul: Option[PostProcMultiplier] = if (needPostMul) { Some(
    Module(new PostProcMultiplier(
      fracW     = bitW.outBitW,
      manW      = manW,
      roundSpec = RoundSpec.roundToEven,
      stage     = stages.sub("postMul"),
      mulImpl   = postMulImpl,
    ))
  )} else { None }

  private def postMulOperands =
    new MathFunctions.MulOperands(1 + bitW.outBitW, 1 + manW)

  // Post-mult operands, one set per function, driven the same way as the
  // pre-mult ones; the selection itself sits at the "mulArgs.out" boundary
  // below, where the function signal has reached the formed operands.
  val postMulArgs: Option[Map[FuncKind, Valid[MathFunctions.MulOperands]]] =
    if (needPostMul) {
      Some(cfg.funcs.map(fn => {
        fn -> WireDefault(0.U.asTypeOf(Valid(postMulOperands)))
      }).toMap)
    } else {
      None
    }

  // ==========================================================================
  // The latency of each function
  //
  // The table input and the output selection are shared by every function;
  // between them runs one chain. On either side each function has its own
  // path, of its own length when `perFunctionLatency`.
  // ==========================================================================

  val nTable = stages.at("table")
  val nCoef  = stages.at("coef")
  val nPoly  = stages.latency("poly")

  val nTabToPoly = nTable + nCoef + nPoly

  private def skipsPreMul(fn: FuncKind): Boolean = {
    perFunctionLatency && !FuncKind.needPreMult(fn)
  }

  private def skipsPostMul(fn: FuncKind): Boolean = {
    perFunctionLatency && !FuncKind.needPostMult(fn)
  }

  /** registers `fn` crosses between the `"in"` cut and its argument reduction. */
  private def nPreMulOf(fn: FuncKind): Int = {
    if (skipsPreMul(fn)) { 0 } else { nPreMul }
  }

  /** registers `fn` crosses between the polynomial output and its result
   *  assembly. */
  private def nMulToPostOf(fn: FuncKind): Int = {
    if (skipsPostMul(fn)) { 0 } else { nMulToPost }
  }

  /** latency of the argument reduction of `fn`. */
  private def nPreOf(fn: FuncKind): Int = {
    stages.latency("pre." + FuncKind.pointName(fn))
  }

  /** latency of the result assembly of `fn`. */
  private def nPostOf(fn: FuncKind): Int = {
    stages.latency("post." + FuncKind.pointName(fn))
  }

  private def nInToTab(fn: FuncKind): Int = {
    nIn + nPreMulOf(fn) + nPreOf(fn)
  }

  private def nTabToPostOf(fn: FuncKind): Int = {
    nTabToPoly + nMulToPostOf(fn)
  }

  private def nPolyToOut(fn: FuncKind): Int = {
    nMulToPostOf(fn) + nPostOf(fn)
  }

  /** Latency of `fn`: the cycles between `io.in` and its result on `io.out`. */
  def getStage(fn: FuncKind): Int = {
    nInToTab(fn) + nTabToPoly + nPolyToOut(fn) + nOut
  }

  /** The largest latency of the configured functions. */
  def getStage: Int = {
    cfg.funcs.map(fn => getStage(fn)).max
  }

  // ==========================================================================
  // stage "pre": the argument reduction of every function
  // ==========================================================================

  /** the largest latency from `io.in` to the table input. */
  private val nTabMax   = cfg.funcs.map(fn => nInToTab(fn)).max
  private val nPreChain = nTabMax - nIn - nPreMul

  require(nPreChain >= 0)

  private val cutPre = new SharedPipe("pre", cutPreMul.valid, nPreChain)

  private val selTabTaps = cutPre.taps(selPre)
  private val xTabTaps   = cutTaps(cutPre, xPre)
  private val tagTabTaps = tagPreTaps.map(t => cutPre.taps(t.last))

  val xTab = xTabTaps.last

  /** the valid flag at every age between `io.in` and the table input. */
  private val validToTab: Seq[Bool] =
    cutIn.valids ++ cutPreMul.valids.tail ++ cutPre.valids.tail

  /** the function signal at every age between `io.in` and the table input. */
  private val selToTab: Seq[UInt] =
    selInTaps ++ selPreTaps.tail ++ selTabTaps.tail

  /** the tag at every age between `io.in` and the table input. */
  private val tagToTab: Option[Seq[UInt]] = tagInTaps.map(t => {
    t ++ tagPreTaps.get.tail ++ tagTabTaps.get.tail
  })

  /** the argument at every age between the `"in"` cut and the table input:
   *  `xToTab(k)` is `xIn` delayed by `k` registers. */
  private val xToTab: Seq[Valid[DecomposedReal]] = xPreTaps ++ xTabTaps.tail

  /** the argument of `fn` where its argument reduction reads it. */
  private def xPreOf(fn: FuncKind): Valid[DecomposedReal] = {
    xToTab(nPreMulOf(fn))
  }

  /** the second argument of `fn` where its argument reduction reads it. */
  private def yPreOf(fn: FuncKind): Valid[DecomposedReal] = {
    yPreTaps.get(nPreMulOf(fn))
  }

  /** the argument of `fn` at the table input. */
  private def xTabOf(fn: FuncKind): Valid[DecomposedReal] = {
    xToTab(nInToTab(fn) - nIn)
  }

  /** the function signal where the argument reduction of `fn` reads it. */
  private def selPreOf(fn: FuncKind): UInt = {
    selToTab(nIn + nPreMulOf(fn))
  }

  /** the function signal where `fn` reaches the table input. */
  private def selTabOf(fn: FuncKind): UInt = {
    selToTab(nInToTab(fn))
  }

  /** whether an operation of `fn` reaches the table input this cycle.
   *
   * When every function reaches it at the same boundary, the signal there
   * alone says so: it is the one the arriving operation brought. When each
   * function has its own boundary, the signal of `fn` also holds the last
   * operation that passed that one, so its valid has to agree.
   */
  private def arrivesTabOf(fn: FuncKind): Bool = {
    if (perFunctionLatency) {
      xTabOf(fn).valid && selTabOf(fn) === cfg.signal(fn)
    } else {
      selTabOf(fn) === cfg.signal(fn)
    }
  }

  // --------------------------------------------------------------------------
  // The table input, where the argument reductions rejoin: the operation that
  // arrives there this cycle names itself, so the shared chain below carries
  // its signal and its tag.

  /** whether an operation that took `d` cycles to get here reaches the table
   *  input this cycle, one entry per distinct `d`.
   */
  private val arrivesTab: Seq[(Int, Bool)] =
    cfg.funcs.map(fn => nInToTab(fn)).distinct.map(d => {
      val arrivals = cfg.funcs.filter(fn => nInToTab(fn) == d)
      if (arrivals.length == cfg.funcs.length) {
        d -> validToTab(d)
      } else {
        d -> (validToTab(d) &&
          arrivals.map(fn => selToTab(d) === cfg.signal(fn)).reduce(_ || _))
      }
    })

  if (perFunctionLatency) {
    assert(PopCount(arrivesTab.map { case (_, a) => a }) <= 1.U,
      "two operations reached the coefficient table in the same cycle; " +
      "io.issueOk was not honoured")
  }

  private val validTab = arrivesTab.map { case (_, a) => a }.reduce(_ || _)

  val selTab = Mux1H(arrivesTab.map { case (d, a) => a -> selToTab(d) })

  private val tagTab = tagToTab.map(taps => {
    Mux1H(arrivesTab.map { case (d, a) => a -> taps(d) })
  })

  // ==========================================================================
  // stage "table": the coefficient table lookup of every function
  // ==========================================================================

  private val cutTable = new SharedPipe("table", validTab, nTable)

  val selCalc = cutTable(selTab)

  private val tagCalc = tagTab.map(t => cutTable(t))

  // ==========================================================================
  // Polynomial Evaluator (shared), stages "poly.*"
  // ==========================================================================

  val polynomialEval = Module(
    new PolynomialEval(bitW, stages.sub("poly"), mulImpl = polyMulImpl))

  val romLayout = cfg.romLayout

  /** the address of its own table each function offers, driven by its
   *  argument reduction. */
  val tableAdrs = cfg.funcs.map(fn => {
    fn -> WireDefault(0.U.asTypeOf(Valid(UInt(romLayout.adrW(fn).W))))
  }).toMap

  val polynomialDxs = if (bitW.nOrder == 0) { None } else { Some(
    cfg.funcs.map(fn => {
      fn -> WireDefault(0.U.asTypeOf(Valid(UInt(bitW.inBitW.W))))
    }).toMap
  )}

  /** the address `fn` reads the merged ROM at: the address of its own table,
   *  under the constant rows that table starts at. */
  private def romAdrOf(fn: FuncKind): Valid[UInt] = {
    val local = tableAdrs(fn)
    val hiW   = romLayout.addrW - romLayout.adrW(fn)
    val adr   = Wire(Valid(UInt(romLayout.addrW.W)))
    adr.valid := local.valid
    if (hiW == 0) {
      adr.bits := local.bits
    } else {
      adr.bits := Cat((romLayout.offset(fn) >> romLayout.adrW(fn)).U(hiW.W),
        local.bits)
    }
    adr
  }

  val tableRom = Module(new TableCoeffRom(romLayout, nTable))

  tableRom.io.adr := MuxLookup(selTab, 0.U.asTypeOf(tableRom.io.adr))(
    cfg.funcs.map(fn => cfg.signal(fn) -> romAdrOf(fn))
  )
  val polynomialCoef = tableRom.io.cs
  val polynomialDx   = if (bitW.nOrder == 0) { None } else { Some(
    Pipe(MuxLookup(selTab, 0.U.asTypeOf(Valid(UInt(bitW.inBitW.W))))(
      cfg.funcs.map(fn => cfg.signal(fn) -> polynomialDxs.get(fn))
    ), nTable)
  )}

  // ==========================================================================
  // stage "coef": the selected coefficients, on their way to the polynomial
  // ==========================================================================

  private val cutCoef = new SharedPipe("coef", cutTable.valid, nCoef)

  val selPoly = cutCoef(selCalc)

  private val tagPoly = tagCalc.map(t => cutCoef(t))

  private val polyInValid = if (bitW.nOrder == 0) {
    polynomialCoef.valid
  } else {
    polynomialCoef.valid && polynomialDx.get.valid
  }

  private val cutCoefIn = new SharedPipe("coefIn", polyInValid, nCoef)

  polynomialEval.io.in.valid       := cutCoefIn.valid
  polynomialEval.io.in.bits.coeffs := cutCoefIn(polynomialCoef.bits).asTypeOf(
    MixedVec(bitW.coefBits.map(w => SInt(w.W)))
  )
  if(bitW.nOrder != 0) {
    polynomialEval.io.in.bits.x.get := cutCoefIn(polynomialDx.get.bits).asSInt
  }

  // ==========================================================================
  // stage "poly": everything that bypasses the polynomial
  // ==========================================================================

  private val cutPoly = new SharedPipe("poly", cutCoef.valid, nPoly)

  val selPost = cutPoly(selPoly)

  private val tagPost = tagPoly.map(t => cutPoly(t))

  // ==========================================================================
  // stages "mulArgs" and "postMul": everything that bypasses the post-process
  // multiplier
  // ==========================================================================

  private val cutMulArgs = new SharedPipe("mulArgs", cutPoly.valid, nMulArgs)

  private val selMulTaps = cutMulArgs.taps(selPost)
  private val tagMulTaps = tagPost.map(t => cutMulArgs.taps(t))

  val selMul = selMulTaps.last

  /** the function signal where the formed operands are selected. */
  private val selMulArgsForm = selMulTaps(nMulArgsForm)

  // stage "mulArgs.form": the operand set every function offers, each one
  // carried by the valid of its own function; a function that never uses
  // the multiplier offers its zero set as it is
  private val mulArgsForm: Option[Map[FuncKind, Valid[MathFunctions.MulOperands]]] =
    postMulArgs.map(args => cfg.funcs.map(fn => {
      if (FuncKind.needPostMult(fn)) {
        val name = "mulArgsForm" + FuncKind.pointName(fn).capitalize
        fn -> cutValid(new SharedPipe(name, args(fn).valid, nMulArgsForm), args(fn))
      } else {
        fn -> args(fn)
      }
    }).toMap)

  private val postMulIn = mulArgsForm.map(args => {
    MuxLookup(selMulArgsForm, 0.U.asTypeOf(Valid(postMulOperands)))(
      cfg.funcs.map(fn => cfg.signal(fn) -> args(fn)))
  })

  // stage "mulArgs.out": the selected operands
  private val cutMulArgsIn = postMulIn.map(in => {
    new SharedPipe("mulArgsIn", in.valid, nMulArgsOut)
  })

  if (postMul.isDefined) {
    postMul.get.io.in.valid    := cutMulArgsIn.get.valid
    postMul.get.io.in.bits.lhs := cutMulArgsIn.get(postMulIn.get.bits.lhs)
    postMul.get.io.in.bits.rhs := cutMulArgsIn.get(postMulIn.get.bits.rhs)
  }

  private val cutPostMul = new SharedPipe("postMul", cutMulArgs.valid, nPostMul)

  private val selPostMulTaps = cutPostMul.taps(selMul)
  private val tagPostMulTaps = tagMulTaps.map(t => cutPostMul.taps(t.last))

  val selPostMul = selPostMulTaps.last

  // the polynomial result at the boundary the post-process multiplier ends
  // at, for the functions that do not use it but wait for it anyway
  private val cutPostMulRes = new SharedPipe("postMulRes",
    polynomialEval.io.out.valid, nMulToPost)

  val zresPostMul = cutValid(cutPostMulRes, polynomialEval.io.out)

  /** the polynomial result where the result assembly of `fn` reads it. */
  private def zresAtPostOf(fn: FuncKind): Valid[UInt] = {
    if (skipsPostMul(fn)) {
      polynomialEval.io.out
    } else {
      zresPostMul
    }
  }

  // ==========================================================================
  // stage "post": the result assembly of every function
  // ==========================================================================

  /** the largest latency from the polynomial output to the output selection. */
  private val nOutMax = cfg.funcs.map(fn => nPolyToOut(fn)).max

  private val nPostChain = nOutMax - nMulToPost

  require(nPostChain >= 0)

  private val cutPost = new SharedPipe("post", cutPostMul.valid, nPostChain)

  private val selOutTaps = cutPost.taps(selPostMul)
  private val tagOutTaps = tagPostMulTaps.map(t => cutPost.taps(t.last))

  /** the valid flag at every age between the polynomial output and the output
   *  selection. */
  private val validToOut: Seq[Bool] =
    cutMulArgs.valids ++ cutPostMul.valids.tail ++ cutPost.valids.tail

  /** the function signal at every age between the polynomial output and the
   *  output selection. */
  private val selToOut: Seq[UInt] =
    selMulTaps ++ selPostMulTaps.tail ++ selOutTaps.tail

  /** the tag at every age between the polynomial output and the output
   *  selection. */
  private val tagToOut: Option[Seq[UInt]] = tagMulTaps.map(t => {
    t ++ tagPostMulTaps.get.tail ++ tagOutTaps.get.tail
  })

  // ==========================================================================
  // Output selection, stage "out"
  // ==========================================================================

  /** whether a result whose assembly took `e` cycles reaches the output
   *  selection this cycle, one entry per distinct `e`.
   */
  private val arrivesOut: Seq[(Int, Bool)] =
    cfg.funcs.map(fn => nPolyToOut(fn)).distinct.map(e => {
      val arrivals = cfg.funcs.filter(fn => nPolyToOut(fn) == e)
      if (arrivals.length == cfg.funcs.length) {
        e -> validToOut(e)
      } else {
        e -> (validToOut(e) &&
          arrivals.map(fn => selToOut(e) === cfg.signal(fn)).reduce(_ || _))
      }
    })

  if (perFunctionLatency) {
    assert(PopCount(arrivesOut.map { case (_, a) => a }) <= 1.U,
      "two results reached the output selection in the same cycle; " +
      "io.issueOk was not honoured")
  }

  private val validOut = arrivesOut.map { case (_, a) => a }.reduce(_ || _)

  val selOut = Mux1H(arrivesOut.map { case (e, a) => a -> selToOut(e) })

  private val tagOut = tagToOut.map(taps => {
    Mux1H(arrivesOut.map { case (e, a) => a -> taps(e) })
  })

  val zs = cfg.funcs.map(fn => {
    fn -> WireDefault(0.U.asTypeOf(Valid(UInt(spec.W.W))))
  }).toMap

  // `io.out.valid` is `validOut`, which is low in every cycle in which the
  // selected function's `zs` is, so the result is selected on `selOut` alone.
  val zMux = MuxLookup(selOut, 0.U(spec.W.W))(
    cfg.funcs.map(fn => cfg.signal(fn) -> zs(fn).bits))

  private val cutOut = new SharedPipe("out", validOut, nOut)

  io.out.valid  := cutOut.valid
  io.out.bits.z := cutOut(zMux)
  io.out.bits.tag.foreach(t => {
    t := cutOut(tagOut.get)
  })

  // ==========================================================================
  // io.issueOk: the functions that may be issued this cycle
  //
  // Every function passes the table input and the output selection, the two
  // boundaries they share, and reaches each of them a known number of cycles
  // after it was issued. A reservation table on those two latencies says which
  // functions may be issued now.
  // ==========================================================================

  io.issueOk.foreach(ok => {
    val tableInput = cfg.funcs.map(fn => {
      cfg.signal(fn).litValue.toInt -> nInToTab(fn)
    }).toMap

    val output = cfg.funcs.map(fn => {
      cfg.signal(fn).litValue.toInt -> getStage(fn)
    }).toMap

    val table = Module(new ReservationTable(cfg.signalW,
      1 + cfg.funcs.length, Seq(tableInput, output)))

    table.io.in.valid := io.in.valid && cfg.isValid(io.in.bits.sel)
    table.io.in.bits  := io.in.bits.sel

    ok := table.io.issueOk
  })

  // ==========================================================================
  // Per-function pre/tab/nontab/post modules
  //
  // `.collect` so that some funcs does not have the corresponding modules
  // ==========================================================================

  /** the boundaries of the pre-process of `fn`, as the module sees them. */
  private def preStage(fn: FuncKind): PipelineStageConfig = {
    stages.sub("pre." + FuncKind.pointName(fn))
  }

  /** the boundaries of the post-process of `fn`, as the module sees them. */
  private def postStage(fn: FuncKind): PipelineStageConfig = {
    stages.sub("post." + FuncKind.pointName(fn))
  }

  val pres: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtPreProcess       (spec, bitWOf(Sqrt),        preStage(Sqrt       )))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtPreProcess    (spec, bitWOf(InvSqrt),     preStage(InvSqrt    )))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalPreProcess (spec, bitWOf(Reciprocal),  preStage(Reciprocal )))
    case Sin         => Sin         -> Module(new SinCosPreProcess     (spec, bitWOf(Sin),         preStage(Sin        ), cfg.sinCos))
    case SinPi       => SinPi       -> Module(new SinCosPreProcess     (spec, bitWOf(SinPi),       preStage(SinPi      ), cfg.sinCos, scaled = true))
    case SinPiFixed  => SinPiFixed  -> Module(new SinPiFixedPreProcess (spec, bitWOf(SinPiFixed),  preStage(SinPiFixed ), cfg.sinCos))
    case Exp         => Exp         -> Module(new ExpPreProcess        (spec, bitWOf(Exp),         preStage(Exp        )))
    case Pow2        => Pow2        -> Module(new ExpPreProcess        (spec, bitWOf(Pow2),        preStage(Pow2       )))
    case Log         => Log         -> Module(new LogPreProcess        (spec, bitWOf(Log),         preStage(Log        )))
    case Log2Fixed   => Log2Fixed   -> Module(new Log2FixedPreProcess  (spec, bitWOf(Log2Fixed),   preStage(Log2Fixed  )))
    case ACosPhase1  => ACosPhase1  -> Module(new ACosPhase1PreProcess (spec, bitWOf(ACosPhase1),  preStage(ACosPhase1 )))
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2PreProcess (spec, bitWOf(ACosPhase2),  preStage(ACosPhase2 )))
    case ATan2Phase1 => ATan2Phase1 -> Module(new ATan2Phase1PreProcess(spec, bitWOf(ATan2Phase1), preStage(ATan2Phase1)))
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2PreProcess(spec, bitWOf(ATan2Phase2), preStage(ATan2Phase2)))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidPreProcess    (spec, bitWOf(Sigmoid),     preStage(Sigmoid    )))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusPreProcess   (spec, bitWOf(SoftPlus),    preStage(SoftPlus   )))
    case ScaleMixtureGaussian => {
      val (sA, sB) = cfg.scaleMixtureGaussianSigma.get
      ScaleMixtureGaussian -> Module(new ScaleMixtureGaussianPreProcess(
        sA, sB, spec, bitWOf(ScaleMixtureGaussian),
        preStage(ScaleMixtureGaussian)))
    }
    case ATanFixed   => ATanFixed   -> Module(new ATanFixedPreProcess  (spec, bitWOf(ATanFixed),   preStage(ATanFixed  )))
    case ACosFixed   => ACosFixed   -> Module(new ACosFixedPreProcess  (spec, bitWOf(ACosFixed),   preStage(ACosFixed  )))
  }.toMap[FuncKind, Module]

  // the non-table modules run beside the table, the coefficient selection and
  // the polynomial, and place their registers themselves
  private val nonTableStage = PipelineStageConfig.default(nTabToPoly)

  val nontabs: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt       => Sqrt       -> Module(new SqrtNonTable      (spec, bitWOf(Sqrt),        nonTableStage))
    case InvSqrt    => InvSqrt    -> Module(new InvSqrtNonTable   (spec, bitWOf(InvSqrt),     nonTableStage))
    case Reciprocal => Reciprocal -> Module(new ReciprocalNonTable(spec, bitWOf(Reciprocal),  nonTableStage))
    // sincos / acos do not have a dedicated non-table path
    case Exp        => Exp        -> Module(new ExpNonTable       (spec, bitWOf(Exp),         nonTableStage))
    case Pow2       => Pow2       -> Module(new ExpNonTable       (spec, bitWOf(Pow2),        nonTableStage))
    case Log        => Log        -> Module(new LogNonTable       (spec, bitWOf(Log),         nonTableStage))
    // ACosPhase1 reuses SqrtNonTable
    // ACosPhase2 doesn't need one.
  }.toMap

  val posts: Map[FuncKind, Module] = cfg.funcs.collect {
    case Sqrt        => Sqrt        -> Module(new SqrtPostProcess       (spec, bitWOf(Sqrt),        postStage(Sqrt       )))
    case InvSqrt     => InvSqrt     -> Module(new InvSqrtPostProcess    (spec, bitWOf(InvSqrt),     postStage(InvSqrt    )))
    case Reciprocal  => Reciprocal  -> Module(new ReciprocalPostProcess (spec, bitWOf(Reciprocal),  postStage(Reciprocal )))
    case Sin         => Sin         -> Module(new SinCosPostProcess     (spec, bitWOf(Sin),         postStage(Sin        )))
    case SinPi       => SinPi       -> Module(new SinCosPostProcess     (spec, bitWOf(SinPi),       postStage(SinPi      )))
    case SinPiFixed  => SinPiFixed  -> Module(new SinPiFixedPostProcess (spec, bitWOf(SinPiFixed),  postStage(SinPiFixed )))
    case Exp         => Exp         -> Module(new ExpPostProcess        (spec, bitWOf(Exp),         postStage(Exp        )))
    case Pow2        => Pow2        -> Module(new ExpPostProcess        (spec, bitWOf(Pow2),        postStage(Pow2       )))
    case Log         => Log         -> Module(new LogPostProcess        (spec, bitWOf(Log),         postStage(Log        )))
    case Log2Fixed   => Log2Fixed   -> Module(new Log2FixedPostProcess  (spec, bitWOf(Log2Fixed),   postStage(Log2Fixed  )))
    case ACosPhase1  => ACosPhase1  -> Module(new ACosPhase1PostProcess (spec, bitWOf(ACosPhase1),  postStage(ACosPhase1 )))
    case ACosPhase2  => ACosPhase2  -> Module(new ACosPhase2PostProcess (spec, bitWOf(ACosPhase2),  postStage(ACosPhase2 )))
    case ATan2Phase1 => ATan2Phase1 -> Module(new ATan2Phase1PostProcess(spec, bitWOf(ATan2Phase1), postStage(ATan2Phase1)))
    case ATan2Phase2 => ATan2Phase2 -> Module(new ATan2Phase2PostProcess(spec, bitWOf(ATan2Phase2), postStage(ATan2Phase2)))
    case Sigmoid     => Sigmoid     -> Module(new SigmoidPostProcess    (spec, bitWOf(Sigmoid),     postStage(Sigmoid    )))
    case SoftPlus    => SoftPlus    -> Module(new SoftPlusPostProcess   (spec, bitWOf(SoftPlus),    postStage(SoftPlus   )))
    case ScaleMixtureGaussian => ScaleMixtureGaussian ->
      Module(new ScaleMixtureGaussianPostProcess(spec,
        bitWOf(ScaleMixtureGaussian), postStage(ScaleMixtureGaussian)))
    case ATanFixed   => ATanFixed   -> Module(new ATanFixedPostProcess  (spec, bitWOf(ATanFixed),   postStage(ATanFixed  )))
    case ACosFixed   => ACosFixed   -> Module(new ACosFixedPostProcess  (spec, bitWOf(ACosFixed),   postStage(ACosFixed  )))
  }.toMap

  private def toModule[T <: Module](x: Module)(implicit ct: ClassTag[T]): T = {
    require(ct.runtimeClass.isInstance(x),
      s"expected ${ct.runtimeClass.getSimpleName} but got ${x.getClass.getSimpleName}")
    x.asInstanceOf[T]
  }

  // --------------------------------------------------------------------------

  /** [[SqrtNonTable]]'s output at the post-process input.  Sqrt and ACosPhase1
   *  read the same net, so it is delayed once. */
  private val sqrtNtabAtPost = if (cfg.has(Sqrt)) {
    if (cfg.has(ACosPhase1)) {
      require(nMulToPostOf(ACosPhase1) == nMulToPostOf(Sqrt),
        "Sqrt and ACosPhase1 read the same net, so their result assemblies " +
        "have to start at the same boundary")
    }
    Some(Pipe(toModule[SqrtNonTable](nontabs(Sqrt)).io.z, nMulToPostOf(Sqrt)))
  } else {
    None
  }

  if(cfg.has(Sqrt)) {
    val pre  = toModule[SqrtPreProcess ](pres   (Sqrt))
    val ntab = toModule[SqrtNonTable   ](nontabs(Sqrt))
    val post = toModule[SqrtPostProcess](posts  (Sqrt))
    pre.io.x.valid := xPreOf(Sqrt).valid && selPreOf(Sqrt) === cfg.signal(Sqrt)
    pre.io.x.bits  := xPreOf(Sqrt).bits
    tableAdrs(Sqrt) := pre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(Sqrt) := pre.io.dx.get
    }

    ntab.io.x.valid := xTabOf(Sqrt).valid && selTabOf(Sqrt) === cfg.signal(Sqrt)
    ntab.io.x.bits  := xTabOf(Sqrt).bits

    val ntabAtPost = sqrtNtabAtPost.get

    post.io.in.valid     := ntabAtPost.valid && zresAtPostOf(Sqrt).valid
    post.io.in.bits.info := ntabAtPost.bits
    post.io.in.bits.zres := zresAtPostOf(Sqrt).bits

    zs(Sqrt) := post.io.z
  }

  // --------------------------------------------------------------------------

  if(cfg.has(InvSqrt)) {
    val pre  = toModule[InvSqrtPreProcess ](pres   (InvSqrt))
    val ntab = toModule[InvSqrtNonTable   ](nontabs(InvSqrt))
    val post = toModule[InvSqrtPostProcess](posts  (InvSqrt))
    pre.io.x.valid := xPreOf(InvSqrt).valid && selPreOf(InvSqrt) === cfg.signal(InvSqrt)
    pre.io.x.bits  := xPreOf(InvSqrt).bits
    tableAdrs(InvSqrt) := pre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(InvSqrt) := pre.io.dx.get
    }

    ntab.io.x.valid := xTabOf(InvSqrt).valid && selTabOf(InvSqrt) === cfg.signal(InvSqrt)
    ntab.io.x.bits  := xTabOf(InvSqrt).bits

    val ntabAtPost = Pipe(ntab.io.z, nMulToPostOf(InvSqrt))

    post.io.in.valid     := ntabAtPost.valid && zresAtPostOf(InvSqrt).valid
    post.io.in.bits.info := ntabAtPost.bits
    post.io.in.bits.zres := zresAtPostOf(InvSqrt).bits

    zs(InvSqrt) := post.io.z
  }

  // --------------------------------------------------------------------------

  if(cfg.has(Reciprocal)) {
    val pre  = toModule[ReciprocalPreProcess ](pres   (Reciprocal))
    val ntab = toModule[ReciprocalNonTable   ](nontabs(Reciprocal))
    val post = toModule[ReciprocalPostProcess](posts  (Reciprocal))
    pre.io.x.valid := xPreOf(Reciprocal).valid && selPreOf(Reciprocal) === cfg.signal(Reciprocal)
    pre.io.x.bits  := xPreOf(Reciprocal).bits
    tableAdrs(Reciprocal) := pre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(Reciprocal) := pre.io.dx.get
    }

    ntab.io.x.valid := xTabOf(Reciprocal).valid && selTabOf(Reciprocal) === cfg.signal(Reciprocal)
    ntab.io.x.bits  := xTabOf(Reciprocal).bits

    val ntabAtPost = Pipe(ntab.io.z, nMulToPostOf(Reciprocal))

    post.io.in.valid     := ntabAtPost.valid && zresAtPostOf(Reciprocal).valid
    post.io.in.bits.info := ntabAtPost.bits
    post.io.in.bits.zres := zresAtPostOf(Reciprocal).bits

    zs(Reciprocal) := post.io.z
  }

  // --------------------------------------------------------------------------
  // Sin / Cos

  if(cfg.has(Sin) || cfg.has(Cos)) {
    assert(preMulProd.isDefined, "Sin/Cos requires the pre-process multiplier")
    assert(postMul.isDefined, "Sin/Cos requires the post-process multiplier")

    val sincosPre  = toModule[SinCosPreProcess ](pres (Sin))
    val sincosPost = toModule[SinCosPostProcess](posts(Sin))

    // ----- pre-mult operands: 1/π * (1 + x.man), the same set for both -----
    Seq(Sin, Cos).foreach(fn => {
      preMulArgs.get(fn).valid    := xIn.valid
      preMulArgs.get(fn).bits.lhs := SinCosPreMulArgs.lhs(spec, preMulLhsW)
      preMulArgs.get(fn).bits.rhs := Cat(1.U(1.W), xIn.bits.man)
    })

    // ----- pre stage: x and sel on the far side of the pre-mult -----
    val isSinAtPre     = selPreOf(Sin) === cfg.signal(Sin)
    val sincosSelAtPre = isSinAtPre || (selPreOf(Cos) === cfg.signal(Cos))

    sincosPre.io.isSin := isSinAtPre
    sincosPre.io.x.valid := xPreOf(Sin).valid && sincosSelAtPre
    sincosPre.io.x.bits  := xPreOf(Sin).bits

    sincosPre.io.prod.valid := preMulValid.get && sincosSelAtPre
    sincosPre.io.prod.bits  := SinCosPreMulArgs.prod(spec, preMulProd.get)

    // ----- table -----
    tableAdrs(Sin) := sincosPre.io.adr
    tableAdrs(Cos) := sincosPre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(Sin) := sincosPre.io.dx.get
      polynomialDxs.get(Cos) := sincosPre.io.dx.get
    }

    // ----- post-mult operands: polynomial * (1 + yman) -----
    val sincosInfoAtPost = Pipe(sincosPre.io.info, nTabToPoly)

    Seq(Sin, Cos).foreach(fn => {
      postMulArgs.get(fn).valid    :=
        polynomialEval.io.out.valid && sincosInfoAtPost.valid
      postMulArgs.get(fn).bits.lhs := Cat(1.U(1.W), polynomialEval.io.out.bits)
      postMulArgs.get(fn).bits.rhs := Cat(1.U(1.W), sincosInfoAtPost.bits.yman)
    })

    // ----- post -----
    sincosPost.io.pre := Pipe(sincosInfoAtPost, nMulToPostOf(Sin))
    sincosPost.io.mul := postMul.get.io.out

    zs(Sin) := sincosPost.io.z
    zs(Cos) := sincosPost.io.z
  }

  // --------------------------------------------------------------------------
  // SinPi / CosPi -- the caller scaled the argument, so the constant of the
  // pre-multiplication is one and the pre-process multiplier stays out of it.

  if(cfg.has(SinPi) || cfg.has(CosPi)) {
    assert(postMul.isDefined, "SinPi/CosPi requires the post-process multiplier")

    val sinpiPre  = toModule[SinCosPreProcess ](pres (SinPi))
    val sinpiPost = toModule[SinCosPostProcess](posts(SinPi))

    // ----- pre stage: the argument and its product with the constant one -----
    val isSinAtPre    = selPreOf(SinPi) === cfg.signal(SinPi)
    val sinpiSelAtPre = isSinAtPre || (selPreOf(CosPi) === cfg.signal(CosPi))

    sinpiPre.io.isSin := isSinAtPre
    sinpiPre.io.x.valid := xPreOf(SinPi).valid && sinpiSelAtPre
    sinpiPre.io.x.bits  := xPreOf(SinPi).bits

    sinpiPre.io.prod.valid := xPreOf(SinPi).valid && sinpiSelAtPre
    sinpiPre.io.prod.bits  :=
      SinCosPreMulArgs.prodOfOne(spec, xPreOf(SinPi).bits.man)

    // ----- table -----
    tableAdrs(SinPi) := sinpiPre.io.adr
    tableAdrs(CosPi) := sinpiPre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(SinPi) := sinpiPre.io.dx.get
      polynomialDxs.get(CosPi) := sinpiPre.io.dx.get
    }

    // ----- post-mult operands: polynomial * (1 + yman) -----
    val sinpiInfoAtPost = Pipe(sinpiPre.io.info, nTabToPoly)

    Seq(SinPi, CosPi).foreach(fn => {
      postMulArgs.get(fn).valid    :=
        polynomialEval.io.out.valid && sinpiInfoAtPost.valid
      postMulArgs.get(fn).bits.lhs := Cat(1.U(1.W), polynomialEval.io.out.bits)
      postMulArgs.get(fn).bits.rhs := Cat(1.U(1.W), sinpiInfoAtPost.bits.yman)
    })

    // ----- post -----
    sinpiPost.io.pre := Pipe(sinpiInfoAtPost, nMulToPostOf(SinPi))
    sinpiPost.io.mul := postMul.get.io.out

    zs(SinPi) := sinpiPost.io.z
    zs(CosPi) := sinpiPost.io.z
  }

  // --------------------------------------------------------------------------
  // SinPiFixed / CosPiFixed -- one region of the table per binade of the
  // reduced argument, each holding sin(π y) itself, so the result assembly
  // normalizes it and neither shared multiplier takes part.

  if(cfg.has(SinPiFixed) || cfg.has(CosPiFixed)) {
    val pre  = toModule[SinPiFixedPreProcess ](pres (SinPiFixed))
    val post = toModule[SinPiFixedPostProcess](posts(SinPiFixed))

    // ----- pre stage: the argument and its product with the constant one -----
    val isSinAtPre = selPreOf(SinPiFixed) === cfg.signal(SinPiFixed)
    val selAtPre   = isSinAtPre || (selPreOf(CosPiFixed) === cfg.signal(CosPiFixed))

    pre.io.isSin := isSinAtPre
    pre.io.x.valid := xPreOf(SinPiFixed).valid && selAtPre
    pre.io.x.bits  := xPreOf(SinPiFixed).bits

    pre.io.prod.valid := xPreOf(SinPiFixed).valid && selAtPre
    pre.io.prod.bits  :=
      SinCosPreMulArgs.prodOfOne(spec, xPreOf(SinPiFixed).bits.man)

    // ----- table -----
    tableAdrs(SinPiFixed) := pre.io.adr
    tableAdrs(CosPiFixed) := pre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(SinPiFixed) := pre.io.dx.get
      polynomialDxs.get(CosPiFixed) := pre.io.dx.get
    }

    // ----- post -----
    post.io.pre  := Pipe(pre.io.info, nTabToPostOf(SinPiFixed))
    post.io.zres := zresAtPostOf(SinPiFixed)

    zs(SinPiFixed) := post.io.z
    zs(CosPiFixed) := post.io.z
  }

  // --------------------------------------------------------------------------
  // Exp

  if(cfg.has(Exp)) {
    assert(preMulProd.isDefined, "Exp requires the pre-process multiplier")

    val pre  = toModule[ExpPreProcess ](pres   (Exp))
    val ntab = toModule[ExpNonTable   ](nontabs(Exp))
    val post = toModule[ExpPostProcess](posts  (Exp))

    // ----- pre-mult operands: log2(e) * (1 + x.man) -----
    preMulArgs.get(Exp).valid    := xIn.valid
    preMulArgs.get(Exp).bits.lhs := ExpPreMulArgs.lhs(spec, preMulLhsW)
    preMulArgs.get(Exp).bits.rhs := Cat(1.U(1.W), xIn.bits.man)

    // ----- pre stage -----
    val expSelAtPre = selPreOf(Exp) === cfg.signal(Exp)

    pre.io.x.valid := xPreOf(Exp).valid && expSelAtPre
    pre.io.x.bits  := xPreOf(Exp).bits
    pre.io.prod.valid := preMulValid.get && expSelAtPre
    pre.io.prod.bits  := ExpPreMulArgs.prod(spec, preMulProd.get)

    // ----- table -----
    tableAdrs(Exp) := pre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(Exp) := pre.io.dx.get
    }

    // ----- non-table ------
    ntab.io.x.valid := xTabOf(Exp).valid && selTabOf(Exp) === cfg.signal(Exp)
    ntab.io.x.bits  := xTabOf(Exp).bits
    ntab.io.info := pre.io.info

    // ----- post -----
    val ntabAtPost = Pipe(ntab.io.z, nMulToPostOf(Exp))

    post.io.in.valid     := ntabAtPost.valid && zresAtPostOf(Exp).valid
    post.io.in.bits.info := ntabAtPost.bits
    post.io.in.bits.zres := zresAtPostOf(Exp).bits

    zs(Exp) := post.io.z
  }

  // --------------------------------------------------------------------------
  // Pow2 -- the caller scaled the argument, so the constant of the
  // pre-multiplication is one and the pre-process multiplier stays out of it.

  if(cfg.has(Pow2)) {
    val pre  = toModule[ExpPreProcess ](pres   (Pow2))
    val ntab = toModule[ExpNonTable   ](nontabs(Pow2))
    val post = toModule[ExpPostProcess](posts  (Pow2))

    // ----- pre stage: the argument and its product with the constant one -----
    val pow2SelAtPre = selPreOf(Pow2) === cfg.signal(Pow2)

    pre.io.x.valid := xPreOf(Pow2).valid && pow2SelAtPre
    pre.io.x.bits  := xPreOf(Pow2).bits
    pre.io.prod.valid := xPreOf(Pow2).valid && pow2SelAtPre
    pre.io.prod.bits  := ExpPreMulArgs.prodOfOne(spec, xPreOf(Pow2).bits.man)

    // ----- table -----
    tableAdrs(Pow2) := pre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(Pow2) := pre.io.dx.get
    }

    // ----- non-table ------
    ntab.io.x.valid := xTabOf(Pow2).valid && selTabOf(Pow2) === cfg.signal(Pow2)
    ntab.io.x.bits  := xTabOf(Pow2).bits
    ntab.io.info := pre.io.info

    // ----- post -----
    val ntabAtPost = Pipe(ntab.io.z, nMulToPostOf(Pow2))

    post.io.in.valid     := ntabAtPost.valid && zresAtPostOf(Pow2).valid
    post.io.in.bits.info := ntabAtPost.bits
    post.io.in.bits.zres := zresAtPostOf(Pow2).bits

    zs(Pow2) := post.io.z
  }

  // --------------------------------------------------------------------------
  // Log -- uses the post-process multiplier (poly * const), no pre-mult.

  if(cfg.has(Log)) {
    assert(postMul.isDefined, "Log requires the post-process multiplier")

    val pre  = toModule[LogPreProcess ](pres   (Log))
    val ntab = toModule[LogNonTable   ](nontabs(Log))
    val post = toModule[LogPostProcess](posts  (Log))

    pre.io.x.valid := xPreOf(Log).valid && selPreOf(Log) === cfg.signal(Log)
    pre.io.x.bits  := xPreOf(Log).bits
    tableAdrs(Log) := pre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(Log) := pre.io.dx.get
    }

    ntab.io.x.valid := xTabOf(Log).valid && selTabOf(Log) === cfg.signal(Log)
    ntab.io.x.bits  := xTabOf(Log).bits

    // ----- post-mult arg builder (combinational, at the polynomial output) -----
    val mulArgs = Module(new LogPostMultArgs(spec, bitW))
    mulArgs.io.info := ntab.io.z
    mulArgs.io.zres := polynomialEval.io.out.bits

    postMulArgs.get(Log).valid    :=
      mulArgs.io.mul.valid && polynomialEval.io.out.valid
    postMulArgs.get(Log).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(Log).bits.rhs := mulArgs.io.mul.bits.rhs

    // ----- post -----
    post.io.info := Pipe(ntab.io.z, nMulToPostOf(Log))
    post.io.mul  := postMul.get.io.out

    zs(Log) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ACosPhase1 -- reads the table of Sqrt, and reuses its NonTable.

  if(cfg.has(ACosPhase1)) {
    assert(cfg.has(Sqrt),
      "ACosPhase1 requires Sqrt (cfg.normalize should pull it in automatically)")

    val pre  = toModule[ACosPhase1PreProcess ](pres   (ACosPhase1))
    val post = toModule[ACosPhase1PostProcess](posts  (ACosPhase1))
    val sqrtNtab = toModule[SqrtNonTable  ](nontabs(Sqrt))

    pre.io.x.valid := xPreOf(ACosPhase1).valid && selPreOf(ACosPhase1) === cfg.signal(ACosPhase1)
    pre.io.x.bits  := xPreOf(ACosPhase1).bits

    tableAdrs(ACosPhase1) := pre.io.adr

    when (arrivesTabOf(ACosPhase1)) {
      sqrtNtab.io.x := pre.io.y
    }

    if (bitW.nOrder != 0) {
      polynomialDxs.get(ACosPhase1) := pre.io.dx.get
    }

    val flagsAtPost = Pipe(pre.io.flags, nTabToPostOf(ACosPhase1))
    val ntabAtPost  = sqrtNtabAtPost.get

    post.io.in.valid       := ntabAtPost.valid && zresAtPostOf(ACosPhase1).valid && flagsAtPost.valid
    post.io.in.bits.info   := ntabAtPost.bits
    post.io.in.bits.zres   := zresAtPostOf(ACosPhase1).bits
    post.io.in.bits.flags  := flagsAtPost.bits

    zs(ACosPhase1) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ACosPhase2 -- own table, uses the post-process multiplier (poly * x).

  if(cfg.has(ACosPhase2)) {
    assert(postMul.isDefined, "ACosPhase2 requires the post-process multiplier")

    val pre  = toModule[ACosPhase2PreProcess ](pres (ACosPhase2))
    val post = toModule[ACosPhase2PostProcess](posts(ACosPhase2))

    pre.io.x.valid := xPreOf(ACosPhase2).valid && selPreOf(ACosPhase2) === cfg.signal(ACosPhase2)
    pre.io.x.bits  := xPreOf(ACosPhase2).bits
    tableAdrs(ACosPhase2) := pre.io.adr

    if (bitW.nOrder != 0) {
      polynomialDxs.get(ACosPhase2) := pre.io.dx.get
    }

    // post-mult args at the polynomial output: lhs = polynomial w/ hidden 1,
    // rhs = x.manW1.
    val xManW1AtPost = Pipe(pre.io.xManW1, nTabToPoly)
    val mulArgs = Module(new ACosPhase2PostMultArgs(spec, bitW))
    mulArgs.io.zres   := polynomialEval.io.out.bits
    mulArgs.io.xManW1 := xManW1AtPost

    postMulArgs.get(ACosPhase2).valid    :=
      mulArgs.io.mul.valid && polynomialEval.io.out.valid
    postMulArgs.get(ACosPhase2).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(ACosPhase2).bits.rhs := mulArgs.io.mul.bits.rhs

    val flagsAtAfterMul = Pipe(pre.io.flags, nTabToPostOf(ACosPhase2))
    val xexAtAfterMul   = Pipe(pre.io.xex,   nTabToPostOf(ACosPhase2))

    post.io.flags := flagsAtAfterMul
    post.io.xex   := xexAtAfterMul
    post.io.mul   := postMul.get.io.out

    zs(ACosPhase2) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ATan2Phase1 -- reads the table of Reciprocal.
  //
  // Uses the post-process multiplier to assemble (1+min) * (2/(1+max)).

  if(cfg.has(ATan2Phase1)) {
    assert(cfg.has(Reciprocal),
      "ATan2Phase1 requires Reciprocal (cfg.normalize should pull it in automatically)")
    assert(postMul.isDefined, "ATan2Phase1 requires the post-process multiplier")
    assert(yIn.isDefined,     "ATan2Phase1 requires y as an additional input")

    val pre    = toModule[ATan2Phase1PreProcess ](pres (ATan2Phase1))
    val post   = toModule[ATan2Phase1PostProcess](posts(ATan2Phase1))

    pre.io.x.valid := xPreOf(ATan2Phase1).valid && selPreOf(ATan2Phase1) === cfg.signal(ATan2Phase1)
    pre.io.x.bits  := xPreOf(ATan2Phase1).bits
    pre.io.y.valid := yPreOf(ATan2Phase1).valid && selPreOf(ATan2Phase1) === cfg.signal(ATan2Phase1)
    pre.io.y.bits  := yPreOf(ATan2Phase1).bits

    tableAdrs(ATan2Phase1) := pre.io.adr

    if (bitW.nOrder != 0) {
      polynomialDxs.get(ATan2Phase1) := pre.io.dx.get
    }

    // ---- post-mult args at the polynomial output: lhs (poly+hidden 1), rhs (1+minMan) ----
    val infoAtPost = Pipe(pre.io.info, nTabToPoly)
    val mulArgs = Module(new ATan2Phase1PostMultArgs(spec, bitW))
    mulArgs.io.zres      := polynomialEval.io.out.bits
    mulArgs.io.maxXYMan0 := infoAtPost.bits.maxXYMan0
    mulArgs.io.minMan    := infoAtPost.bits.minMan
    mulArgs.io.valid     := infoAtPost.valid && polynomialEval.io.out.valid

    postMulArgs.get(ATan2Phase1).valid    := mulArgs.io.mul.valid
    postMulArgs.get(ATan2Phase1).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(ATan2Phase1).bits.rhs := mulArgs.io.mul.bits.rhs

    // the post-process reads only zex and xySameMan
    val postInfo = Wire(Valid(new ATan2Phase1PostInfo(spec)))
    postInfo.valid          := infoAtPost.valid
    postInfo.bits.zex       := infoAtPost.bits.zex
    postInfo.bits.xySameMan := infoAtPost.bits.xySameMan

    val infoAtAfterMul  = Pipe(postInfo,     nMulToPostOf(ATan2Phase1))
    val flagsAtAfterMul = Pipe(pre.io.flags, nTabToPostOf(ATan2Phase1))

    post.io.info  := infoAtAfterMul
    post.io.flags := flagsAtAfterMul
    post.io.mul   := postMul.get.io.out

    zs(ATan2Phase1) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ATan2Phase2 -- own table; uses the post-process multiplier (atan(t)/t * t).

  if(cfg.has(ATan2Phase2)) {
    assert(postMul.isDefined, "ATan2Phase2 requires the post-process multiplier")

    val pre  = toModule[ATan2Phase2PreProcess ](pres (ATan2Phase2))
    val post = toModule[ATan2Phase2PostProcess](posts(ATan2Phase2))

    pre.io.x.valid := xPreOf(ATan2Phase2).valid && selPreOf(ATan2Phase2) === cfg.signal(ATan2Phase2)
    pre.io.x.bits  := xPreOf(ATan2Phase2).bits
    tableAdrs(ATan2Phase2) := pre.io.adr

    if (bitW.nOrder != 0) {
      polynomialDxs.get(ATan2Phase2) := pre.io.dx.get
    }

    // ---- post-mult args at the polynomial output ----
    val xManW1AtPost   = Pipe(pre.io.xManW1,   nTabToPoly)
    val isLinearAtPost = Pipe(pre.io.isLinear, nTabToPoly)
    val mulArgs = Module(new ATan2Phase2PostMultArgs(spec, bitW))
    mulArgs.io.zres     := polynomialEval.io.out.bits
    mulArgs.io.isLinear := isLinearAtPost.bits
    mulArgs.io.xManW1   := xManW1AtPost

    postMulArgs.get(ATan2Phase2).valid    :=
      mulArgs.io.mul.valid && polynomialEval.io.out.valid
    postMulArgs.get(ATan2Phase2).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(ATan2Phase2).bits.rhs := mulArgs.io.mul.bits.rhs

    val flagsAtAfterMul    = Pipe(pre.io.flags,    nTabToPostOf(ATan2Phase2))
    val xexAtAfterMul      = Pipe(pre.io.xex,      nTabToPostOf(ATan2Phase2))
    val isLinearAtAfterMul = Pipe(isLinearAtPost,  nMulToPostOf(ATan2Phase2))

    post.io.flags    := flagsAtAfterMul
    post.io.xex      := xexAtAfterMul
    post.io.isLinear := isLinearAtAfterMul
    post.io.mul      := postMul.get.io.out

    zs(ATan2Phase2) := post.io.z
  }

  // --------------------------------------------------------------------------
  // Sigmoid -- own table, no pre/post-mult (the post-process fold is a shift).

  if(cfg.has(Sigmoid)) {
    val pre  = toModule[SigmoidPreProcess ](pres (Sigmoid))
    val post = toModule[SigmoidPostProcess](posts(Sigmoid))

    pre.io.x.valid := xPreOf(Sigmoid).valid && selPreOf(Sigmoid) === cfg.signal(Sigmoid)
    pre.io.x.bits  := xPreOf(Sigmoid).bits
    tableAdrs(Sigmoid) := pre.io.adr

    if (bitW.nOrder != 0) {
      polynomialDxs.get(Sigmoid) := pre.io.dx.get
    }

    val infoAtPost = Pipe(pre.io.info, nTabToPostOf(Sigmoid))
    post.io.in.valid     := infoAtPost.valid && zresAtPostOf(Sigmoid).valid
    post.io.in.bits.info := infoAtPost.bits
    post.io.in.bits.zres := zresAtPostOf(Sigmoid).bits

    zs(Sigmoid) := post.io.z
  }

  // --------------------------------------------------------------------------
  // SoftPlus -- own table (concat pos+neg, Mux'd by xsgn), no pre/post-mult.

  if(cfg.has(SoftPlus)) {
    val pre  = toModule[SoftPlusPreProcess ](pres (SoftPlus))
    val post = toModule[SoftPlusPostProcess](posts(SoftPlus))

    pre.io.x.valid := xPreOf(SoftPlus).valid && selPreOf(SoftPlus) === cfg.signal(SoftPlus)
    pre.io.x.bits  := xPreOf(SoftPlus).bits
    tableAdrs(SoftPlus) := pre.io.adr

    if (bitW.nOrder != 0) {
      polynomialDxs.get(SoftPlus) := pre.io.dx.get
    }

    val infoAtPost = Pipe(pre.io.info, nTabToPostOf(SoftPlus))
    post.io.in.valid     := infoAtPost.valid && zresAtPostOf(SoftPlus).valid
    post.io.in.bits.info := infoAtPost.bits
    post.io.in.bits.zres := zresAtPostOf(SoftPlus).bits

    zs(SoftPlus) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ScaleMixtureGaussian -- uses BOTH the pre-process and post-process
  // multipliers.
  //   pre-mult: (1/sgmA^2).manW1 * x.manW1   (= x/sgmA^2)
  //   post-mult: z0.manW1 * xsgmA2.manW1    (where z0 = polynomial*2^maxBitDigit + 1)

  if(cfg.has(ScaleMixtureGaussian)) {
    assert(preMulProd.isDefined, "ScaleMixtureGaussian requires the pre-process multiplier")
    assert(postMul.isDefined, "ScaleMixtureGaussian requires the post-process multiplier")
    val (sgmA, sgmB) = cfg.scaleMixtureGaussianSigma.get

    val pre  = toModule[ScaleMixtureGaussianPreProcess ](pres (ScaleMixtureGaussian))
    val post = toModule[ScaleMixtureGaussianPostProcess](posts(ScaleMixtureGaussian))

    // ----- pre-mult operands: (1/sgmA^2) * (1 + x.man) -----
    preMulArgs.get(ScaleMixtureGaussian).valid    := xIn.valid
    preMulArgs.get(ScaleMixtureGaussian).bits.lhs :=
      ScaleMixtureGaussianPreMulArgs.lhs(sgmA, spec, preMulLhsW)
    preMulArgs.get(ScaleMixtureGaussian).bits.rhs := Cat(1.U(1.W), xIn.bits.man)

    // ----- pre stage: x and sel on the far side of the pre-mult -----
    val smgSelAtPre = selPreOf(ScaleMixtureGaussian) === cfg.signal(ScaleMixtureGaussian)
    pre.io.x.valid := xPreOf(ScaleMixtureGaussian).valid && smgSelAtPre
    pre.io.x.bits  := xPreOf(ScaleMixtureGaussian).bits
    pre.io.prod.valid := preMulValid.get && smgSelAtPre
    pre.io.prod.bits  := ScaleMixtureGaussianPreMulArgs.prod(spec, preMulProd.get)

    // ----- table -----
    tableAdrs(ScaleMixtureGaussian) := pre.io.adr
    if (bitW.nOrder != 0) {
      polynomialDxs.get(ScaleMixtureGaussian) := pre.io.dx.get
    }

    // ----- post-mult args at the polynomial output -----
    val infoAtPost = Pipe(pre.io.info, nTabToPoly)
    val mulArgs    = Module(new ScaleMixtureGaussianPostMulArgs(sgmA, sgmB, spec, bitW))
    mulArgs.io.useTable  := infoAtPost.bits.useTable
    mulArgs.io.xsgmA2Man := infoAtPost.bits.xsgmA2Man
    mulArgs.io.zres      := polynomialEval.io.out.bits
    mulArgs.io.valid     := infoAtPost.valid && polynomialEval.io.out.valid

    postMulArgs.get(ScaleMixtureGaussian).valid    := mulArgs.io.mul.valid
    postMulArgs.get(ScaleMixtureGaussian).bits.lhs := mulArgs.io.mul.bits.lhs
    postMulArgs.get(ScaleMixtureGaussian).bits.rhs := mulArgs.io.mul.bits.rhs

    // ----- post stage: pipe info & z0Ex to the input of the post-process -----
    val infoAtAfterMul = Pipe(infoAtPost, nMulToPostOf(ScaleMixtureGaussian))
    val z0exAtAfterMul = Pipe(mulArgs.io.mul.valid, mulArgs.io.mul.bits.z0ex,
                              nMulToPostOf(ScaleMixtureGaussian))

    post.io.info := infoAtAfterMul
    post.io.z0ex := z0exAtAfterMul
    post.io.mul  := postMul.get.io.out
    zs(ScaleMixtureGaussian) := post.io.z
  }

  // --------------------------------------------------------------------------
  // Log2Fixed -- one region of the table per binade of the distance of the
  // argument from one, on each side of it, each holding |log2 x| itself, so
  // the result assembly normalizes it and neither shared multiplier takes
  // part.

  if(cfg.has(Log2Fixed)) {
    val pre  = toModule[Log2FixedPreProcess ](pres (Log2Fixed))
    val post = toModule[Log2FixedPostProcess](posts(Log2Fixed))

    pre.io.x.valid :=
      xPreOf(Log2Fixed).valid && selPreOf(Log2Fixed) === cfg.signal(Log2Fixed)
    pre.io.x.bits  := xPreOf(Log2Fixed).bits

    // ----- table -----
    tableAdrs(Log2Fixed) := pre.io.adr

    if(bitW.nOrder != 0) {
      polynomialDxs.get(Log2Fixed) := pre.io.dx.get
    }

    // ----- post -----
    post.io.pre  := Pipe(pre.io.info, nTabToPostOf(Log2Fixed))
    post.io.zres := zresAtPostOf(Log2Fixed)

    zs(Log2Fixed) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ATanFixed -- one region of the table per binade of |t|, each holding
  // atan(t) itself, so the result assembly normalizes it and neither shared
  // multiplier takes part.

  if(cfg.has(ATanFixed)) {
    val pre  = toModule[ATanFixedPreProcess ](pres (ATanFixed))
    val post = toModule[ATanFixedPostProcess](posts(ATanFixed))

    pre.io.x.valid := xPreOf(ATanFixed).valid && selPreOf(ATanFixed) === cfg.signal(ATanFixed)
    pre.io.x.bits  := xPreOf(ATanFixed).bits

    tableAdrs(ATanFixed) := pre.io.adr

    if (bitW.nOrder != 0) {
      polynomialDxs.get(ATanFixed) := pre.io.dx.get
    }

    post.io.pre  := Pipe(pre.io.info, nTabToPostOf(ATanFixed))
    post.io.zres := zresAtPostOf(ATanFixed)

    zs(ATanFixed) := post.io.z
  }

  // --------------------------------------------------------------------------
  // ACosFixed -- one region of the table per binade of 1 - |x|, each holding
  // the angle itself, so the result assembly normalizes it and neither shared
  // multiplier takes part.

  if(cfg.has(ACosFixed)) {
    val pre  = toModule[ACosFixedPreProcess ](pres (ACosFixed))
    val post = toModule[ACosFixedPostProcess](posts(ACosFixed))

    pre.io.x.valid := xPreOf(ACosFixed).valid && selPreOf(ACosFixed) === cfg.signal(ACosFixed)
    pre.io.x.bits  := xPreOf(ACosFixed).bits

    tableAdrs(ACosFixed) := pre.io.adr

    if (bitW.nOrder != 0) {
      polynomialDxs.get(ACosFixed) := pre.io.dx.get
    }

    post.io.pre  := Pipe(pre.io.info, nTabToPostOf(ACosFixed))
    post.io.zres := zresAtPostOf(ACosFixed)

    zs(ACosFixed) := post.io.z
  }
}

/** Input and output bundles and pipeline stage points of
 * [[rial.math.MathFunctions]].
 *
 * @see [[rial.util.StagePoint]]
 * @see [[rial.util.StageGroup]]
 */
object MathFunctions {

  /** Input of [[rial.math.MathFunctions]]: the selected function and its
   * argument(s).
   *
   * `sel` is a signal built with [[rial.math.MathFuncConfig.signal]]; `y` is
   * present only when [[rial.math.FuncKind.ATan2Phase1]] is configured, and
   * `tag` only when the module was given a tag width. The tag comes back with
   * the result on [[rial.math.MathFunctions.Output]].
   *
   * The arguments are in the floating point format of `cfg`.
   *
   * @param tagW width of the tag that rides along with the operation; 0 for
   *             no tag
   */
  class Input(val cfg: MathFuncConfig, val tagW: Int = 0) extends Bundle {
    val spec = cfg.spec
    val sel = UInt(cfg.signalW.W)
    val x   = UInt(spec.W.W)
    val y   = if (cfg.has(FuncKind.ATan2Phase1)) {
      Some(UInt(spec.W.W))
    } else {
      None
    }
    val tag = if (tagW > 0) {
      Some(UInt(tagW.W))
    } else {
      None
    }
  }

  /** Output of [[rial.math.MathFunctions]]: the result of the function and the
   * tag it was issued with.
   *
   * `tag` is present only when the module was given a tag width; it is the
   * `tag` of [[rial.math.MathFunctions.Input]] the result belongs to.
   *
   * @param spec floating point spec of the result
   * @param tagW width of the tag that rides along with the operation; 0 for
   *             no tag
   */
  class Output(val spec: RealSpec, val tagW: Int = 0) extends Bundle {
    val z   = UInt(spec.W.W)
    val tag = if (tagW > 0) {
      Some(UInt(tagW.W))
    } else {
      None
    }
  }

  /** The operands one function offers to a shared multiplier. */
  class MulOperands(val lhsW: Int, val rhsW: Int) extends Bundle {
    val lhs = UInt(lhsW.W)
    val rhs = UInt(rhsW.W)
  }

  /** Boundaries at which pipeline registers can be placed; see
   * [[rial.math.MathFunctions]] for the points and their names.
   *
   * `"preMul"`, `"mulArgs"` and `"postMul"` are declared only when one of the
   * functions of `cfg` needs the corresponding multiplier.
   *
   * @param preMulImpl  structure of the pre-process multiplier's product
   * @param postMulImpl structure of the post-process multiplier's product
   * @param polyMulImpl structure of the polynomial's products
   * @param perFunctionLatency whether each function keeps the latency of its
   *                    own path, in which case the paths of `"pre"` and
   *                    `"post"` are not padded to a common latency
   */
  def stagePoints(cfg: MathFuncConfig,
    preMulImpl:  MulImpl = MulImpl.plain,
    postMulImpl: MulImpl = MulImpl.plain,
    polyMulImpl: MulImpl = MulImpl.plain,
    perFunctionLatency: Boolean = false): Seq[StageElem] = {
    val spec = cfg.spec
    val preMul: Seq[StageElem] = if (PreProcMultiplier.isNeeded(cfg.funcs)) {
      Seq(StageGroup("preMul", if (preMulImpl.isSplit) {
        SplitMultiplier.stagePoints(SplitMultiplier.chainPoints(
          PreProcMultiplier.lhsW(cfg.funcs, spec), 1 + spec.manW,
          preMulImpl.sliceW.get))
      } else {
        Seq(StagePoint("out", 4.0))
      }))
    } else {
      Seq.empty[StageElem]
    }
    val postMul: Seq[StageElem] = if (PostProcMultiplier.isNeeded(cfg.funcs)) {
      Seq(StageGroup("mulArgs",
            Seq(StagePoint("form", 1.0), StagePoint("out", 0.5))),
          StageGroup("postMul", PostProcMultiplier.stagePoints(postMulImpl)))
    } else {
      Seq.empty[StageElem]
    }
    preMul ++
    Seq(StageParallel("pre", prePaths(cfg), !perFunctionLatency),
        StagePoint("table", 2.0),
        StagePoint("coef",  1.0),
        StageGroup("poly",  PolynomialEval.stagePoints(cfg.bitW.nOrder,
          polyMulImpl))) ++
    postMul ++
    Seq(StageParallel("post", postPaths(cfg), !perFunctionLatency),
        StagePoint("out",  1.0))
  }

  /** The argument reduction of each configured function, in `cfg.funcs` order,
   *  named with [[rial.math.FuncKind.pointName]]: the functions that share a
   *  name contribute one path. */
  def prePaths(cfg: MathFuncConfig): Seq[StageGroup] = {
    import FuncKind._
    cfg.funcs.collect {
      case Sqrt        => path(Sqrt,        SqrtPreProcess               .stagePoints)
      case InvSqrt     => path(InvSqrt,     InvSqrtPreProcess            .stagePoints)
      case Reciprocal  => path(Reciprocal,  ReciprocalPreProcess         .stagePoints)
      case Sin         => path(Sin,         SinCosPreProcess             .stagePoints)
      case SinPi       => path(SinPi,       SinCosPreProcess             .stagePoints)
      case SinPiFixed  => path(SinPiFixed,  SinPiFixedPreProcess         .stagePoints)
      case ACosPhase1  => path(ACosPhase1,  ACosPhase1PreProcess         .stagePoints)
      case ACosPhase2  => path(ACosPhase2,  ACosPhase2PreProcess         .stagePoints)
      case ATan2Phase1 => path(ATan2Phase1, ATan2Phase1PreProcess        .stagePoints)
      case ATan2Phase2 => path(ATan2Phase2, ATan2Phase2PreProcess        .stagePoints)
      case Exp         => path(Exp,         ExpPreProcess                .stagePoints)
      case Pow2        => path(Pow2,        ExpPreProcess                .stagePoints)
      case Log         => path(Log,         LogPreProcess                .stagePoints)
      case Log2Fixed   => path(Log2Fixed,   Log2FixedPreProcess          .stagePoints)
      case Sigmoid     => path(Sigmoid,     SigmoidPreProcess            .stagePoints)
      case SoftPlus    => path(SoftPlus,    SoftPlusPreProcess           .stagePoints)
      case ScaleMixtureGaussian =>
                          path(ScaleMixtureGaussian,
                               ScaleMixtureGaussianPreProcess.stagePoints)
      case ATanFixed   => path(ATanFixed,   ATanFixedPreProcess          .stagePoints)
      case ACosFixed   => path(ACosFixed,   ACosFixedPreProcess          .stagePoints)
    }
  }

  /** The result assembly of each configured function, in `cfg.funcs` order,
   *  named with [[rial.math.FuncKind.pointName]]: the functions that share a
   *  name contribute one path. */
  def postPaths(cfg: MathFuncConfig): Seq[StageGroup] = {
    import FuncKind._
    cfg.funcs.collect {
      case Sqrt        => path(Sqrt,        SqrtPostProcess               .stagePoints)
      case InvSqrt     => path(InvSqrt,     InvSqrtPostProcess            .stagePoints)
      case Reciprocal  => path(Reciprocal,  ReciprocalPostProcess         .stagePoints)
      case Sin         => path(Sin,         SinCosPostProcess             .stagePoints)
      case SinPi       => path(SinPi,       SinCosPostProcess             .stagePoints)
      case SinPiFixed  => path(SinPiFixed,  SinPiFixedPostProcess         .stagePoints)
      case ACosPhase1  => path(ACosPhase1,  ACosPhase1PostProcess         .stagePoints)
      case ACosPhase2  => path(ACosPhase2,  ACosPhase2PostProcess         .stagePoints)
      case ATan2Phase1 => path(ATan2Phase1, ATan2Phase1PostProcess        .stagePoints)
      case ATan2Phase2 => path(ATan2Phase2, ATan2Phase2PostProcess        .stagePoints)
      case Exp         => path(Exp,         ExpPostProcess                .stagePoints)
      case Pow2        => path(Pow2,        ExpPostProcess                .stagePoints)
      case Log         => path(Log,         LogPostProcess                .stagePoints)
      case Log2Fixed   => path(Log2Fixed,   Log2FixedPostProcess          .stagePoints)
      case Sigmoid     => path(Sigmoid,     SigmoidPostProcess            .stagePoints)
      case SoftPlus    => path(SoftPlus,    SoftPlusPostProcess           .stagePoints)
      case ScaleMixtureGaussian =>
                          path(ScaleMixtureGaussian,
                               ScaleMixtureGaussianPostProcess.stagePoints)
      case ATanFixed   => path(ATanFixed,   ATanFixedPostProcess          .stagePoints)
      case ACosFixed   => path(ACosFixed,   ACosFixedPostProcess          .stagePoints)
    }
  }

  /** One path of the `"pre"` or `"post"` parallel, named after its function. */
  private def path(fn: FuncKind.FuncKind, points: Seq[StageElem]): StageGroup = {
    StageGroup(FuncKind.pointName(fn), points)
  }
}
