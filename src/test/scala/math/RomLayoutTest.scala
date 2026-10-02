package rial.tests

import chisel3._
import chisel3.util.log2Ceil

import circt.stage.ChiselStage

import org.scalatest.funspec.AnyFunSpec

import rial.math._
import rial.math.FuncKind._
import rial.arith._
import rial.table.PolynomialGenerator
import rial.util.PipelineStageConfig

/** The coefficient tables of a [[rial.math.MathFuncConfig]] in the ROM they
 *  share: where each one lands, what the ROM holds, and the hardware built
 *  from it.
 */
class RomLayoutTest extends AnyFunSpec {

  /** the address width the table of each function is given in the
   *  configuration whose tables do not share one; every function left out
   *  keeps the default.
   */
  val widths = Map(
    Sqrt        -> 5,
    InvSqrt     -> 6,
    Reciprocal  -> 7,
    Exp         -> 6,
    Sin         -> 6,
    Log         -> 6,
    ACosPhase2  -> 7,
    ATan2Phase2 -> 6,
  )

  val defaultW = 8

  private def config(addrW: TableAddrWidth,
    funcs: Seq[FuncKind] = MathFuncConfig.standard): MathFuncConfig = {
    new MathFuncConfig(
      funcs   = funcs,
      spec    = RealSpec.Float32Spec,
      nOrder  = 2,
      addrW   = addrW,
      inBitW  = RealSpec.Float32Spec.manW - defaultW,
      outBitW = RealSpec.Float32Spec.manW + 3,
    )
  }

  /** the `standard` set, each table at a width of its own. */
  val perFunc = config(PerFuncAddrWidth(defaultW, widths))

  val configs = Seq(
    "standardFP32"      -> MathFuncConfig.standardFP32,
    "standardFixedFP32" -> MathFuncConfig.standardFixedFP32,
    "allFP32"           -> MathFuncConfig.allFP32,
    "standardBF16"      -> MathFuncConfig.standardBF16,
    "perFuncAddrWidth"  -> perFunc,
  )

  // ==========================================================================
  // The layout
  // ==========================================================================

  describe("the merged coefficient ROM of a configuration") {

    for ((name, cfg) <- configs) {

      val layout = cfg.romLayout

      it(f"holds the table of every function of ${name}%s once") {
        assert(layout.owners.distinct == layout.owners)
        assert(layout.owners.toSet ==
          cfg.funcs.map(fn => TableAddrWidth.owner(fn)).toSet)
      }

      it(f"starts the table of each function of ${name}%s at a multiple of " +
         f"its own rows") {
        for (fn <- cfg.funcs) {
          assert(layout.offset(fn) % layout.depth(fn) == 0,
            f"the table of ${FuncKind.getString(fn)} starts at row " +
            f"${layout.offset(fn)}%d, which its ${layout.depth(fn)}%d rows " +
            f"do not divide")
        }
      }

      it(f"gives the tables of ${name}%s rows of their own") {
        val ranges = layout.owners.map(fn =>
          (layout.offset(fn), layout.offset(fn) + layout.depth(fn)))
        for (lhs <- ranges; rhs <- ranges if lhs != rhs) {
          assert(lhs._2 <= rhs._1 || rhs._2 <= lhs._1,
            f"the rows ${lhs} and ${rhs} overlap")
        }
        for ((begin, end) <- ranges) {
          assert(0 <= begin && end <= layout.rows)
        }
      }

      it(f"holds as many rows as the tables of ${name}%s hold together") {
        assert(layout.rows == layout.owners.map(fn => layout.depth(fn)).sum)
        assert(layout.addrW == log2Ceil(layout.rows))
        assert(layout.packed.length == layout.rows)
      }

      it(f"holds each table of ${name}%s at the rows it starts at") {
        for (fn <- layout.owners) {
          val table = PolynomialGenerator.packTable(cfg.coeffs(fn))
          assert(table.length == layout.depth(fn))
          assert(layout.packed.slice(layout.offset(fn),
            layout.offset(fn) + layout.depth(fn)) == table,
            f"the table of ${FuncKind.getString(fn)} is not at row " +
            f"${layout.offset(fn)}%d of the ROM")
        }
      }

      it(f"reads a table of ${name}%s at the width of that table") {
        for (fn <- cfg.funcs) {
          assert((1 << layout.adrW(fn)) == layout.depth(fn))
          assert(layout.adrW(fn) <= layout.addrW)
        }
      }

      it(f"holds one table of ${name}%s alone at the width of that table") {
        for (fn <- layout.owners) {
          val one = new RomLayout(Seq(fn), cfg.coeffs)
          assert(one.addrW == layout.adrW(fn))
          assert(one.rows  == layout.depth(fn))
          assert(one.rowW  == layout.rowW)
          assert(one.offset(fn) == 0)
        }
      }
    }

    it("sends a function that reads another one's table to the same rows") {
      val layout = perFunc.romLayout
      assert(layout.offset(Cos        ) == layout.offset(Sin       ))
      assert(layout.offset(ACosPhase1 ) == layout.offset(Sqrt      ))
      assert(layout.offset(ATan2Phase1) == layout.offset(Reciprocal))
      assert(layout.depth (Cos        ) == layout.depth (Sin       ))
      assert(layout.depth (ACosPhase1 ) == layout.depth (Sqrt      ))
      assert(layout.depth (ATan2Phase1) == layout.depth (Reciprocal))
    }

    it("refuses a configuration that would hold one table twice") {
      assertThrows[IllegalArgumentException] {
        config(defaultW, FuncKind.normalize(Seq(Exp, Pow2)))
      }
      assertThrows[IllegalArgumentException] {
        config(defaultW, FuncKind.normalize(Seq(Sin, SinPi)))
      }
      assertThrows[IllegalArgumentException] {
        config(defaultW, FuncKind.normalize(Seq(Cos, CosPi)))
      }
    }
  }

  // ==========================================================================
  // The hardware built from it
  // ==========================================================================

  describe("the hardware that reads the merged coefficient ROM") {

    for ((name, cfg) <- configs) {
      it(f"elaborates every function of ${name}%s") {
        ChiselStage.emitCHIRRTL(new MathFunctions(cfg, MathTestUtil.allPoints(cfg)))
        ChiselStage.emitCHIRRTL(new MathFunctions(cfg, PipelineStageConfig.none))
      }
    }
  }
}
