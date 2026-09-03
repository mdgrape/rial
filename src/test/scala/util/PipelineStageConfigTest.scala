package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.util.{PipelineStageConfig, PipelineStages, StagePoint}

/** Pure-Scala unit tests of the pipeline stage placement.
 *
 * No RTL is elaborated here: [[rial.util.PipelineStageConfig.resolve]] is a
 * plain function from a configuration and a module's boundary list to a
 * register count per boundary.
 */
class PipelineStageConfigTest extends AnyFunSpec {

  /** The same boundary list as the generic floating point adder. */
  val adder: Seq[StagePoint] = Seq(
    StagePoint("unpack", 3), StagePoint("swap",  1), StagePoint("align", 3),
    StagePoint("add",    3), StagePoint("norm",  2), StagePoint("round", 4),
    StagePoint("out",    1))

  val allNames: Seq[String] = "in" +: adder.map(_.name)

  /** names of the points that carry at least one register, in declared order */
  def cutsOf(s: PipelineStages): Seq[String] = {
    allNames.filter(n => s.at(n) > 0)
  }

  /** counts of the points that carry at least one register, in declared order */
  def countsOf(s: PipelineStages): Seq[(String, Int)] = {
    allNames.map(n => n -> s.at(n)).filter(_._2 > 0)
  }

  /** lengths of the segments of logic left by the cuts, in datapath order */
  def segmentsOf(s: PipelineStages): Seq[Double] = {
    val b = Seq.newBuilder[Double]
    var acc = 0.0
    for (p <- adder) {
      acc += p.weight
      if ((p.name != "out") && (s.at(p.name) > 0)) {
        b += acc
        acc = 0.0
      }
    }
    b += acc
    b.result()
  }

  describe("PipelineStageConfig.resolve: automatic placement") {

    it("with no register leaves the datapath combinational") {
      val s = PipelineStageConfig.none.resolve(adder)
      assert(s.total == 0)
      assert(cutsOf(s).isEmpty)
      assert(segmentsOf(s) == Seq(17.0))
    }

    it("default(1) cuts the datapath in the middle") {
      val s = PipelineStageConfig.default(1).resolve(adder)
      assert(s.total == 1)
      assert(cutsOf(s) == Seq("align"))
      // 3+1+3 | 3+2+4+1 : the best achievable longest segment is 10
      assert(segmentsOf(s) == Seq(7.0, 10.0))
      assert(segmentsOf(s).max == 10.0)
    }

    it("default(2) balances three segments") {
      val s = PipelineStageConfig.default(2).resolve(adder)
      assert(s.total == 2)
      assert(cutsOf(s) == Seq("align", "norm"))
      assert(segmentsOf(s) == Seq(7.0, 5.0, 5.0))
    }

    it("default(3) balances four segments") {
      val s = PipelineStageConfig.default(3).resolve(adder)
      assert(s.total == 3)
      assert(cutsOf(s) == Seq("unpack", "align", "norm"))
      assert(segmentsOf(s) == Seq(3.0, 4.0, 5.0, 5.0))
    }

    it("default(4) balances five segments") {
      val s = PipelineStageConfig.default(4).resolve(adder)
      assert(s.total == 4)
      assert(cutsOf(s) == Seq("unpack", "align", "add", "norm"))
      assert(segmentsOf(s) == Seq(3.0, 4.0, 3.0, 2.0, 5.0))
    }

    it("default(5) leaves the shortest piece of logic uncut") {
      val s = PipelineStageConfig.default(5).resolve(adder)
      assert(s.total == 5)
      // 'swap' is the cheapest boundary, so it is the one that stays merged
      assert(cutsOf(s) == Seq("unpack", "align", "add", "norm", "round"))
      assert(segmentsOf(s) == Seq(3.0, 4.0, 3.0, 2.0, 4.0, 1.0))
    }

    it("default(6) covers every internal point") {
      val s = PipelineStageConfig.default(6).resolve(adder)
      assert(s.total == 6)
      assert(cutsOf(s) == Seq("unpack", "swap", "align", "add", "norm", "round"))
      assert(countsOf(s).forall(_._2 == 1))
      assert(s.at("in") == 0 && s.at("out") == 0)
    }

    it("default(7) stacks the surplus at the end of the longest segment") {
      val s = PipelineStageConfig.default(7).resolve(adder)
      assert(s.total == 7)
      // every internal point is cut; the longest remaining segment is 'round' (4)
      assert(s.at("round") == 2)
      assert(countsOf(s) == Seq("unpack" -> 1, "swap" -> 1, "align" -> 1,
        "add" -> 1, "norm" -> 1, "round" -> 2))
    }

    it("default(9) stacks all the surplus deterministically") {
      val s = PipelineStageConfig.default(9).resolve(adder)
      assert(s.total == 9)
      assert(s.at("round") == 4)
      assert(s.at("in") == 0 && s.at("out") == 0)
    }

    it("keeps total invariant for every n") {
      for (n <- 0 to 12) {
        val c = PipelineStageConfig.default(n)
        val s = c.resolve(adder)
        assert(s.total == c.total, s"total mismatch for n=$n")
        assert(s.total == n)
      }
    }

    it("is deterministic") {
      for (n <- 0 to 8) {
        val a = PipelineStageConfig.default(n).resolve(adder)
        val b = PipelineStageConfig.default(n).resolve(adder)
        assert(a.counts == b.counts)
      }
    }
  }

  describe("PipelineStageConfig.resolve: explicit placement") {

    it("atOut(3) puts every register at the output") {
      val s = PipelineStageConfig.atOut(3).resolve(adder)
      assert(s.total == 3)
      assert(s.at("out") == 3)
      assert(countsOf(s) == Seq("out" -> 3))
      // output registers do not split any logic
      assert(segmentsOf(s) == Seq(17.0))
    }

    it("atIn(2) puts every register at the input") {
      val s = PipelineStageConfig.atIn(2).resolve(adder)
      assert(s.total == 2)
      assert(s.at("in") == 2)
      assert(countsOf(s) == Seq("in" -> 2))
      assert(segmentsOf(s) == Seq(17.0))
    }

    it("specify places registers at the named points") {
      val s = PipelineStageConfig.specify("add" -> 1, "norm" -> 1, "round" -> 1).resolve(adder)
      assert(s.total == 3)
      assert(cutsOf(s) == Seq("add", "norm", "round"))
      assert(segmentsOf(s) == Seq(10.0, 2.0, 4.0, 1.0))
    }

    it("specify accepts a map as well as pairs") {
      val a = PipelineStageConfig.specify("add" -> 1, "norm" -> 2)
      val b = PipelineStageConfig.specify(Map("add" -> 1, "norm" -> 2))
      assert(a == b)
      assert(a.resolve(adder).counts == b.resolve(adder).counts)
    }

    it("allows more than one register at a point") {
      val s = PipelineStageConfig.specify("norm" -> 3).resolve(adder)
      assert(s.total == 3)
      assert(s.at("norm") == 3)
    }
  }

  describe("PipelineStageConfig.resolve: mixing explicit and automatic") {

    it("places the automatic registers around the explicit ones") {
      val c = PipelineStageConfig(Map("add" -> 1), auto = 1)
      assert(c.total == 2)
      val s = c.resolve(adder)
      assert(s.total == 2)
      // 'add' is already cut; the best remaining cut is 'swap' (4 | 6 | 7)
      assert(cutsOf(s) == Seq("swap", "add"))
      assert(segmentsOf(s) == Seq(4.0, 6.0, 7.0))
    }

    it("ignores input and output registers when balancing") {
      val c = PipelineStageConfig(Map("in" -> 1, "out" -> 2), auto = 1)
      val s = c.resolve(adder)
      assert(s.total == 4)
      assert(s.at("in") == 1)
      assert(s.at("out") == 2)
      // the automatic register lands where it would have without in/out ones
      assert(cutsOf(s) == Seq("in", "align", "out"))
    }

    it("keeps total invariant when mixing") {
      for (n <- 0 to 8) {
        val c = PipelineStageConfig(Map("in" -> 1, "align" -> 1, "out" -> 1), auto = n)
        assert(c.resolve(adder).total == c.total)
      }
    }
  }

  describe("PipelineStageConfig.resolve: error reporting") {

    it("rejects an unknown point name and lists the valid ones") {
      val e = intercept[IllegalArgumentException] {
        PipelineStageConfig.specify("shift" -> 1).resolve(adder)
      }
      val m = e.getMessage
      assert(m.contains("shift"))
      for (n <- allNames) {
        assert(m.contains(n), s"message should list '$n': $m")
      }
      assert(m.contains("default"))
      assert(m.contains("atOut"))
    }

    it("rejects a negative explicit count") {
      intercept[IllegalArgumentException] {
        PipelineStageConfig.specify("add" -> -1).resolve(adder)
      }
    }

    it("rejects a negative automatic count") {
      intercept[IllegalArgumentException] {
        PipelineStageConfig.default(-1).resolve(adder)
      }
    }

    it("rejects a malformed point list") {
      intercept[IllegalArgumentException] { // no 'out' at the end
        PipelineStageConfig.none.resolve(Seq(StagePoint("a", 1), StagePoint("b", 1)))
      }
      intercept[IllegalArgumentException] { // empty
        PipelineStageConfig.none.resolve(Seq.empty)
      }
      intercept[IllegalArgumentException] { // duplicated name
        PipelineStageConfig.none.resolve(
          Seq(StagePoint("a", 1), StagePoint("a", 1), StagePoint("out", 1)))
      }
      intercept[IllegalArgumentException] { // 'in' is implicit
        PipelineStageConfig.none.resolve(Seq(StagePoint("in", 1), StagePoint("out", 1)))
      }
      intercept[IllegalArgumentException] { // 'out' only at the end
        PipelineStageConfig.none.resolve(
          Seq(StagePoint("out", 1), StagePoint("a", 1), StagePoint("out", 1)))
      }
    }

    it("rejects an unknown name in PipelineStages.at") {
      val s = PipelineStageConfig.none.resolve(adder)
      val e = intercept[IllegalArgumentException] {
        s.at("nope")
      }
      assert(e.getMessage.contains("unpack"))
    }
  }

  describe("PipelineStageConfig: shorthands and strings") {

    it("converts an Int implicitly into an automatic configuration") {
      def f(c: PipelineStageConfig): Int = {
        c.auto
      }
      assert(f(3) == 3)
      assert(PipelineStageConfig.default(3) == PipelineStageConfig(auto = 3))
    }

    it("prints the configuration") {
      assert(PipelineStageConfig.none.getString == "none")
      assert(PipelineStageConfig.default(3).getString == "auto:3")
      assert(PipelineStageConfig.atOut(3).getString == "out:3")
      assert(PipelineStageConfig.atIn(2).getString == "in:2")
      assert(PipelineStageConfig.specify("add" -> 1, "norm" -> 1).getString == "add:1 norm:1")
      assert(PipelineStageConfig(Map("add" -> 1), auto = 2).getString == "auto:2 add:1")
    }

    it("prints the resolved stages in datapath order") {
      val s = PipelineStageConfig.default(2).resolve(adder)
      assert(s.getString == "in:0 unpack:0 swap:0 align:1 add:0 norm:1 round:0 out:0")
    }
  }

  describe("PipelineStageConfig.resolve: other point lists") {

    it("works for a module with a single point") {
      val one = Seq(StagePoint("out", 1))
      assert(PipelineStageConfig.default(2).resolve(one).at("out") == 2)
      assert(PipelineStageConfig.default(0).resolve(one).total == 0)
    }

    it("works for uniform weights") {
      val p = Seq.tabulate(8)(i => StagePoint(s"p$i", 1.0)) :+ StagePoint("out", 1.0)
      val s = PipelineStageConfig.default(2).resolve(p)
      assert(s.total == 2)
      assert(segmentsOf0(p, s) == Seq(3.0, 3.0, 3.0))
    }
  }

  /** [[segmentsOf]] for an arbitrary point list. */
  def segmentsOf0(points: Seq[StagePoint], s: PipelineStages): Seq[Double] = {
    val b = Seq.newBuilder[Double]
    var acc = 0.0
    for (p <- points) {
      acc += p.weight
      if ((p.name != "out") && (s.at(p.name) > 0)) {
        b += acc
        acc = 0.0
      }
    }
    b += acc
    b.result()
  }
}
