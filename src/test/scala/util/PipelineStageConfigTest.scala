package rial.tests

import org.scalatest.funspec.AnyFunSpec

import rial.math.{MathFuncConfig, MathFunctions}
import rial.util.{PipelineStageConfig, PipelineStages, StageElem, StageGroup, StageParallel, StagePoint}

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
    flatSegments(points.map(p => p.name -> p.weight), s)
  }

  // ==========================================================================
  // hierarchy
  // ==========================================================================

  /** The boundary list of a submodule, in the style of a polynomial. */
  val child: Seq[StagePoint] = Seq(
    StagePoint("step1", 3), StagePoint("step2", 3), StagePoint("out", 1))

  /** A parent that instantiates that submodule on its datapath. */
  val parent: Seq[StageElem] = Seq(
    StagePoint("addr", 2), StageGroup("poly", child), StagePoint("out", 2))

  /** The flattened form of [[parent]]: name and weight, in datapath order. */
  val parentFlat: Seq[(String, Double)] = Seq(
    "addr" -> 2.0, "poly.in" -> 0.0, "poly.step1" -> 3.0,
    "poly.step2" -> 3.0, "poly.out" -> 1.0, "out" -> 2.0)

  /** lengths of the segments left by the cuts of a flattened point list */
  def flatSegments(flat: Seq[(String, Double)], s: PipelineStages): Seq[Double] = {
    val b = Seq.newBuilder[Double]
    var acc = 0.0
    for ((n, w) <- flat) {
      acc += w
      if ((n != "out") && (s.at(n) > 0)) {
        b += acc
        acc = 0.0
      }
    }
    b += acc
    b.result()
  }

  /** names of the flattened points that carry a register, in datapath order */
  def flatCuts(flat: Seq[(String, Double)], s: PipelineStages): Seq[String] = {
    flat.map(_._1).filter(n => s.at(n) > 0)
  }

  describe("PipelineStageConfig.resolve: a group of points") {

    it("flattens a group into dot separated names") {
      val s = PipelineStageConfig.none.resolve(parent)
      assert(s.getString ==
        "in:0 addr:0 poly.in:0 poly.step1:0 poly.step2:0 poly.out:0 out:0")
      assert(s.total == 0)
    }

    it("balances the whole hierarchy at once") {
      val s1 = PipelineStageConfig.default(1).resolve(parent)
      // 2+0+3 | 3+1+2 is the best single cut
      assert(flatCuts(parentFlat, s1) == Seq("poly.step1"))
      assert(flatSegments(parentFlat, s1) == Seq(5.0, 6.0))

      val s2 = PipelineStageConfig.default(2).resolve(parent)
      assert(flatCuts(parentFlat, s2) == Seq("poly.step1", "poly.step2"))
      assert(flatSegments(parentFlat, s2) == Seq(5.0, 3.0, 3.0))

      val s3 = PipelineStageConfig.default(3).resolve(parent)
      assert(flatCuts(parentFlat, s3) ==
        Seq("addr", "poly.step1", "poly.step2"))
      assert(flatSegments(parentFlat, s3) == Seq(2.0, 3.0, 3.0, 3.0))
    }

    it("never picks the 'in' of a group automatically") {
      for (n <- 0 to 4) {
        val s = PipelineStageConfig.default(n).resolve(parent)
        assert(s.at("poly.in") == 0, s"poly.in was cut for n=$n")
      }
      // every other internal point is cut before the surplus is stacked
      val s = PipelineStageConfig.default(4).resolve(parent)
      assert(flatCuts(parentFlat, s) ==
        Seq("addr", "poly.step1", "poly.step2", "poly.out"))
    }

    it("stacks the surplus inside the group when it is the longest segment") {
      val s = PipelineStageConfig.default(5).resolve(parent)
      assert(s.total == 5)
      assert(s.at("poly.step1") == 2)
    }

    it("accepts a dotted name in specify") {
      val s = PipelineStageConfig.specify("poly.step2" -> 1, "addr" -> 1).resolve(parent)
      assert(s.total == 2)
      assert(flatCuts(parentFlat, s) == Seq("addr", "poly.step2"))
    }

    it("reads the name of a group as the group's out") {
      val s = PipelineStageConfig.specify("poly" -> 2).resolve(parent)
      assert(s.total == 2)
      assert(s.at("poly.out") == 2)
      assert(s.at("poly") == 2)
      assert(flatCuts(parentFlat, s) == Seq("poly.out"))
    }

    it("places the automatic registers around a register inside the group") {
      val c = PipelineStageConfig(Map("poly.step2" -> 1), auto = 1)
      val s = c.resolve(parent)
      assert(s.total == 2)
      assert(flatCuts(parentFlat, s) == Seq("poly.step1", "poly.step2"))
    }

    it("keeps total invariant for every n") {
      for (n <- 0 to 8) {
        val c = PipelineStageConfig.default(n)
        assert(c.resolve(parent).total == n, s"total mismatch for n=$n")
      }
    }

    it("lists the dotted names when a name is unknown") {
      val e = intercept[IllegalArgumentException] {
        PipelineStageConfig.specify("step1" -> 1).resolve(parent)
      }
      val m = e.getMessage
      for ((n, _) <- parentFlat) {
        assert(m.contains(n), s"message should list '$n': $m")
      }
      assert(m.contains("poly"))
    }

    it("rejects a malformed child list") {
      intercept[IllegalArgumentException] { // the child does not end with 'out'
        PipelineStageConfig.none.resolve(Seq(
          StageGroup("poly", Seq(StagePoint("step1", 1))), StagePoint("out", 1)))
      }
      intercept[IllegalArgumentException] { // '.' is reserved
        PipelineStageConfig.none.resolve(Seq(
          StagePoint("a.b", 1), StagePoint("out", 1)))
      }
    }
  }

  describe("PipelineStages.sub: handing a group down to a submodule") {

    it("reproduces the parent's decision in the submodule") {
      for (n <- 0 to 6) {
        val s = PipelineStageConfig.default(n).resolve(parent)
        val c = s.sub("poly")
        assert(c.auto == 0)
        val inner = c.resolve(child)
        for (p <- child) {
          assert(inner.at(p.name) == s.at("poly." + p.name),
            s"mismatch at '${p.name}' for n=$n")
        }
        assert(inner.at("in") == s.at("poly.in"))
        assert(inner.total == s.at("poly.in") + child.map(p => s.at("poly." + p.name)).sum)
      }
    }

    it("works through two levels of hierarchy") {
      val outer: Seq[StageElem] = Seq(
        StagePoint("pre", 2),
        StageGroup("a", Seq(StagePoint("x", 3), StageGroup("b", child),
          StagePoint("out", 1))),
        StagePoint("out", 1))
      val s = PipelineStageConfig.default(4).resolve(outer)
      assert(namesOf(s) == Seq("in", "pre", "a.in", "a.x", "a.b.in",
        "a.b.step1", "a.b.step2", "a.b.out", "a.out", "out"))
      val inner = s.sub("a")
      val innerPoints: Seq[StageElem] =
        Seq(StagePoint("x", 3), StageGroup("b", child), StagePoint("out", 1))
      val ri = inner.resolve(innerPoints)
      assert(ri.at("x") == s.at("a.x"))
      assert(ri.at("b.step1") == s.at("a.b.step1"))
      assert(ri.at("b.out") == s.at("a.b.out"))
      assert(ri.at("out") == s.at("a.out"))
      // and one level further down
      val leaf = ri.sub("b").resolve(child)
      assert(leaf.at("step1") == s.at("a.b.step1"))
      assert(leaf.at("out") == s.at("a.b.out"))
    }

    it("rejects an unknown group") {
      val s = PipelineStageConfig.default(2).resolve(parent)
      val e = intercept[IllegalArgumentException] {
        s.sub("nope")
      }
      assert(e.getMessage.contains("poly"))
    }
  }

  // ==========================================================================
  // the exact placement, checked against a brute force
  // ==========================================================================

  /** All the ways to place `n` cuts among `avail`, best first.
   *
   * The cost of a placement is (longest segment, sum of squared segment
   * lengths, cut indices), compared lexicographically: exactly what
   * [[rial.util.PipelineStageConfig.resolve]] is expected to minimise.
   */
  def bruteForce(weights: Seq[Double], forced: Set[Int], avail: Seq[Int],
    n: Int): (Double, Double, Seq[Int]) = {
    val nPoints = weights.length
    def segs(cuts: Set[Int]): Seq[Double] = {
      val b = Seq.newBuilder[Double]
      var acc = 0.0
      for (i <- 0 until nPoints) {
        acc += weights(i)
        if (cuts.contains(i)) {
          b += acc
          acc = 0.0
        }
      }
      b += acc
      b.result()
    }
    var best: (Double, Double, Seq[Int]) = null
    for (sel <- avail.combinations(n)) {
      val s = segs(forced ++ sel)
      val c = (s.max, s.map(x => x * x).sum, sel.sorted.toSeq)
      if ((best == null) || cheaper(c, best)) {
        best = c
      }
    }
    best
  }

  /** Lexicographic comparison of two brute force costs. */
  def cheaper(a: (Double, Double, Seq[Int]), b: (Double, Double, Seq[Int])): Boolean = {
    if (a._1 != b._1) {
      a._1 < b._1
    } else if (a._2 != b._2) {
      a._2 < b._2
    } else {
      a._3.zip(b._3).find { case (x, y) => x != y } match {
        case Some((x, y)) => {
          x < y
        }
        case None => {
          a._3.length < b._3.length
        }
      }
    }
  }

  /** The declared point names, in datapath order, from a resolved result. */
  def namesOf(s: PipelineStages): Seq[String] = {
    s.getString.split(" ").toSeq.map(t => t.substring(0, t.lastIndexOf(':')))
  }

  describe("PipelineStageConfig.resolve: exact placement") {

    it("matches a brute force on random point lists") {
      val rng = new scala.util.Random(20240517)
      for (trial <- 0 until 400) {
        val k = 2 + rng.nextInt(9) // 2 .. 10 points, 'out' included
        val weights = Seq.fill(k)(1.0 + rng.nextInt(5))
        val points = (0 until k - 1).map(i => StagePoint(s"p$i", weights(i))) :+
          StagePoint("out", weights(k - 1))
        val internal = (0 until k - 1)
        // cut a few points explicitly, place the rest automatically
        val forced = internal.filter(_ => rng.nextInt(4) == 0).toSet
        val avail  = internal.filterNot(forced.contains)
        val n = if (avail.isEmpty) {
          0
        } else {
          rng.nextInt(avail.length + 1)
        }
        val cfg = PipelineStageConfig(
          explicit = forced.map(i => points(i).name -> 1).toMap, auto = n)
        val s = cfg.resolve(points)
        assert(s.total == forced.size + n, s"total mismatch in trial $trial")
        val got = internal.filter(i => s.at(points(i).name) > 0)
          .filterNot(forced.contains)
        val (bestMax, bestSq, bestSel) = bruteForce(weights, forced, avail, n)
        val segs = flatSegments(points.map(p => p.name -> p.weight), s)
        assert(segs.max == bestMax,
          s"trial $trial: longest segment ${segs.max} != $bestMax " +
          s"(weights $weights, forced $forced, n=$n)")
        assert(segs.map(x => x * x).sum == bestSq,
          s"trial $trial: sum of squares mismatch (weights $weights)")
        assert(got == bestSel,
          s"trial $trial: cuts $got != $bestSel (weights $weights)")
      }
    }

    it("matches a brute force on random hierarchies") {
      val rng = new scala.util.Random(19660809)
      for (trial <- 0 until 200) {
        val nInner = 1 + rng.nextInt(4)
        val innerW = Seq.fill(nInner + 1)(1.0 + rng.nextInt(5))
        val inner: Seq[StagePoint] =
          (0 until nInner).map(i => StagePoint(s"q$i", innerW(i))) :+
            StagePoint("out", innerW(nInner))
        val nPre  = 1 + rng.nextInt(3)
        val preW  = Seq.fill(nPre)(1.0 + rng.nextInt(5))
        val postW = 1.0 + rng.nextInt(5)
        val points: Seq[StageElem] =
          (0 until nPre).map(i => StagePoint(s"p$i", preW(i))) ++
          Seq(StageGroup("g", inner), StagePoint("out", postW))

        val flatNames = (0 until nPre).map(i => s"p$i") ++ Seq("g.in") ++
          (0 until nInner).map(i => s"g.q$i") ++ Seq("g.out", "out")
        val flatW = preW ++ Seq(0.0) ++ innerW.init ++
          Seq(innerW(nInner), postW)
        val k = flatW.length
        // 'g.in' is never available to the automatic placement
        val avail = (0 until k - 1).filterNot(_ == nPre)

        val n = rng.nextInt(avail.length + 1)
        val s = PipelineStageConfig.default(n).resolve(points)
        assert(s.total == n)
        val got = avail.filter(i => s.at(flatNames(i)) > 0)
        val (bestMax, bestSq, bestSel) = bruteForce(flatW, Set.empty, avail, n)
        val segs = flatSegments(flatNames.zip(flatW), s)
        assert(segs.max == bestMax, s"trial $trial: $segs vs $bestMax")
        assert(segs.map(x => x * x).sum == bestSq, s"trial $trial: $segs")
        assert(got == bestSel, s"trial $trial: $got != $bestSel")
      }
    }
  }

  // ==========================================================================
  // parallel paths
  // ==========================================================================

  /** The boundary list of a short alternative path. */
  val fast: Seq[StageElem] = Seq(StagePoint("a", 3), StagePoint("out", 3))

  /** The boundary list of a long alternative path. */
  val slow: Seq[StageElem] = Seq(
    StagePoint("p", 4), StagePoint("q", 4), StagePoint("out", 2))

  /** A module whose datapath splits into those two paths. */
  val split: Seq[StageElem] = Seq(
    StagePoint("head", 2),
    StageParallel("mid", Seq(StageGroup("slow", slow), StageGroup("fast", fast))),
    StagePoint("out", 1))

  /** The names [[split]] declares, in datapath order. */
  val splitNames: Seq[String] = Seq(
    "in", "head", "mid.in", "mid.slow.in", "mid.slow.p", "mid.slow.q",
    "mid.slow.out", "mid.fast.in", "mid.fast.a", "mid.fast.out", "out")

  /** counts of the points that carry a register, in declared order */
  def splitCuts(s: PipelineStages): Seq[(String, Int)] = {
    splitNames.map(n => n -> s.at(n)).filter(_._2 > 0)
  }

  describe("PipelineStageConfig.resolve: parallel paths") {

    it("declares every point of every path") {
      val s = PipelineStageConfig.none.resolve(split)
      assert(namesOf(s) == splitNames)
      assert(s.total == 0)
      assert(s.latency("mid") == 0)
    }

    it("places the automatic registers on the longest path") {
      val s = PipelineStageConfig.default(1).resolve(split)
      // head 2 | mid.slow.p 4 | 4+2+1 : the best single cut of the long path
      assert(s.at("mid.slow.p") == 1)
      assert(s.latency("mid.slow") == 1)
    }

    it("pads the short path to the latency of the long one") {
      val s = PipelineStageConfig.default(1).resolve(split)
      assert(s.latency("mid.fast") == 1)
      assert(s.latency("mid") == 1)
      assert(s.total == 1)
      // the short path is balanced on its own points
      assert(s.at("mid.fast.a") == 1)
    }

    it("keeps every path at the same latency for every n") {
      for (n <- 0 to 10) {
        val s = PipelineStageConfig.default(n).resolve(split)
        assert(s.total == n, s"total mismatch for n=$n")
        assert(s.latency("mid.slow") == s.latency("mid.fast"),
          s"paths differ for n=$n: ${s.getString}")
        assert(s.latency("mid") == s.latency("mid.slow"))
      }
    }

    it("never cuts the short path in the placement of the parent") {
      // 'mid.fast.a' is not on the chain the placement sees, so an explicit
      // register there does not move the automatic one
      val a = PipelineStageConfig.default(1).resolve(split)
      val b = PipelineStageConfig(Map("mid.fast.a" -> 1), auto = 1).resolve(split)
      assert(a.at("mid.slow.p") == b.at("mid.slow.p"))
      assert(b.latency("mid") == 1)
      assert(b.total == 1)
    }

    it("raises the latency when a short path is cut explicitly") {
      val s = PipelineStageConfig.specify("mid.fast.a" -> 2).resolve(split)
      assert(s.latency("mid") == 2)
      assert(s.latency("mid.slow") == 2)
      assert(s.at("mid.slow.p") == 1)
      assert(s.at("mid.slow.q") == 1)
      assert(s.total == 2)
    }

    it("reads the name of a parallel as the out of every path") {
      val s = PipelineStageConfig.specify("mid" -> 2).resolve(split)
      assert(s.at("mid.slow.out") == 2)
      assert(s.at("mid.fast.out") == 2)
      assert(s.latency("mid") == 2)
      assert(s.total == 2)
      assert(splitCuts(s) == Seq("mid.slow.out" -> 2, "mid.fast.out" -> 2))
    }

    it("reads the in of a parallel as the in of every path") {
      val s = PipelineStageConfig.specify("mid.in" -> 1).resolve(split)
      assert(s.at("mid.slow.in") == 1)
      assert(s.at("mid.fast.in") == 1)
      assert(s.at("mid.in") == 0)
      assert(s.latency("mid") == 1)
      assert(s.total == 1)
    }

    it("gives latency for a point, a path and a parallel") {
      val s = PipelineStageConfig.default(3).resolve(split)
      assert(s.latency("head") == s.at("head"))
      assert(s.latency("mid.slow") ==
        Seq("mid.slow.in", "mid.slow.p", "mid.slow.q", "mid.slow.out")
          .map(s.at).sum)
      assert(s.latency("mid") == s.latency("mid.slow"))
      assert(s.total == s.at("in") + s.at("head") + s.latency("mid") + s.at("out"))
    }

    it("hands one path down to its submodule") {
      for (n <- 0 to 6) {
        val s = PipelineStageConfig.default(n).resolve(split)
        val c = s.sub("mid.fast")
        assert(c.auto == 0)
        val inner = c.resolve(fast)
        assert(inner.at("in")  == s.at("mid.fast.in"))
        assert(inner.at("a")   == s.at("mid.fast.a"))
        assert(inner.at("out") == s.at("mid.fast.out"))
        assert(inner.total == s.latency("mid"))
      }
    }

    it("prints every path") {
      val s = PipelineStageConfig.specify("mid" -> 1).resolve(split)
      assert(s.getString == "in:0 head:0 mid.in:0 mid.slow.in:0 mid.slow.p:0 " +
        "mid.slow.q:0 mid.slow.out:1 mid.fast.in:0 mid.fast.a:0 " +
        "mid.fast.out:1 out:0")
    }

    it("rejects the name of a parallel where a single point is expected") {
      val s = PipelineStageConfig.default(2).resolve(split)
      val e = intercept[IllegalArgumentException] {
        s.at("mid")
      }
      assert(e.getMessage.contains("latency"))
      val f = intercept[IllegalArgumentException] {
        s.sub("mid")
      }
      assert(f.getMessage.contains("mid.<path>"))
    }

    it("rejects an unknown name in latency") {
      val s = PipelineStageConfig.default(2).resolve(split)
      intercept[IllegalArgumentException] {
        s.latency("nope")
      }
    }

    it("rejects a malformed parallel") {
      intercept[IllegalArgumentException] { // no path
        PipelineStageConfig.none.resolve(
          Seq(StageParallel("p", Seq.empty), StagePoint("out", 1)))
      }
      intercept[IllegalArgumentException] { // duplicated path name
        PipelineStageConfig.none.resolve(Seq(
          StageParallel("p", Seq(StageGroup("a", fast), StageGroup("a", fast))),
          StagePoint("out", 1)))
      }
      intercept[IllegalArgumentException] { // a path named 'out'
        PipelineStageConfig.none.resolve(Seq(
          StageParallel("p", Seq(StageGroup("out", fast))), StagePoint("out", 1)))
      }
      intercept[IllegalArgumentException] { // a path that does not end with 'out'
        PipelineStageConfig.none.resolve(Seq(
          StageParallel("p", Seq(StageGroup("a", Seq(StagePoint("x", 1))))),
          StagePoint("out", 1)))
      }
      intercept[IllegalArgumentException] { // nested parallels
        val inner = StageParallel("q", Seq(StageGroup("a", fast)))
        PipelineStageConfig.none.resolve(Seq(
          StageParallel("p", Seq(StageGroup("a", Seq(inner, StagePoint("out", 1))))),
          StagePoint("out", 1)))
      }
    }
  }

  // ==========================================================================
  // parallel paths of their own latencies
  // ==========================================================================

  /** The same two paths, but each ending on its own boundary. */
  val ownSplit: Seq[StageElem] = Seq(
    StagePoint("head", 2),
    StageParallel("mid", Seq(StageGroup("slow", slow), StageGroup("fast", fast)),
      equalLatency = false),
    StagePoint("out", 1))

  /** Three paths of the same length but of very different granularity, each
   *  ending on its own boundary.
   */
  val grains: Seq[StageElem] = Seq(
    StageParallel("mid", Seq(
      StageGroup("long", Seq(StagePoint("p1", 3), StagePoint("out", 3))),
      StageGroup("fine", Seq(
        StagePoint("p1", 2), StagePoint("p2", 2), StagePoint("out", 2))),
      StageGroup("coarse", Seq(StagePoint("a", 5), StagePoint("out", 0)))),
      equalLatency = false),
    StagePoint("out", 0))

  describe("PipelineStageConfig.resolve: parallel paths of their own latencies") {

    it("asks for equal latencies unless told otherwise") {
      assert(StageParallel("mid", Seq(StageGroup("fast", fast))).equalLatency)
      val s = PipelineStageConfig.specify("mid.fast.a" -> 1).resolve(split)
      assert(s.latency("mid.slow") == s.latency("mid.fast"))
    }

    it("leaves the other paths alone when a path is cut explicitly") {
      val s = PipelineStageConfig.specify("mid.slow.p" -> 1).resolve(ownSplit)
      assert(s.latency("mid.slow") == 1)
      assert(s.latency("mid.fast") == 0)
      assert(s.total == 1)
    }

    it("counts the longest path in the total") {
      val s = PipelineStageConfig.specify("mid.fast.a" -> 2).resolve(ownSplit)
      assert(s.latency("mid.fast") == 2)
      assert(s.latency("mid.slow") == 0)
      assert(s.total == 2)
    }

    it("refuses to give one latency for the whole parallel") {
      val s = PipelineStageConfig.specify("mid.slow.p" -> 1).resolve(ownSplit)
      val e = intercept[IllegalArgumentException] {
        s.latency("mid")
      }
      assert(e.getMessage.contains("mid.<path>"))
      assert(e.getMessage.contains("mid.slow"))
    }

    it("reads the name of a parallel as the out of every path") {
      val s = PipelineStageConfig.specify("mid" -> 1).resolve(ownSplit)
      assert(s.at("mid.slow.out") == 1)
      assert(s.at("mid.fast.out") == 1)
      assert(s.latency("mid.slow") == 1)
      assert(s.latency("mid.fast") == 1)
      assert(s.total == 1)
    }

    it("leaves a path alone when its logic is already short enough") {
      // head 2 | mid.slow.p 4 | 4+2+1 : the bound is 7, the fast path is 6
      val s = PipelineStageConfig.default(1).resolve(ownSplit)
      assert(s.at("mid.slow.p") == 1)
      assert(s.latency("mid.slow") == 1)
      assert(s.latency("mid.fast") == 0)
      assert(s.total == 1)
    }

    it("cuts every other path down to the longest stretch of logic left") {
      // the long path is cut once, into 3 and 3, so the bound is 3
      val s = PipelineStageConfig.default(1).resolve(grains)
      assert(s.at("mid.long.p1") == 1)
      assert(s.latency("mid.long") == 1)
      // 2+2+2 needs two cuts to stay within 3, one would leave 4
      assert(s.at("mid.fine.p1") == 1)
      assert(s.at("mid.fine.p2") == 1)
      assert(s.latency("mid.fine") == 2)
      // 5 is longer than the bound, but there is only one point to cut at
      assert(s.at("mid.coarse.a") == 1)
      assert(s.latency("mid.coarse") == 1)
      assert(s.total == 2)
    }

    it("keeps an explicit register and adds only the missing ones") {
      val s = PipelineStageConfig(Map("mid.fine.p1" -> 1), auto = 1).resolve(grains)
      assert(s.at("mid.long.p1") == 1)
      assert(s.at("mid.fine.p1") == 1)
      assert(s.at("mid.fine.p2") == 1)
      assert(s.latency("mid.fine") == 2)
    }

    it("places nothing at all without an automatic register") {
      for (c <- Seq(PipelineStageConfig.none, PipelineStageConfig.default(0))) {
        val s = c.resolve(grains)
        assert(s.total == 0, s"${c.getString}: ${s.getString}")
        assert(s.latency("mid.long") == 0)
        assert(s.latency("mid.fine") == 0)
        assert(s.latency("mid.coarse") == 0)
      }
    }

    it("hands one path down to its submodule") {
      val s = PipelineStageConfig.default(1).resolve(grains)
      val c = s.sub("mid.fine")
      assert(c.auto == 0)
      val inner = c.resolve(Seq(
        StagePoint("p1", 2), StagePoint("p2", 2), StagePoint("out", 2)))
      assert(inner.total == s.latency("mid.fine"))
    }
  }

  // ==========================================================================
  // the parallels of the math datapath
  // ==========================================================================

  describe("MathFunctions.stagePoints: one path per function") {

    val cfg    = MathFuncConfig.standardFP32
    val points = MathFunctions.stagePoints(cfg)

    /** latency of every path of the parallel `name`. */
    def pathLatencies(s: PipelineStages, name: String,
      paths: Seq[StageGroup]): Seq[(String, Int)] = {
      paths.map(c => c.name -> s.latency(name + "." + c.name))
    }

    it("resolves for every n with every path at the same latency") {
      for (n <- 1 to 12) {
        val s = PipelineStageConfig.default(n).resolve(points)
        assert(s.total == n, s"total mismatch for n=$n: ${s.getString}")

        val pre = pathLatencies(s, "pre", MathFunctions.prePaths(cfg))
        assert(pre.map(_._2).distinct.length == 1,
          s"pre paths differ for n=$n: $pre")
        assert(s.latency("pre") == pre.head._2)

        val post = pathLatencies(s, "post", MathFunctions.postPaths(cfg))
        assert(post.map(_._2).distinct.length == 1,
          s"post paths differ for n=$n: $post")
        assert(s.latency("post") == post.head._2)
      }
    }

    it("hands each path down to its submodule") {
      val s = PipelineStageConfig.default(8).resolve(points)
      val sinCos = s.sub("pre.sinCos").resolve(rial.math.SinCosPreProcess.stagePoints)
      assert(sinCos.total == s.latency("pre"))
      assert(sinCos.at("out") == s.at("pre.sinCos.out"))
      val sqrt = s.sub("post.sqrt").resolve(rial.math.SqrtPostProcess.stagePoints)
      assert(sqrt.total == s.latency("post"))
    }

    it("adds a register to every path when the parallel is named") {
      val s = PipelineStageConfig.specify("pre" -> 1, "post" -> 1).resolve(points)
      assert(s.latency("pre") == 1)
      assert(s.latency("post") == 1)
      assert(s.total == 2)
      for (c <- MathFunctions.prePaths(cfg)) {
        assert(s.at("pre." + c.name + ".out") == 1)
      }
    }
  }
}
