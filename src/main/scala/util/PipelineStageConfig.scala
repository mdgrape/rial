package rial.util

import scala.language.implicitConversions

/** An element of the ordered list of pipeline boundaries a module declares.
 *
 * A module declares its boundaries as an ordered `Seq[StageElem]` whose last
 * element must be a [[rial.util.StagePoint]] named `"out"` (the module
 * output). The point `"in"` (the module input) is implicit and always
 * available; neither `"in"` nor `"out"` may name an internal boundary.
 *
 * An element is either a plain boundary ([[rial.util.StagePoint]]) or the
 * whole boundary list of a submodule ([[rial.util.StageGroup]]). A group is
 * flattened by [[rial.util.PipelineStageConfig.resolve]] into dot-separated
 * names, so a name must not contain `'.'`.
 */
sealed trait StageElem {

  /** Name of the element, unique among its siblings. */
  def name: String
}

/** A named boundary in a module's datapath at which pipeline registers may be
 *  inserted.
 *
 * `weight` is the relative delay of the logic between the previous point (or
 * the module input, for the first point) and this point. It is used only by
 * the automatic placement in [[rial.util.PipelineStageConfig.resolve]].
 *
 * @param name   name of the boundary, as used by `PipelineStageConfig.specify`
 * @param weight relative delay of the logic ending at this boundary
 */
case class StagePoint(name: String, weight: Double) extends StageElem

/** The boundary list of a submodule, included in the parent's list.
 *
 * A module that instantiates a submodule on its datapath declares the
 * submodule's own point list as a group, so that the automatic placement sees
 * the whole hierarchy at once. [[rial.util.PipelineStageConfig.resolve]]
 * flattens a group into dot-separated names: a group named `"poly"` whose
 * child list is `Seq(StagePoint("step1", w1), StagePoint("out", w2))`
 * contributes, in this order, the points `"poly.in"`, `"poly.step1"` and
 * `"poly.out"`.
 *
 * `"poly.out"` is an ordinary internal boundary of the parent. `"poly.in"`
 * sits at the group's input, where the parent's preceding boundary already is;
 * it can be named explicitly but is never chosen by the automatic placement.
 *
 * The child list follows the same rule as a module's own list: its last
 * element must be a [[rial.util.StagePoint]] named `"out"`, and its `"in"` is
 * implicit. Groups may be nested.
 *
 * A resolved configuration is handed down to the submodule with
 * [[rial.util.PipelineStages.sub]].
 *
 * @param name   name of the group, used as the prefix of its children
 * @param points the submodule's own ordered boundary list
 */
case class StageGroup(name: String, points: Seq[StageElem]) extends StageElem

/** A boundary of the flattened (fully dotted) point list.
 *
 * @param name       full dot-separated name
 * @param weight     relative delay of the logic ending at this boundary
 * @param selectable whether the automatic placement may cut here; the `"in"`
 *                   pseudo-point of a group is not selectable
 */
private[util] case class FlatPoint(name: String, weight: Double, selectable: Boolean)

/** Resolved number of pipeline registers per named point, owned by the module.
 *
 * Obtained from [[rial.util.PipelineStageConfig.resolve]]. A module inserts
 * `ShiftRegister(sig, stages.at("<point>"))` on every signal that crosses the
 * boundary `"<point>"`, and hands a submodule its share of the configuration
 * with [[sub]].
 *
 * @param counts number of registers at each point, keyed by full dotted name;
 *               contains an entry for `"in"` and for every flattened point
 * @param points the ordered boundary list the counts were resolved against
 */
case class PipelineStages(counts: Map[String, Int], points: Seq[StageElem]) {

  private lazy val names: Seq[String] = {
    PipelineStageConfig.validNames(points)
  }

  private lazy val aliases: Map[String, String] = {
    PipelineStageConfig.groupAliases(points)
  }

  /** Number of registers at the named point.
   *
   * The name of a group stands for the group's `"out"`, so `at("poly")` and
   * `at("poly.out")` are the same boundary.
   *
   * @throws IllegalArgumentException if `name` is not a declared point.
   */
  def at(name: String): Int = {
    val key = aliases.getOrElse(name, name)
    counts.getOrElse(key,
      throw new IllegalArgumentException(
        s"unknown pipeline stage point '$name'; " +
        s"valid names are ${names.mkString(", ")}"))
  }

  /** Total number of pipeline registers, i.e. the latency of the module. */
  def total: Int = {
    counts.values.sum
  }

  /** The part of this result that belongs to the group `name`.
   *
   * The counts of the group's children, with the `"name."` prefix stripped,
   * as an explicit configuration with no automatic register left. A parent
   * elaborates its submodule with `Module(new Child(..., stages.sub("poly")))`;
   * resolving that configuration against the child's own point list reproduces
   * exactly the placement the parent decided.
   *
   * @throws IllegalArgumentException if `name` is not a declared group.
   */
  def sub(name: String): PipelineStageConfig = {
    val prefix = name + "."
    val groups = PipelineStageConfig.groupNames(points)
    if (!groups.contains(name)) {
      val known = if (groups.isEmpty) {
        "none"
      } else {
        groups.mkString(", ")
      }
      throw new IllegalArgumentException(
        s"unknown pipeline stage group '$name'; declared groups are $known")
    }
    PipelineStageConfig(explicit = counts.collect {
      case (k, v) if k.startsWith(prefix) => k.drop(prefix.length) -> v
    })
  }

  /** Human readable form, e.g. `"in:0 unpack:0 swap:1 align:0 out:0"`. */
  def getString: String = {
    names.map(n => s"$n:${counts.getOrElse(n, 0)}").mkString(" ")
  }
}

/** Pipeline configuration given by the user of a module.
 *
 * A configuration is a number of registers at explicitly named boundaries
 * (`explicit`) plus a number of registers whose position is left to the module
 * (`auto`). The module resolves the configuration against the boundaries it
 * declares, via [[resolve]]; the automatic registers are then placed so that
 * the longest remaining piece of combinational logic is as short as possible.
 *
 * {{{
 * // three registers, placed by the module
 * new AddFPGeneric(spec, spec, spec, round, 3)
 * new AddFPGeneric(spec, spec, spec, round, PipelineStageConfig.default(3))
 *
 * // three registers at the output / two registers at the input
 * new AddFPGeneric(spec, spec, spec, round, PipelineStageConfig.atOut(3))
 * new AddFPGeneric(spec, spec, spec, round, PipelineStageConfig.atIn(2))
 *
 * // hand-picked positions
 * new AddFPGeneric(spec, spec, spec, round,
 *   PipelineStageConfig.specify("add" -> 1, "norm" -> 1, "round" -> 1))
 *
 * // combinational
 * new AddFPGeneric(spec, spec, spec, round, PipelineStageConfig.none)
 * }}}
 *
 * A module whose datapath contains a submodule declares the submodule's points
 * as a [[rial.util.StageGroup]]; the names of the group's children are then
 * dot-separated, and the whole hierarchy is balanced at once:
 *
 * {{{
 * // the parent's list
 * val stagePoints = Seq(
 *   StagePoint("addr", 2),
 *   StageGroup("poly", PolynomialEval.stagePoints(nOrder)),
 *   StagePoint("out", 1))
 *
 * val stages = stage.resolve(stagePoints)
 * val poly   = Module(new PolynomialEval(bitW, stages.sub("poly")))
 *
 * // a register inside the submodule, chosen by the user of the parent
 * PipelineStageConfig.specify("poly.step1" -> 1)
 * // two registers at the submodule's output
 * PipelineStageConfig.specify("poly" -> 2)
 * }}}
 *
 * A module that does not declare any boundary puts all the registers at its
 * output, i.e. it treats any configuration as "[[total]] registers at the
 * output". Only [[explicit]] names that a module actually declares are
 * accepted by [[resolve]]; consult the module's scaladoc for its list of
 * points.
 *
 * @param explicit number of registers at explicitly named points; the name
 *                 `"in"` (module input) is always accepted, the other names
 *                 must be declared by the module. A name may be dotted, and
 *                 the name of a group stands for the group's `"out"`
 * @param auto     number of registers whose position is chosen by the module
 */
case class PipelineStageConfig(explicit: Map[String, Int] = Map.empty, auto: Int = 0) {

  /** Total number of pipeline registers, i.e. the latency of the module. */
  def total: Int = {
    auto + explicit.values.sum
  }

  /** Resolves this configuration against the boundaries declared by a module.
   *
   * `points` is the ordered list of the module's internal boundaries; its last
   * element must be a [[rial.util.StagePoint]] named `"out"`, `"in"` is
   * implicit and always present in the result. A [[rial.util.StageGroup]] is
   * flattened into dot-separated names first, so the placement sees the whole
   * hierarchy as one flat datapath.
   *
   * The [[explicit]] counts are taken as they are; the [[auto]] registers are
   * then placed at the internal points (everything but `"in"`, `"out"` and the
   * `"in"` pseudo-point of a group, and excluding the points already cut
   * explicitly) so that the longest stretch of logic between two registers is
   * as short as possible, and as evenly spread as possible among the
   * placements that achieve it. The result is deterministic. Registers at
   * `"in"` and `"out"` do not split any logic.
   *
   * If [[auto]] exceeds the number of still uncut internal points, every
   * internal point is cut and the surplus registers are stacked at the point
   * that ends the longest remaining segment. A point whose count is `n > 1` is
   * simply a `ShiftRegister(x, n)`.
   *
   * @throws IllegalArgumentException on a negative count, on a malformed
   *         `points` list, or if [[explicit]] names a point the module does not
   *         declare.
   */
  def resolve(points: Seq[StageElem]): PipelineStages = {
    PipelineStageConfig.checkPoints(points)

    val flat    = PipelineStageConfig.flatten(points)
    val valid   = PipelineStageConfig.validNames(points)
    val aliases = PipelineStageConfig.groupAliases(points)

    val unknown = explicit.keys.filterNot(k =>
      valid.contains(k) || aliases.contains(k)).toSeq.sorted
    if (unknown.nonEmpty) {
      val groups = PipelineStageConfig.groupNames(points)
      val hint = if (groups.isEmpty) {
        ""
      } else {
        s" A group name (${groups.mkString(", ")}) stands for the group's 'out'."
      }
      throw new IllegalArgumentException(
        s"unknown pipeline stage point(s) ${unknown.map("'" + _ + "'").mkString(", ")}; " +
        s"valid names are ${valid.mkString(", ")}." + hint + " " +
        "Use PipelineStageConfig.default(n) to let the module place n registers, " +
        "or PipelineStageConfig.atOut(n) / atIn(n) to put them at the module boundary.")
    }
    explicit.foreach { case (k, v) =>
      if (v < 0) {
        throw new IllegalArgumentException(
          s"negative number of pipeline registers ($v) at point '$k'")
      }
    }
    if (auto < 0) {
      throw new IllegalArgumentException(
        s"negative number of automatically placed pipeline registers ($auto)")
    }

    val nPoints = flat.length
    val weights = flat.map(_.weight).toVector
    // index of "out" is nPoints-1; the internal points are 0 until nPoints-1
    val internal = (0 until nPoints - 1)

    val counts = scala.collection.mutable.Map(valid.map(n => n -> 0): _*)
    explicit.foreach { case (k, v) =>
      val key = aliases.getOrElse(k, k)
      counts(key) = counts(key) + v
    }

    val alreadyCut = internal.filter(i => counts(flat(i).name) > 0).toSet
    val available  = internal.filter(i =>
      flat(i).selectable && !alreadyCut.contains(i))

    val nAuto = math.min(auto, available.length)
    val chosen = PipelineStageConfig.placeAuto(
      weights, nPoints, alreadyCut, available, nAuto)

    chosen.foreach(i => counts(flat(i).name) += 1)

    // surplus: every internal point is cut, stack the rest at the point that
    // ends the longest remaining segment
    val surplus = auto - nAuto
    if (surplus > 0) {
      val cuts = alreadyCut ++ chosen
      val segs = PipelineStageConfig.segments(weights, nPoints, cuts)
      // segment i ends at point `ends(i)`
      val ends = (cuts.toSeq.sorted :+ (nPoints - 1))
      val longest = segs.zipWithIndex.maxBy { case (w, i) => (w, -i) }._2
      counts(flat(ends(longest)).name) += surplus
    }

    PipelineStages(counts.toMap, points)
  }

  /** Human readable form, e.g. `"auto:2"`, `"add:1 norm:1 round:1"`. */
  def getString: String = {
    val a = if (auto > 0) {
      Seq(s"auto:$auto")
    } else {
      Seq.empty[String]
    }
    val e = explicit.toSeq.sortBy { case (k, _) =>
      val order = if (k == "in") {
        0
      } else if (k == "out") {
        2
      } else {
        1
      }
      (order, k)
    }.map { case (k, v) => s"$k:$v" }
    val s = a ++ e
    if (s.isEmpty) {
      "none"
    } else {
      s.mkString(" ")
    }
  }
}

/** Constructor shorthands and the automatic placement for
 *  [[rial.util.PipelineStageConfig]].
 */
object PipelineStageConfig {

  /** `n` registers, placed by the module.
   *
   * Also used implicitly, so a plain `Int` can be given wherever a
   * [[rial.util.PipelineStageConfig]] is expected.
   */
  implicit def default(n: Int): PipelineStageConfig = {
    PipelineStageConfig(auto = n)
  }

  /** Registers at explicitly named points, e.g.
   *  `specify("add" -> 1, "norm" -> 1, "round" -> 1)`.
   */
  def specify(kv: (String, Int)*): PipelineStageConfig = {
    specify(kv.toMap)
  }

  /** Registers at explicitly named points. */
  def specify(m: Map[String, Int]): PipelineStageConfig = {
    PipelineStageConfig(explicit = m)
  }

  /** `n` registers at the module output. */
  def atOut(n: Int): PipelineStageConfig = {
    specify("out" -> n)
  }

  /** `n` registers at the module input. */
  def atIn(n: Int): PipelineStageConfig = {
    specify("in" -> n)
  }

  /** No pipeline register: the module is purely combinational. */
  def none: PipelineStageConfig = {
    PipelineStageConfig()
  }

  /** The declared points as one flat list of dot-separated boundaries. */
  private[util] def flatten(points: Seq[StageElem]): Seq[FlatPoint] = {
    def walk(elems: Seq[StageElem], prefix: String): Seq[FlatPoint] = {
      elems.flatMap {
        case StagePoint(n, w) => {
          Seq(FlatPoint(prefix + n, w, true))
        }
        case StageGroup(n, ps) => {
          FlatPoint(prefix + n + ".in", 0.0, false) +: walk(ps, prefix + n + ".")
        }
      }
    }
    walk(points, "")
  }

  /** Full dotted names of every declared group, outermost first. */
  private[util] def groupNames(points: Seq[StageElem]): Seq[String] = {
    def walk(elems: Seq[StageElem], prefix: String): Seq[String] = {
      elems.flatMap {
        case _: StagePoint => {
          Seq.empty[String]
        }
        case StageGroup(n, ps) => {
          (prefix + n) +: walk(ps, prefix + n + ".")
        }
      }
    }
    walk(points, "")
  }

  /** Maps the name of a group onto the boundary it stands for, its `"out"`. */
  private[util] def groupAliases(points: Seq[StageElem]): Map[String, String] = {
    groupNames(points).map(g => g -> (g + ".out")).toMap
  }

  /** All the names accepted for a module declaring `points`: `"in"` first. */
  private[util] def validNames(points: Seq[StageElem]): Seq[String] = {
    "in" +: flatten(points).map(_.name)
  }

  /** Checks that a module's point list is well formed. */
  private[util] def checkPoints(points: Seq[StageElem]): Unit = {
    if (points.isEmpty) {
      throw new IllegalArgumentException(
        "a module's pipeline stage point list must not be empty; " +
        "its last element must be StagePoint(\"out\", w)")
    }
    val names = points.map(_.name)
    points.last match {
      case StagePoint("out", _) => {
        // well formed
      }
      case other => {
        throw new IllegalArgumentException(
          "the last pipeline stage point must be a StagePoint named 'out', " +
          s"got '${other.name}'")
      }
    }
    if (names.distinct.length != names.length) {
      throw new IllegalArgumentException(
        s"duplicated pipeline stage point name(s) in ${names.mkString(", ")}")
    }
    if (names.contains("in")) {
      throw new IllegalArgumentException(
        "'in' is implicit and must not be declared as a pipeline stage point")
    }
    if (names.init.contains("out")) {
      throw new IllegalArgumentException(
        "'out' must only be used for the last pipeline stage point")
    }
    names.foreach { n =>
      if (n.isEmpty || n.contains(".")) {
        throw new IllegalArgumentException(
          s"'.' separates the levels of a pipeline stage point hierarchy and " +
          s"must not appear in a point name, got '$n'")
      }
    }
    points.foreach {
      case _: StagePoint => {
        // nothing to check
      }
      case StageGroup(_, ps) => {
        checkPoints(ps)
      }
    }
  }

  /** Lengths of the segments of combinational logic left by `cuts`.
   *
   * `cuts` are indices of internal points (registers there split the logic);
   * the last segment always ends at the output.
   */
  private[util] def segments(weights: Seq[Double], nPoints: Int, cuts: Set[Int]): Seq[Double] = {
    val b = Seq.newBuilder[Double]
    var acc = 0.0
    for (i <- 0 until nPoints) {
      acc += weights(i)
      if (cuts.contains(i)) {
        b += acc
        acc = 0.0
      }
    }
    b += acc // segment ending at the output
    b.result()
  }

  /** Tolerance for comparing segment lengths, which are sums of weights. */
  private val tol = 1e-9

  /** Lexicographic comparison of two equally long cut index lists. */
  private def lexLess(a: Seq[Int], b: Seq[Int]): Boolean = {
    val n = math.min(a.length, b.length)
    var i = 0
    var res = a.length < b.length
    var decided = false
    while ((i < n) && !decided) {
      if (a(i) != b(i)) {
        res = a(i) < b(i)
        decided = true
      }
      i += 1
    }
    res
  }

  /** Places `nAuto` cuts among `available`, given the cuts already `fixed`.
   *
   * This is the linear partition problem: the points are a sequence of
   * weighted pieces of logic, a cut may only be made at a point of
   * `available` (and is already made at every point of `fixed`), and the
   * longest resulting segment is to be minimised. It is solved exactly by a
   * dynamic programme over (boundary, number of cuts used) in two passes: the
   * first computes the smallest achievable longest segment, the second picks,
   * among the placements that reach it, the one with the smallest sum of
   * squared segment lengths and, on a tie, the lexicographically first list of
   * cuts. The last point (the module output) always ends the last segment.
   *
   * @return the chosen cut indices, in increasing order
   */
  private[util] def placeAuto(weights: Seq[Double], nPoints: Int, fixed: Set[Int],
    available: Seq[Int], nAuto: Int): Seq[Int] = {

    if (nAuto == 0) {
      return Seq.empty[Int]
    }

    val prefix = Array.ofDim[Double](nPoints + 1)
    for (i <- 0 until nPoints) {
      prefix(i + 1) = prefix(i) + weights(i)
    }

    // boundary 0 is the module input, boundary j (1 <= j <= nCand) is the
    // candidate point cands(j-1), boundary nCand+1 is the module output
    val cands = (fixed.toSeq ++ available).distinct.sorted.toVector
    val nCand = cands.length
    val last  = nCand + 1

    def endOf(j: Int): Int = {
      if (j == 0) {
        0
      } else if (j == last) {
        nPoints
      } else {
        cands(j - 1) + 1
      }
    }
    def lengthOf(i: Int, j: Int): Double = {
      prefix(endOf(j)) - prefix(endOf(i))
    }
    def isFixed(j: Int): Boolean = {
      (j >= 1) && (j <= nCand) && fixed.contains(cands(j - 1))
    }
    // a cut at boundary j costs one automatic register unless it is already
    // fixed; the module output is not a cut at all
    def costOf(j: Int): Int = {
      if ((j == last) || isFixed(j)) {
        0
      } else {
        1
      }
    }
    // a transition i -> j must not skip a boundary that is cut anyway
    val lastFixed = Array.ofDim[Int](last + 1)
    var seen = 0
    for (j <- 0 to last) {
      lastFixed(j) = seen
      if (isFixed(j)) {
        seen = j
      }
    }

    val inf = Double.PositiveInfinity

    // pass 1: the smallest achievable longest segment
    val longest = Array.fill(last + 1, nAuto + 1)(inf)
    longest(0)(0) = 0.0
    for (j <- 1 to last) {
      val c0 = costOf(j)
      for (c <- c0 to nAuto) {
        var best = inf
        for (i <- lastFixed(j) until j) {
          val prev = longest(i)(c - c0)
          if (prev < inf) {
            val cur = math.max(prev, lengthOf(i, j))
            if (cur < best) {
              best = cur
            }
          }
        }
        longest(j)(c) = best
      }
    }
    val limit = longest(last)(nAuto)
    if (limit == inf) {
      throw new IllegalArgumentException(
        s"cannot place $nAuto pipeline registers on ${available.length} points")
    }

    // pass 2: among the placements reaching `limit`, the most even one
    val square = Array.fill(last + 1, nAuto + 1)(inf)
    val picked = Array.fill(last + 1, nAuto + 1)(Vector.empty[Int])
    square(0)(0) = 0.0
    for (j <- 1 to last) {
      val c0 = costOf(j)
      for (c <- c0 to nAuto) {
        var best: Double = inf
        var bestSel = Vector.empty[Int]
        var found = false
        for (i <- lastFixed(j) until j) {
          val prev = square(i)(c - c0)
          val len  = lengthOf(i, j)
          if ((prev < inf) && (len <= limit + tol)) {
            val cur = prev + len * len
            val sel = if (c0 == 0) {
              picked(i)(c - c0)
            } else {
              picked(i)(c - c0) :+ cands(j - 1)
            }
            val takes = !found || (cur < best - tol) ||
              ((cur <= best + tol) && lexLess(sel, bestSel))
            if (takes) {
              best = cur
              bestSel = sel
              found = true
            }
          }
        }
        square(j)(c) = best
        picked(j)(c) = bestSel
      }
    }
    picked(last)(nAuto)
  }
}
