package rial.util

import scala.language.implicitConversions

/** A named boundary in a module's datapath at which pipeline registers may be
 *  inserted.
 *
 * `weight` is the relative delay of the logic between the previous point (or
 * the module input, for the first point) and this point. Weights are rough,
 * target-independent estimates; they are used only by the default (automatic)
 * placement in [[rial.util.PipelineStageConfig.resolve]] and never leave the
 * elaboration-time Scala code.
 *
 * A module declares its boundaries as an ordered `Seq[StagePoint]` whose last
 * element must be named `"out"` (the module output). The point `"in"` (the
 * module input) is implicit and always available; neither `"in"` nor `"out"`
 * may name an internal boundary.
 *
 * @param name   name of the boundary, as used by
 *               [[rial.util.PipelineStageConfig.specify]]
 * @param weight relative delay of the logic ending at this boundary
 */
case class StagePoint(name: String, weight: Double)

/** Resolved number of pipeline registers per named point, owned by the module.
 *
 * Obtained from [[rial.util.PipelineStageConfig.resolve]]. A module inserts
 * `ShiftRegister(sig, stages.at("<point>"))` on every signal that crosses the
 * boundary `"<point>"`.
 *
 * @param counts number of registers at each point; contains an entry for
 *               `"in"` and for every element of `points`
 * @param points the ordered boundary list the counts were resolved against
 */
case class PipelineStages(counts: Map[String, Int], points: Seq[StagePoint]) {

  /** Number of registers at the named point.
   *
   * @throws IllegalArgumentException if `name` is not a declared point.
   */
  def at(name: String): Int = counts.getOrElse(name,
    throw new IllegalArgumentException(
      s"unknown pipeline stage point '$name'; " +
      s"valid names are ${PipelineStageConfig.validNames(points).mkString(", ")}"))

  /** Total number of pipeline registers, i.e. the latency of the module. */
  def total: Int = counts.values.sum

  /** Human readable form, e.g. `"in:0 unpack:0 swap:1 align:0 out:0"`. */
  def getString: String =
    PipelineStageConfig.validNames(points).map(n => s"$n:${counts.getOrElse(n, 0)}").mkString(" ")
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
 * A module that does not declare any boundary (for example the multiplier, the
 * FMA, the fixed-point converters and the random number generators) simply
 * uses [[total]] and puts all the registers as a single `ShiftRegister` at its
 * output, leaving the placement to register balancing / retiming; such a module
 * therefore treats any configuration as "[[total]] registers at the output".
 * Only [[explicit]] names that a module actually declares are accepted by
 * [[resolve]]; consult the module's scaladoc for its list of points.
 *
 * @param explicit number of registers at explicitly named points; the name
 *                 `"in"` (module input) is always accepted, the other names
 *                 must be declared by the module
 * @param auto     number of registers whose position is chosen by the module
 */
case class PipelineStageConfig(explicit: Map[String, Int] = Map.empty, auto: Int = 0) {

  /** Total number of pipeline registers, i.e. the latency of the module. */
  def total: Int = auto + explicit.values.sum

  /** Resolves this configuration against the boundaries declared by a module.
   *
   * `points` is the ordered list of the module's internal boundaries; its last
   * element must be named `"out"`, `"in"` is implicit and always present in the
   * result. The [[explicit]] counts are taken as they are; the [[auto]]
   * registers are then distributed over the internal points (everything but
   * `"in"` and `"out"`, and excluding the points already cut explicitly) so
   * that the longest resulting segment of combinational logic - the sum of the
   * [[StagePoint.weight]]s between two consecutive cuts - is minimal. Ties are
   * broken by comparing the segment lengths sorted in decreasing order, and
   * then by the order in which the points are declared, so the result is
   * deterministic. Registers at `"in"` and `"out"` do not split any logic.
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
  def resolve(points: Seq[StagePoint]): PipelineStages = {
    PipelineStageConfig.checkPoints(points)

    val valid = PipelineStageConfig.validNames(points)
    val unknown = explicit.keys.filterNot(valid.contains).toSeq.sorted
    if (unknown.nonEmpty) {
      throw new IllegalArgumentException(
        s"unknown pipeline stage point(s) ${unknown.map("'" + _ + "'").mkString(", ")}; " +
        s"valid names are ${valid.mkString(", ")}. " +
        "Use PipelineStageConfig.default(n) to let the module place n registers, " +
        "or PipelineStageConfig.atOut(n) / atIn(n) to put them at the module boundary.")
    }
    explicit.foreach { case (k, v) =>
      if (v < 0) throw new IllegalArgumentException(
        s"negative number of pipeline registers ($v) at point '$k'")
    }
    if (auto < 0) throw new IllegalArgumentException(
      s"negative number of automatically placed pipeline registers ($auto)")

    val nPoints  = points.length
    val weights  = points.map(_.weight).toVector
    // index of "out" is nPoints-1; the internal points are 0 until nPoints-1
    val internal = (0 until nPoints - 1)

    val counts = scala.collection.mutable.Map(valid.map(n => n -> 0): _*)
    explicit.foreach { case (k, v) => counts(k) = v }

    val alreadyCut = internal.filter(i => counts(points(i).name) > 0).toSet
    val available  = internal.filterNot(alreadyCut.contains)

    val nAuto = math.min(auto, available.length)
    val chosen =
      if (nAuto == 0) Seq.empty[Int]
      else if (available.length <= PipelineStageConfig.bruteForceLimit)
        PipelineStageConfig.bestSubset(weights, nPoints, alreadyCut, available, nAuto)
      else
        PipelineStageConfig.greedySubset(weights, nPoints, alreadyCut, available, nAuto)

    chosen.foreach(i => counts(points(i).name) += 1)

    // surplus: every internal point is cut, stack the rest at the point that
    // ends the longest remaining segment
    var surplus = auto - nAuto
    if (surplus > 0) {
      val cuts = alreadyCut ++ chosen
      val segs = PipelineStageConfig.segments(weights, nPoints, cuts)
      // segment i ends at point `ends(i)`
      val ends = (cuts.toSeq.sorted :+ (nPoints - 1))
      val longest = segs.zipWithIndex.maxBy { case (w, i) => (w, -i) }._2
      counts(points(ends(longest)).name) += surplus
      surplus = 0
    }

    PipelineStages(counts.toMap, points)
  }

  /** Human readable form, e.g. `"auto:2"`, `"add:1 norm:1 round:1"`. */
  def getString: String = {
    val a = if (auto > 0) Seq(s"auto:$auto") else Seq.empty[String]
    val e = explicit.toSeq.sortBy { case (k, _) =>
      (if (k == "in") 0 else if (k == "out") 2 else 1, k)
    }.map { case (k, v) => s"$k:$v" }
    val s = a ++ e
    if (s.isEmpty) "none" else s.mkString(" ")
  }
}

/** Constructor shorthands and the automatic placement for
 *  [[rial.util.PipelineStageConfig]].
 */
object PipelineStageConfig {

  /** Number of internal points up to which the automatic placement is exact. */
  private[util] val bruteForceLimit = 16

  /** `n` registers, placed by the module.
   *
   * Also used implicitly, so a plain `Int` can be given wherever a
   * [[rial.util.PipelineStageConfig]] is expected.
   */
  implicit def default(n: Int): PipelineStageConfig = PipelineStageConfig(auto = n)

  /** Registers at explicitly named points, e.g.
   *  `specify("add" -> 1, "norm" -> 1, "round" -> 1)`.
   */
  def specify(kv: (String, Int)*): PipelineStageConfig = specify(kv.toMap)

  /** Registers at explicitly named points. */
  def specify(m: Map[String, Int]): PipelineStageConfig = PipelineStageConfig(explicit = m)

  /** `n` registers at the module output. */
  def atOut(n: Int): PipelineStageConfig = specify("out" -> n)

  /** `n` registers at the module input. */
  def atIn(n: Int): PipelineStageConfig = specify("in" -> n)

  /** No pipeline register: the module is purely combinational. */
  def none: PipelineStageConfig = PipelineStageConfig()

  /** All the names accepted for a module declaring `points`: `"in"` first. */
  private[util] def validNames(points: Seq[StagePoint]): Seq[String] =
    "in" +: points.map(_.name)

  /** Checks that a module's point list is well formed. */
  private[util] def checkPoints(points: Seq[StagePoint]): Unit = {
    if (points.isEmpty) throw new IllegalArgumentException(
      "a module's pipeline stage point list must not be empty; " +
      "its last element must be StagePoint(\"out\", w)")
    if (points.last.name != "out") throw new IllegalArgumentException(
      s"the last pipeline stage point must be named 'out', got '${points.last.name}'")
    val names = points.map(_.name)
    if (names.distinct.length != names.length) throw new IllegalArgumentException(
      s"duplicated pipeline stage point name(s) in ${names.mkString(", ")}")
    if (names.contains("in")) throw new IllegalArgumentException(
      "'in' is implicit and must not be declared as a pipeline stage point")
    if (names.init.contains("out")) throw new IllegalArgumentException(
      "'out' must only be used for the last pipeline stage point")
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
      if (cuts.contains(i)) { b += acc; acc = 0.0 }
    }
    b += acc // segment ending at the output
    b.result()
  }

  /** Segment lengths sorted in decreasing order: the cost to be minimised. */
  private def cost(weights: Seq[Double], nPoints: Int, cuts: Set[Int]): Seq[Double] =
    segments(weights, nPoints, cuts).sorted(Ordering[Double].reverse)

  /** Lexicographic comparison of two decreasing-sorted segment length lists. */
  private def lexLess(a: Seq[Double], b: Seq[Double]): Boolean = {
    val n = math.min(a.length, b.length)
    var i = 0
    while (i < n) {
      if (a(i) != b(i)) return a(i) < b(i)
      i += 1
    }
    a.length < b.length
  }

  /** Exact placement: the best `n`-element subset of `available`. */
  private def bestSubset(weights: Seq[Double], nPoints: Int, fixed: Set[Int],
    available: Seq[Int], n: Int): Seq[Int] = {
    var best: Seq[Int] = Seq.empty
    var bestCost: Seq[Double] = null
    for (sel <- available.combinations(n)) {
      val c = cost(weights, nPoints, fixed ++ sel)
      if ((bestCost == null) || lexLess(c, bestCost)) { best = sel; bestCost = c }
    }
    best
  }

  /** Greedy fallback for modules with a very large number of points. */
  private def greedySubset(weights: Seq[Double], nPoints: Int, fixed: Set[Int],
    available: Seq[Int], n: Int): Seq[Int] = {
    var cuts = fixed
    var rest = available
    val sel  = Seq.newBuilder[Int]
    for (_ <- 0 until n) {
      var best = rest.head
      var bestCost = cost(weights, nPoints, cuts + best)
      for (i <- rest.tail) {
        val c = cost(weights, nPoints, cuts + i)
        if (lexLess(c, bestCost)) { best = i; bestCost = c }
      }
      sel += best
      cuts = cuts + best
      rest = rest.filterNot(_ == best)
    }
    sel.result().sorted
  }
}
