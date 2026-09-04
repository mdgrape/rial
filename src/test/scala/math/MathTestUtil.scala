package rial.math

import rial.util.PipelineStageConfig

/** Pipeline configurations shared by the [[rial.math.MathFunctions]] tests. */
object MathTestUtil {

  /** One pipeline register at every boundary the module declares for `cfg`.
   *
   * The `"in"` pseudo-point of a group sits where the boundary preceding the
   * group already is, so it is left uncut.
   *
   * @param cfg configuration of the module the boundaries are taken from
   */
  def allPoints(cfg: MathFuncConfig): PipelineStageConfig = {
    val points = MathFunctions.stagePoints(cfg)
    val names  = PipelineStageConfig.none.resolve(points).counts.keys.filterNot(n => {
      n.endsWith(".in")
    })
    PipelineStageConfig.specify(names.map(n => n -> 1).toMap)
  }

  /** A short name for `stage`, for the name of a test case.
   *
   * [[allPoints]] names every boundary of the module, which would make a test
   * name (and so a test log file name) longer than the file system accepts.
   *
   * @param cfg   configuration of the module `stage` is for
   * @param stage the configuration to name
   */
  def label(cfg: MathFuncConfig, stage: PipelineStageConfig): String = {
    if (stage == allPoints(cfg)) {
      "allPoints"
    } else {
      stage.getString
    }
  }
}
