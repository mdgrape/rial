package rial.util

import scala.language.implicitConversions

/** Enumerator for pipeline stage strategy. */
sealed trait PipelineStageStrategy 

/** List of available pipeline stage strategies. */
object PipelineStageStrategy {
  case object none    extends PipelineStageStrategy
  case object default extends PipelineStageStrategy
  case object atOut   extends PipelineStageStrategy
  case object atIn    extends PipelineStageStrategy
  case object specify extends PipelineStageStrategy
}

/** Pipeline stage configuration
 *
 * For simplification and customizability, we basically use a single
 * `ShiftRegister` when specifying latency of a module, leaving the specific
 * register location to optimizations such as register balancing / retiming.
 *
 * Normally, it is constructed by [[rial.util.PipelineStageConfig.atOut]].
 */
class PipelineStageConfig (
  strategy : PipelineStageStrategy,
  position : Seq[Int] ) {

  def this() = this( PipelineStageStrategy.none, Seq[Int]() )

  def total = position.sum

  private def posToStr() = position.foldLeft[String](""){
    (s,x) => s + ( if (x==0) '-'
    else if ((x>0) && (x<=9)) x.toChar
    else if ((x>=10)&&(x<36)) (x+'A'.toInt-10).toChar
    else {
      println(f"Illegal value $x%d for pipeline register depth at one place; set to 0")
      '-'
    }
    )
  }

  // get explanation
  def getString = {
    val p0 = if (position.isDefinedAt(0)) position(0) else 0
    strategy match {
      case PipelineStageStrategy.none => "none"
      case PipelineStageStrategy.default => f"$p0%d stages at default (recommended) positions"
      case PipelineStageStrategy.atOut   => f"$p0%d stages at output"
      case PipelineStageStrategy.atIn    => f"$p0%d stages at input"
      case PipelineStageStrategy.specify => "Stage configuration: "+posToStr()
      case _ => "Illegal strategy"
    }
  }

}

/** Constructor shorthands for PipelineStageConfig */
object PipelineStageConfig {
  implicit def default(n: Int ) : PipelineStageConfig = {
    new PipelineStageConfig( PipelineStageStrategy.default, Seq[Int](n) )
  }

  private def char2int( c : Char ) : Int = {
    if (c.isDigit) { c.asDigit }
    else if (c.isLetter) { c.toUpper.toInt - 'A'.toInt + 10 }
    else if ( (c=='-')||(c==' ')) { 0 }
    else { 1 }
  }

  implicit def specify(s: String) : PipelineStageConfig = {
    new PipelineStageConfig( PipelineStageStrategy.specify,
      s.map(c => char2int(c)))
  }

  def atOut(n: Int)  : PipelineStageConfig = {
    new PipelineStageConfig( PipelineStageStrategy.atOut, Seq[Int](n) )
  }

  def atIn(n: Int)  : PipelineStageConfig = {
    new PipelineStageConfig( PipelineStageStrategy.atIn, Seq[Int](n) )
  }

  def none : PipelineStageConfig = {
    new PipelineStageConfig
  }
}
