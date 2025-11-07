package rial.arith

import scala.math._
import spire.math.SafeLong
import spire.implicits._

/** Enumerator to select Rounding Algorithm.
 *
 * `object` defines available rounding algorithms.
 */
sealed trait RoundSpec

/** Enumerator to select Rounding Algorithm.
 */
object RoundSpec {
  case object round          extends RoundSpec
  case object roundToEven    extends RoundSpec
  case object truncate       extends RoundSpec
  case object truncateToZero extends RoundSpec
  case object ceil           extends RoundSpec
  case object ceilToInfinite extends RoundSpec
}

// It seems to be difficult to use the generic type constraint
// for Int, Long, BigInt because there are no child relationships.
// Numeric[T]/Integral[T] does not have shift and logical operators.
// We tried :
// private def mask[T]( pos : Int ) (implicit numeric: Numeric[T]) : T = { ...
//  private def bit[T]( x : T, pos : Int )(implicit ev: T <:< { def >>(x:Int) : T }) : Int = {
// the latter works but too complicated.
//
// So, it seems to be better to define for Int/Long/BigInt separately...

/** Rounding functions for software emulator
 */
object Rounding {

  ////////////////////////////////////////////////////////////////////////
  // SafeLong
  private def tail( pos : Int, x: SafeLong ) : SafeLong = {
    val one = SafeLong(1)
    ((one << pos)-one) & x
  }

  private def bit( pos : Int, x: SafeLong ) : Int = ((x>>pos) & 1).toInt

  def round( pos : Int, x: SafeLong ) : SafeLong = {
    if (pos<1) x else ((x >> pos) + bit(pos-1,x))
  }

  def roundToEven( pos : Int, x: SafeLong ) : SafeLong = {
    if (pos<1) return x
    val s = if (tail(pos-1,x) != 0) 1 else 0
    val l = bit(pos,x)
    val r = bit(pos-1,x)
    val inc = r & (l | s)
    (x >> pos) + inc
  }

  def truncate( pos : Int, x: SafeLong ) : SafeLong = if (pos<1) x else (x >> pos)

  def truncateToZero( pos : Int, x: SafeLong ) : SafeLong = 
    if (pos<1) x else if (x<0) (-((-x) >> pos)) else (x>>pos)

  def ceil( pos : Int, x: SafeLong ) : SafeLong = {
    if (pos<1) x
    else if (tail(pos,x) == 0) (x>>pos)
    else ((x >> pos) + 1)
  }

  def ceilToInfinite( pos : Int, x: SafeLong ) : SafeLong = {
    if (pos<1) x
    else if (tail(pos, x) == 0) (x>>pos)
    else if (x<0) (-(((-x) >> pos) + 1))
    else ((x >> pos) + 1)
  }

  ////////////////////////////////////////////////////////////////////////
  // BigInt
  private def tail( pos : Int, x: BigInt ) : BigInt = {
    val one = BigInt(1)
    ((one << pos)-one) & x
  }

  private def bit( pos : Int, x: BigInt ) : Int = ((x>>pos) & 1).toInt

  def round( pos : Int, x: BigInt ) : BigInt = {
    if (pos<1) x else ((x >> pos) + bit(pos-1,x))
  }

  def roundToEven( pos : Int, x: BigInt ) : BigInt = {
    if (pos<1) return x
    val s = if (tail(pos-1,x) != 0) 1 else 0
    val l = bit(pos,x)
    val r = bit(pos-1,x)
    val inc = r & (l | s)
    (x >> pos) + inc
  }

  def truncate( pos : Int, x: BigInt ) : BigInt = if (pos<1) x else (x >> pos)

  def truncateToZero( pos : Int, x: BigInt ) : BigInt = 
    if (pos<1) x else if (x<0) (-((-x) >> pos)) else (x>>pos)

  def ceil( pos : Int, x: BigInt ) : BigInt = {
    if (pos<1) x
    else if (tail(pos,x) == 0) (x>>pos)
    else ((x >> pos) + 1)
  }

  def ceilToInfinite( pos : Int, x: BigInt ) : BigInt = {
    if (pos<1) x
    else if (tail(pos, x) == 0) (x>>pos)
    else if (x<0) (-(((-x) >> pos) + 1))
    else ((x >> pos) + 1)
  }

  ////////////////////////////////////////////////////////////////////////
  // Long
  private def tail( pos : Int, x: Long ) : Long = {
    if (pos>=64) x
    else if (pos==63) (x & 0x7FFFFFFFFFFFFFFFL)
    else if (pos>0)   ( ((1L << pos)-1L) & x )
    else 0L
  }

  private def bit( pos : Int, x: Long ) : Int =
    ((x>>pos) & 1L).toInt

  def round( pos : Int, x: Long ) : Long = {
    if (pos<1) x else ((x >> pos) + bit(pos-1,x))
  }

  def roundToEven( pos : Int, x: Long ) : Long = {
    if (pos<1) return x
    val s = if (tail(pos-1,x) != 0) 1 else 0
    val l = bit(pos,x)
    val r = bit(pos-1,x)
    val inc = r & (l | s)
    (x >> pos) + inc
  }

  def truncate( pos : Int, x: Long ) : Long = if (pos<1) x else (x >> pos)

  def truncateToZero( pos : Int, x: Long ) : Long = 
    if (pos<1) x else if (x<0L) (-((-x) >> pos)) else (x>>pos)

  def ceil( pos : Int, x: Long ) : Long = {
    if (pos<1) x
    else if (tail(pos,x) == 0L) (x>>pos)
    else ((x >> pos) + 1L)
  }

  def ceilToInfinite( pos : Int, x: Long ) : Long = {
    if (pos<1) x
    else if (tail(pos, x) == 0L) (x>>pos)
    else if (x<0) (-(((-x) >> pos) + 1L))
    else ((x >> pos) + 1L)
  }

  ////////////////////////////////////////////////////////////////////////
  // Int
  private def tail( pos : Int, x: Int ) : Int = {
    ((1 << pos)-1) & x
  }

  private def bit( pos : Int, x: Int ) : Int = ((x>>pos) & 1).toInt

  def round( pos : Int, x: Int ) : Int = {
    if (pos<1) x else ((x >> pos) + bit(pos-1,x))
  }

  def roundToEven( pos : Int, x: Int ) : Int = {
    if (pos<1) return x
    val s = if (tail(pos-1,x) != 0) 1 else 0
    val l = bit(pos,x)
    val r = bit(pos-1,x)
    val inc = r & (l | s)
    (x >> pos) + inc
  }

  def truncate( pos : Int, x: Int ) : Int = if (pos<1) x else (x >> pos)

  def truncateToZero( pos : Int, x: Int ) : Int = 
    if (pos<1) x else if (x<0) (-((-x) >> pos)) else (x>>pos)

  def ceil( pos : Int, x: Int ) : Int = {
    if (pos<1) x
    else if (tail(pos,x) == 0) (x>>pos)
    else ((x >> pos) + 1)
  }

  def ceilToInfinite( pos : Int, x: Int ) : Int = {
    if (pos<1) x
    else if (tail(pos, x) == 0) (x>>pos)
    else if (x<0) (-(((-x) >> pos) + 1))
    else ((x >> pos) + 1)
  }

  ////////////////////////////////////////////////////////////////////////

  def roundBySpec( spec: RoundSpec, pos : Int, x : SafeLong ) : SafeLong = {
    spec match {
      case RoundSpec.round          => round(pos, x)
      case RoundSpec.roundToEven    => roundToEven(pos, x)
      case RoundSpec.truncate       => truncate(pos, x)
      case RoundSpec.truncateToZero => truncateToZero(pos, x)
      case RoundSpec.ceil           => ceil(pos, x)
      case RoundSpec.ceilToInfinite => ceilToInfinite(pos, x)
      case _ => x
    }
  }

  def roundBySpec( spec: RoundSpec, pos : Int, x : BigInt ) : BigInt = {
    spec match {
      case RoundSpec.round          => round(pos, x)
      case RoundSpec.roundToEven    => roundToEven(pos, x)
      case RoundSpec.truncate       => truncate(pos, x)
      case RoundSpec.truncateToZero => truncateToZero(pos, x)
      case RoundSpec.ceil           => ceil(pos, x)
      case RoundSpec.ceilToInfinite => ceilToInfinite(pos, x)
      case _ => x
    }
  }

  def roundBySpec( spec: RoundSpec, pos : Int, x : Long ) : Long = {
    spec match {
      case RoundSpec.round          => round(pos, x)
      case RoundSpec.roundToEven    => roundToEven(pos, x)
      case RoundSpec.truncate       => truncate(pos, x)
      case RoundSpec.truncateToZero => truncateToZero(pos, x)
      case RoundSpec.ceil           => ceil(pos, x)
      case RoundSpec.ceilToInfinite => ceilToInfinite(pos, x)
      case _ => x
    }
  }

  def roundBySpec( spec: RoundSpec, pos : Int, x : Int ) : Int = {
    spec match {
      case RoundSpec.round          => round(pos, x)
      case RoundSpec.roundToEven    => roundToEven(pos, x)
      case RoundSpec.truncate       => truncate(pos, x)
      case RoundSpec.truncateToZero => truncateToZero(pos, x)
      case RoundSpec.ceil           => ceil(pos, x)
      case RoundSpec.ceilToInfinite => ceilToInfinite(pos, x)
      case _ => x
    }
  }


}

