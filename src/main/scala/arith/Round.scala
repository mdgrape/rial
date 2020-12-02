package rial.util

import chisel3._
import chisel3.util._
import scala.math._
import Numeric._

sealed trait RoundSpec

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

object Rounding {

  // BigInt
  private def tail( x: BigInt, pos : Int ) : BigInt = {
    val one = BigInt(1)
    ((one << pos)-one) & x
  }

  private def bit( x: BigInt, pos : Int ) : Int = ((x>>pos) & 1).toInt

  def round(x: BigInt, pos : Int) : BigInt = {
    if (pos<1) x else ((x >> pos) + bit(x, pos-1))
  }

  def roundToEven(x: BigInt, pos : Int) : BigInt = {
    if (pos<1) return x
    val s = if (tail(x,pos-1) != 0) 1 else 0
    val l = bit(x, pos)
    val r = bit(x, pos-1)
    val inc = r & (l | s)
    (x >> pos) + inc
  }

  def truncate(x: BigInt, pos : Int) : BigInt = if (pos<1) x else (x >> pos)

  def truncateToZero(x: BigInt, pos : Int) : BigInt = 
    if (pos<1) x else if (x<0) (-((-x) >> pos)) else (x>>pos)

  def ceil(x: BigInt, pos : Int) : BigInt = {
    if (pos<1) x
    else if (tail(x,pos) == 0) (x>>pos)
    else ((x >> pos) + 1)
  }

  def ceilToInfinite(x: BigInt, pos : Int) : BigInt = {
    if (pos<1) x
    else if (tail(x,pos) == 0) (x>>pos)
    else if (x<0) (-(((-x) >> pos) + 1))
    else ((x >> pos) + 1)
  }

  // Long
  private def tail( x: Long, pos : Int ) : Long = {
    if (pos>=64) x
    else if (pos==63) (x & 0x7FFFFFFFFFFFFFFFL)
    else if (pos>0)   ( ((1L << pos)-1L) & x )
    else 0L
  }

  private def bit( x: Long, pos : Int ) : Int =
    ((x>>pos) & 1L).toInt

  def round(x: Long, pos : Int) : Long = {
    if (pos<1) x else ((x >> pos) + bit(x, pos-1))
  }

}

