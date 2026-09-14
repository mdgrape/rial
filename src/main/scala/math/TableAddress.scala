//% @file TableAddress.scala
//
// The table address and the position inside the sub-interval it selects.
//
// Every function that reads an interval table forms both of them from one
// fraction: the top bits of the fraction name the sub-interval, and the bits
// below them say where inside it the argument lies. The position is the
// domain of the interval polynomial, [-1, 1), so it is signed, and the top
// bit of the fraction carries its sign inverted.
//
// How many bits the position carries is a property of the polynomial, not of
// the fraction, so a fraction that has fewer bits left than the position asks
// for is padded with zeros on the right: the position is then the exact one,
// and the polynomial reads it at whatever precision it was built for.
//
// Copyright (C) Toru Niina RIKEN BDR 2022
//
package rial.math

import scala.math.min

import chisel3._
import chisel3.util._

import spire.math.SafeLong
import spire.implicits._

import rial.util.ScalaUtil._

/** The address a table is read at and the position inside the sub-interval
 *  that address selects, formed from one fraction.
 *
 *  The hardware forms of both take the fraction as a `UInt` and the software
 *  ones take it as a `SafeLong` together with its width; the values they give
 *  are the same. The position the software gives is unsigned, as
 *  [[rial.table.PolynomialEval.eval]] reads it.
 */
object TableAddress {

  /** The table address: the top `addrW` bits of the fraction.
   *
   *  @param frac  the fraction the table is read at
   *  @param addrW width of the address
   */
  def adr(frac: UInt, addrW: Int): UInt = {
    val fracW = frac.getWidth
    require(fracW >= addrW,
      f"a fraction of ${fracW}%d bits cannot hold an address of " +
      f"${addrW}%d bits")
    frac(fracW - 1, fracW - addrW)
  }

  /** The position inside the sub-interval: the `dxW` bits of the fraction
   *  below its top `addrW`, with the top of them inverted so that the result
   *  is the signed position in [-1, 1).
   *
   *  @param frac  the fraction the table is read at
   *  @param addrW width of the address
   *  @param dxW   width of the position
   */
  def dx(frac: UInt, addrW: Int, dxW: Int): UInt = {
    val fracW = frac.getWidth
    require(fracW >= addrW,
      f"a fraction of ${fracW}%d bits cannot hold an address of " +
      f"${addrW}%d bits")
    require(dxW >= 1,
      f"a position of ${dxW}%d bits is none the polynomial can read")

    val belowW = fracW - addrW
    val takenW = min(belowW, dxW)
    val padW   = dxW - takenW

    if (takenW == 0) {
      // the address used every bit the fraction has, so the argument sits at
      // the left end of its sub-interval, which is -1 in the polynomial's
      // domain: the inverted top bit alone spells it
      (BigInt(1) << (dxW - 1)).U(dxW.W)
    } else {
      val taken = if (takenW == 1) {
        ~frac(belowW - 1)
      } else {
        Cat(~frac(belowW - 1), frac(belowW - 2, belowW - takenW))
      }
      if (padW == 0) {
        taken
      } else {
        Cat(taken, 0.U(padW.W))
      }
    }
  }

  /** The step between two positions the fraction can form.
   *
   *  A fraction that has fewer bits left than the position asks for is padded
   *  with zeros, so the positions it reaches are that far apart; a fraction
   *  that has enough reaches every one of them.
   *
   *  @param fracW width of the fraction
   *  @param addrW width of the address
   *  @param dxW   width of the position
   */
  def dxStep(fracW: Int, addrW: Int, dxW: Int): BigInt = {
    require(fracW >= addrW,
      f"a fraction of ${fracW}%d bits cannot hold an address of " +
      f"${addrW}%d bits")
    BigInt(1) << (dxW - min(fracW - addrW, dxW))
  }

  /** The table address, for a software model.
   *
   *  @param frac  the fraction the table is read at
   *  @param fracW width of the fraction
   *  @param addrW width of the address
   */
  def adr(frac: SafeLong, fracW: Int, addrW: Int): Int = {
    require(fracW >= addrW,
      f"a fraction of ${fracW}%d bits cannot hold an address of " +
      f"${addrW}%d bits")
    slice(fracW - addrW, addrW, frac).toInt
  }

  /** The position inside the sub-interval, for a software model. It is the
   *  unsigned form: [[rial.table.PolynomialEval.eval]] takes the sign out of
   *  it itself, and reads nothing at all when the polynomial is of order
   *  zero, where the position is zero bits wide.
   *
   *  @param frac  the fraction the table is read at
   *  @param fracW width of the fraction
   *  @param addrW width of the address
   *  @param dxW   width of the position
   */
  def dx(frac: SafeLong, fracW: Int, addrW: Int, dxW: Int): BigInt = {
    require(fracW >= addrW,
      f"a fraction of ${fracW}%d bits cannot hold an address of " +
      f"${addrW}%d bits")

    val belowW = fracW - addrW
    val takenW = min(belowW, dxW)

    if (dxW == 0) {
      BigInt(0)
    } else if (takenW == 0) {
      BigInt(1) << (dxW - 1)
    } else {
      slice(belowW - takenW, takenW, frac).toBigInt << (dxW - takenW)
    }
  }
}
