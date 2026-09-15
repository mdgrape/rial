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
// The fraction must keep at least two bits below the address: a table that
// has none is read as a plain lookup, of order zero, which asks for no
// position at all. How many bits the position carries beyond those two is a
// property of the polynomial, not of the fraction, so a fraction that has
// fewer bits left than the position asks for is padded with zeros on the
// right: the position is then the exact one, and the polynomial reads it at
// whatever precision it was built for.
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
 *
 *  Only a table read by a polynomial has a position: the fraction must keep
 *  at least two bits below the address, and the position must be at least two
 *  bits wide. A table whose address takes the whole fraction is read as a
 *  plain lookup, of order zero, and asks for no position.
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
   *  is the signed position in [-1, 1). A fraction that has fewer bits left
   *  than `dxW` is padded with zeros on the right.
   *
   *  The fraction must keep at least two bits below the address, and the
   *  position must be at least two bits wide; read a table that has neither
   *  as a plain lookup, of order zero, which needs no position.
   *
   *  @param frac  the fraction the table is read at
   *  @param addrW width of the address
   *  @param dxW   width of the position
   */
  def dx(frac: UInt, addrW: Int, dxW: Int): UInt = {
    val fracW = frac.getWidth
    require(dxW >= 2,
      f"a position of ${dxW}%d bits is none a polynomial can read: read the " +
      f"table as a plain lookup, of order zero, instead")
    require(fracW - addrW >= 2,
      f"an address of ${addrW}%d bits leaves a fraction of ${fracW}%d bits " +
      f"no room for a position: at least two bits of it must stay below the " +
      f"address for a polynomial to read; read the table as a plain lookup, " +
      f"of order zero, instead")

    val belowW = fracW - addrW
    val takenW = min(belowW, dxW)
    val taken  = Cat(~frac(belowW - 1), frac(belowW - 2, belowW - takenW))

    taken << (dxW - takenW)
  }

  /** The step between two positions the fraction can form.
   *
   *  A fraction that has fewer bits left than the position asks for is padded
   *  with zeros, so the positions it reaches are that far apart; a fraction
   *  that has enough reaches every one of them.
   *
   *  It describes the positions [[dx]] gives, so it asks of its widths what
   *  [[dx]] does: at least two bits of the fraction below the address, and a
   *  position at least two bits wide.
   *
   *  @param fracW width of the fraction
   *  @param addrW width of the address
   *  @param dxW   width of the position
   */
  def dxStep(fracW: Int, addrW: Int, dxW: Int): BigInt = {
    require(dxW >= 2,
      f"a position of ${dxW}%d bits is none a polynomial can read: read the " +
      f"table as a plain lookup, of order zero, instead")
    require(fracW - addrW >= 2,
      f"an address of ${addrW}%d bits leaves a fraction of ${fracW}%d bits " +
      f"no room for a position: at least two bits of it must stay below the " +
      f"address for a polynomial to read; read the table as a plain lookup, " +
      f"of order zero, instead")
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
   *  it itself.
   *
   *  The fraction must keep at least two bits below the address, and the
   *  position must be at least two bits wide; read a table that has neither
   *  as a plain lookup, of order zero, which needs no position.
   *
   *  @param frac  the fraction the table is read at
   *  @param fracW width of the fraction
   *  @param addrW width of the address
   *  @param dxW   width of the position
   */
  def dx(frac: SafeLong, fracW: Int, addrW: Int, dxW: Int): BigInt = {
    require(dxW >= 2,
      f"a position of ${dxW}%d bits is none a polynomial can read: read the " +
      f"table as a plain lookup, of order zero, instead")
    require(fracW - addrW >= 2,
      f"an address of ${addrW}%d bits leaves a fraction of ${fracW}%d bits " +
      f"no room for a position: at least two bits of it must stay below the " +
      f"address for a polynomial to read; read the table as a plain lookup, " +
      f"of order zero, instead")

    val belowW = fracW - addrW
    val takenW = min(belowW, dxW)

    slice(belowW - takenW, takenW, frac).toBigInt << (dxW - takenW)
  }
}
