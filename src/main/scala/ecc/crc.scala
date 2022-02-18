//% @file crc.scala
//
// Cyclic Redunduncy Check
// part of RIAL - RIKEN Arithmetic and Logic library
// Copyright (C) Makoto Taiji RIKEN BDR 2022
//
package rial.ecc

import scala.language.reflectiveCalls
import scala.math._
import chisel3._
import chisel3.util._

object CRC {

  def getOrder(g: BigInt, n: Int) : Int = { if (g==0) n else getOrder(g>>1, n+1) }

  def calcCRC(x: BigInt, g: BigInt ) : BigInt = {
    val n = getOrder(g,-1) // Order of generation polynomial : width of g = n+1
    val w = getOrder(x,0)    // bit width of x
    ( w-1 to 0 by -1).foldLeft(x<<n)( (z, i) => 
      if (z.testBit(i+n)) { z ^ (g<<i) } else z )
  }

}


//
// Stateless CRC module 
//   g : coefficients of generator polynomial (with highest order)
//   w : Width of input
//
//   x  : input 
//   ri : reminder in
//   ro : reminder out
//
class CRCcore ( g : BigInt, w : Int ) extends Module {

  val n = CRC.getOrder(g, -1) // polynomial order

  if (n<1) {
    println(s"${this.getClass.getName}: Illegal polynomial order $n ($g)")
    sys.exit(1)
  }

  val io = IO(iodef = new Bundle {
    val x  = Input(UInt(w.W))
    val ri = Input(UInt(n.W))
    val ro = Output(UInt(n.W))
  } )

  val y = io.x ^ ( if (n>w) {
    // Reminder width > input width
    io.ri(n-1, n-w)
  } else {
    // Reminder width <= input width
    io.ri << (w-n)
  } )

  // Calculate CRC of 2^(i+n)
  val r = (0 to w-1).map( i => CRC.calcCRC(BigInt(1)<<i, g) )

  // bitMask(p) : if bit q is 1, y(q) is involved in xor list in r(p)
  val bitMask = (0 to n-1).map( p =>
    (0 to w-1).foldLeft(BigInt(0))( (z, q) => 
      if (r(q).testBit(p)) { z + (BigInt(1)<<q) } else z ) )

  val rv = VecInit(
    (0 to n-1).map( i => (y & bitMask(i).U(w.W)).xorR ) )

  io.ro := rv.asUInt
}

// Parity
object CRCcoreGen_1_32 extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new CRCcore(3, 32)) ) )
}

// CRC16_IBM
object CRCcoreGen_16_16 extends App {
  (new chisel3.stage.ChiselStage).execute(args,
    Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new CRCcore(0x18005, 16)) ) )
}
