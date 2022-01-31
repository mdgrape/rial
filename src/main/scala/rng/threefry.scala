//% @file threefry.scala
//
// Threefry Random Number Generation
//   based on J. Salomon et al.
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.rng

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._

// Fully pipelined implementation, no loop
// r : Number of rotations
// rotStage : Insert pipeline stage every rotStage rotation
class Threefry4_32( r: Int = 20, rotStage: Int = 0 ) extends Module {

  def getParam() = { (r, rotStage) }

  def rotL( x : UInt, n : Int ) = {
    Cat( x(x.getWidth-n-1,0), x(x.getWidth-1,x.getWidth-n) )
  }

  //val Threefish_C240 = 0x1BD11BDAA9FC1A22
  val C240 = 0x1BD11BDA.U

  // Rotation constant
  val R0 = Array[Int]( 10, 11, 13, 23, 6, 17, 25, 18 )
  val R1 = Array[Int]( 26, 21, 27, 5, 20, 11, 10, 20 )

  val io = IO(new Bundle {
    val key   = Input(Vec(4, UInt(32.W)))
    val count = Input(Vec(4, UInt(32.W)))
    val rand  = Output(Vec(4, UInt(32.W)))
  })

  val keys = Wire(Vec(5, UInt(32.W)))
  val X = Wire(Vec(r+1, Vec(4, UInt(32.W))))
  val Xnext = Wire(Vec(r, Vec(4, UInt(32.W))))

  for(i <- 0 to 3) {
    keys(i) := io.key(i)
    X(0)(i) := io.count(i) + io.key(i)
  }

  keys(4) := (0 to 3).toSeq.foldLeft(C240)( (z,n) => z ^ io.key(n) )

  for(j <- 0 to r-1) {
    val keyR = (j/4) + 1
    if (j%2 == 0) {
      Xnext(j)(0) := X(j)(0) + X(j)(1)
      Xnext(j)(1) := rotL( X(j)(1), R0(j%8) ) ^ Xnext(j)(0)
      Xnext(j)(2) := X(j)(2) + X(j)(3)
      Xnext(j)(3) := rotL( X(j)(3), R1(j%8) ) ^ Xnext(j)(2)
    } else {
      Xnext(j)(0) := X(j)(0) + X(j)(3)
      Xnext(j)(3) := rotL( X(j)(3), R0(j%8) ) ^ Xnext(j)(0)
      Xnext(j)(2) := X(j)(2) + X(j)(1)
      Xnext(j)(1) := rotL( X(j)(1), R1(j%8) ) ^ Xnext(j)(2)
    }
    if ( (rotStage == 0) || ( ((j+1)%rotStage) != 0 ) ) { // No staging
      if (j%4 == 3) {
        X(j+1)(0) := Xnext(j)(0) + keys((keyR+0)%5)
        X(j+1)(3) := Xnext(j)(3) + keys((keyR+3)%5) + keyR.U
        X(j+1)(2) := Xnext(j)(2) + keys((keyR+2)%5)
        X(j+1)(1) := Xnext(j)(1) + keys((keyR+1)%5)
      } else {
        X(j+1) := Xnext(j)
      }
    } else {
      if (j%4 == 3) {
        X(j+1)(0) := RegNext(Xnext(j)(0) + keys((keyR+0)%5))
        X(j+1)(3) := RegNext(Xnext(j)(3) + keys((keyR+3)%5) + keyR.U)
        X(j+1)(2) := RegNext(Xnext(j)(2) + keys((keyR+2)%5))
        X(j+1)(1) := RegNext(Xnext(j)(1) + keys((keyR+1)%5))
      } else {
        X(j+1) := RegNext(Xnext(j))
      }
    }
  }

  //for(j <- 0 to r) {
  //printf("%d %x %x %x %x %x\n", j.U, X(j)(0), X(j)(1), X(j)(2), X(j)(3), rotL( X(j)(1), R0(j%8) ))
  //}
  io.rand := X(r)

}

object Threefry4_32Driver extends App {
  (new chisel3.stage.ChiselStage).execute(args,
      Seq(chisel3.stage.ChiselGeneratorAnnotation(() => new Threefry4_32(20,2)))
    )
}
