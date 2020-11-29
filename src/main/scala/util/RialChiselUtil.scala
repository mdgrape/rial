package rial.util

import chisel3._
import chisel3.util._

object RialChiselUtil {

  // OBSOLETE-
  // dropU will be replaced by standard method tail.
  def dropU ( n: Int , x : UInt ) = x(x.getWidth-1-n,0)

  def takeU ( n: Int , x : UInt ) = x(x.getWidth-1, x.getWidth-n)

  def dropULSB ( n: Int , x : UInt ) : UInt = x(x.getWidth-1,n)
  def dropSLSB ( n: Int , x : SInt ) : SInt = x(x.getWidth-1,n).asSInt

  // Separate x into slices with width w[i]
  // Return : Seq[UInt]
  // [0] : MSB side
  def getSlices ( x : UInt, w : Seq[Int] ) = {
    val wSum = w.scanRight(0)( (p,q) => p+q )
    val wTuple = wSum.init.map( p => p-1 ).zip(wSum.tail)
    //println(wTuple)
    wTuple.map(p => x(p._1, p._2))
  }

  def signExtend1( x : UInt, s : Int ) = {
    if (s == -1) Cat(1.U(1.W), x).asSInt
    else         Cat(0.U(1.W), x).asSInt
  }
  
}
