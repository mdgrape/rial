package rial.util

import chisel3._
import chisel3.util._
import ScalaUtil._

object RialChiselUtil {

  // OBSOLETE-
  // dropU will be replaced by standard method tail.
  def dropU ( n: Int , x : UInt ) = x(x.getWidth-1-n,0)

  def takeU ( n: Int , x : UInt ) = x(x.getWidth-1, x.getWidth-n)

  def dropLSB ( n: Int , x : UInt ) : UInt = x(x.getWidth-1,n)
  def dropLSB ( n: Int , x : SInt ) : SInt = x(x.getWidth-1,n).asSInt

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

  def getMSB ( n: Int , x : UInt ) : UInt = { if (n<x.getWidth) x(x.getWidth-1,x.getWidth-n) else x }
  def getMSB ( n: Int , x : SInt ) : SInt = { if (n<x.getWidth) x(x.getWidth-1,x.getWidth-n).asSInt else x }

  def getLSB ( n: Int , x : UInt ) : UInt = x(n-1,0)
  def getLSB ( n: Int , x : SInt ) : SInt = x(n-1,0).asSInt

  def padLSB ( n: Int , x : UInt ) : UInt = { if (x.getWidth<n) Cat(x,0.U((n-x.getWidth).W)) else x }
  def padLSB ( n: Int , x : SInt ) : SInt = { if (x.getWidth<n) Cat(x,0.S((n-x.getWidth).W)).asSInt else x }

  def orRLSB ( n: Int , x : UInt ) : Bool = {
    if (n<=0) false.B
    else x(n.min(x.getWidth)-1,0).orR.asBool
  }

  def orRMSB ( n: Int , x : UInt ) : Bool = {
    if (n<=0) false.B
    else x(x.getWidth-1, (x.getWidth-n).max(0)).orR.asBool
  }

  // take a bit counting from MSB
  // 0 to take MSB
  def getMSB1 ( n: Int, x : UInt ) : Bool = { if (n<x.getWidth) x(x.getWidth-1-n).asBool else false.B }
  def getMSB1 ( n: Int, x : SInt ) : Bool = { if (n<x.getWidth) x(x.getWidth-1-n).asBool else false.B }

  def maskU (n: Int) : UInt = Fill(n, 1.U(1.W))

  // if en === 0, returns 0.U.asTypeOf(x). otherwise, returns x.
  def enableIf[T <: Data](en: UInt, x: T): T = {
    assert(en.getWidth == 1)
    (x.asUInt & Fill(x.asUInt.getWidth, en)).asTypeOf(x)
  }
}
