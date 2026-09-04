package rial.util

import chisel3._
import chisel3.util._
import ScalaUtil._

/** One cut of a pipelined datapath, with the semantics of `chisel3.util.Pipe`:
 * one valid flag per cut, reset to false, travels with the data and enables
 * the data registers of its own stage, so a stage that carries nothing does
 * not clock its registers.
 *
 * A module declares its boundaries with [[rial.util.StagePoint]], places one
 * cut per boundary and pipes every signal that crosses the boundary through
 * it. The valid registers are named after the boundary, as the data ones are
 * named after the signals:
 *
 * {{{
 * private val cutClz = new SharedPipe("clz", validIn, stages.at("clz"))
 *
 * val rndClz   = cutClz(rnd)       // enabled by the valid of its own stage
 * val validClz = cutClz.valid
 * }}}
 *
 * @constructor create the valid registers of one cut.
 * @param name name of the boundary the cut sits at
 * @param v    the valid flag entering the cut
 * @param n    number of pipeline registers of the cut
 *
 * @see [[rial.util.PipelineStageConfig]]
 */
class SharedPipe(name: String, v: Bool, n: Int) {

  /** the valid entering each stage of the cut, then the one leaving it. */
  private val valids = Seq.iterate(v, n+1)(prev => RegNext(prev, false.B))

  // name the registers after the boundary, as the data ones are
  valids.tail.zipWithIndex.foreach { case (b, i) =>
    val suffix = if (n == 1) {
      ""
    } else {
      i.toString
    }
    b.suggestName(s"valid${name.capitalize}${suffix}")
  }

  /** the valid flag leaving the cut. */
  val valid: Bool = valids.last

  /** pipe `x` through the cut, enabled by the valid entering each stage. */
  def apply[T <: Data](x: T): T = {
    valids.init.foldLeft(x)((d, en) => RegEnable(d, en))
  }
}

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
