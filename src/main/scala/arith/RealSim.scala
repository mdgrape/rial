package rial.arith

import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._
import rial.arith._
import rial.arith.RealGeneric._
import rial.arith.Rounding._
import rial.util.ScalaUtil._

import java.lang.Math.scalb
import java.nio.LongBuffer

class RealSpec (
  val exW    : Int,
  val exBias : Int,
  val manW   : Int,
  val disableSign : Boolean = false,
  val disableNaN  : Boolean = false,
  val disableSubnormal : Boolean = true ) {
}

object RealSpec {
  val Float64Spec = new RealSpec(11, 0x3FF, 52, false, false, true)
  val Float32Spec = new RealSpec(8, 0x7F, 23, false, false, true)
}

//
// Again we do not use type generics for shift operators here...
//
// class RealGeneric
//   this uses SafeLong for mantissa representation
//   mantissa includes leading 1
//
class RealGeneric ( val spec : RealSpec, val sgn: Int, val ex : Int, val man: SafeLong  ) {

  def fracIsZero : Boolean = { (man & maskSL(spec.manW)) == 0 }

  def isNaN : Boolean = { (ex == maskI(spec.exW)) && !fracIsZero && !spec.disableNaN }
  def isInfinite : Boolean = { (ex == maskI(spec.exW)) && ( fracIsZero || spec.disableNaN ) }
  def isZero : Boolean = { (ex == 0) && (spec.disableSubnormal || fracIsZero) }
  def isSubnormal : Boolean = { (ex == 0) && (man != 0) }

  // Multiplication
  //
  // Need to check IEEE 754...
  // Zero * (Any value including NaN/Inf) => 0
  // NaN  * (Any value != 0) => NaN
  // Inf  * (Any value != 0 and != NaN) => Inf
  def multiply( resSpec : RealSpec, roundSpec : RoundSpec, that : RealGeneric ) : RealGeneric = {
    val resSgn = this.sgn ^ that.sgn

    if (this.isZero || that.isZero)    return zero(resSpec)
    else if (this.isNaN || that.isNaN) {
      if (resSpec.disableNaN) return inf(resSpec, resSgn)
      else return nan(resSpec)
    } else if (this.isInfinite || that.isInfinite) return inf(resSpec, resSgn)

    // Normal values
    val ex0=(ex-spec.exBias)+(that.ex-that.spec.exBias)
    val prod= man * that.man
    val bp = spec.manW+that.spec.manW
    val moreThan2 = bit(bp+1,prod)
    val roundbits = bp-resSpec.manW+moreThan2
    val prodRound = roundBySpec(roundSpec, roundbits, prod)
    val moreThan2AfterRound = bit(resSpec.manW+1, prodRound)
    val resMan = prodRound >> moreThan2AfterRound
    val exNobias = ex0+(moreThan2|moreThan2AfterRound)
    val resEx = exNobias + resSpec.exBias
    if (resEx <= 0) { // result is zero : currently subnormal not supported
      zero(resSpec)
    } else if (resEx >= maskI(resSpec.exW)) { // result is inf
      inf(resSpec, resSgn)
    } else {
      new RealGeneric(resSpec, resSgn, resEx, resMan)
    }
  }

  // Addition
  //
  // NaN  +/- (Any value) => NaN
  // Inf  +/- (Any value except NaN) => Inf
  //
  def add( resSpec : RealSpec, roundSpec : RoundSpec, that : RealGeneric ) : RealGeneric = {

    if (this.isNaN || that.isNaN) {
      if (resSpec.disableNaN) {
        return ( if (this.isNaN) inf(resSpec, sgn) else inf(resSpec, that.sgn) )
      } else return nan(resSpec)
    }
    if (this.isInfinite) return inf(resSpec, sgn)
    if (that.isInfinite) return inf(resSpec, that.sgn)

    // Normal values
    val thisEx=ex-this.spec.exBias
    val thatEx=that.ex-that.spec.exBias
    val (x,y,exMax) = if (thisEx<thatEx) (that, this, thatEx) else (this, that, thisEx)
    val diffEx = abs(thisEx-thatEx)
    val xshift = x.man << diffEx
    val sum =  if (this.sgn!=that.sgn) xshift-y.man else xshift+y.man
    // In case of diffEx==0, this can be negative
    val (sumPos, resSgn) = if (sum<0) (-sum, x.sgn^1) else (sum,x.sgn)
    if (sumPos==0) return zero(resSpec)
    // Normalization
    val leading1 = sumPos.bitLength-1
    //val ml=sumPos.toLong; println(f"$ml%x $leading1%d $diffEx%d")
    val roundbits = leading1-resSpec.manW
    val (resMan, exInc)  =
      if (roundbits>0) {
        val sumRound = roundBySpec(roundSpec, roundbits, sumPos)
        val moreThan2AfterRound = bit(resSpec.manW+1, sumRound)
        (sumRound >> moreThan2AfterRound, moreThan2AfterRound)
      } else {
        ( sumPos << (-roundbits), 0)
      }
    // no exponent change if leading1 = x.manW+diffEx
    val resEx = exMax+leading1-x.spec.manW-diffEx+exInc+resSpec.exBias
    if (resEx <= 0) { // result is zero : currently subnormal not supported
      zero(resSpec)
    } else if (resEx >= maskI(resSpec.exW)) { // result is inf
      inf(resSpec, resSgn)
    } else {
      new RealGeneric(resSpec, resSgn, resEx, resMan)
    }
  }

  def packToSafeLong : SafeLong = {
    (man & maskSL(spec.manW))+ ( SafeLong( ex + (sgn<<spec.exW) ) << spec.manW )
  }

  def toDouble : Double = {
    //println(sgn, ex, man)
    if (isNaN) return Double.NaN
    if (isInfinite) {
      if (sgn!=0) return Double.NegativeInfinity
      else        return Double.PositiveInfinity
    }
    if (isZero) return 0.0
    val manRound =
      if (spec.manW>52) roundToEven(spec.manW-52,man)
      else ( man << (52-spec.manW) )
    val manD = manRound.toDouble
    val z = scalb(manD, ex-spec.exBias-52)
    if (sgn!=0) -z else z
  }

}

object RealGeneric {
  def zero(spec : RealSpec) : RealGeneric = { new RealGeneric(spec, 0, 0, SafeLong(0)) }
  def nan(spec : RealSpec) : RealGeneric = { new RealGeneric(spec, 0, maskI(spec.exW), SafeLong(1)<<(spec.manW-1)) }
  def inf(spec : RealSpec, sgn: Int) : RealGeneric = { new RealGeneric(spec, sgn, maskI(spec.exW), SafeLong(0)) }

  def fromDouble(spec : RealSpec, x : Double ) : RealGeneric = {
    if (x == 0.0) zero(spec)
    else if (x.isNaN) nan(spec)
    else if (x.isPosInfinity) inf(spec, 0)
    else if (x.isNegInfinity) inf(spec, 1)
    else {
      val bb = java.nio.ByteBuffer.allocate(8).putDouble(x).array()
      //bb.foreach(n=>println(f"$n%02x"))
      val lv = bb.foldLeft(0L){ (acc, x) => (acc<<8)+(x&0xFF) }
      //println(f"$lv%08x")
      val sgn = bit(63,lv).toInt
      val man = SafeLong(slice(0,52,lv) + (1L<<52))
      val manRound =
        if (spec.manW>=52) { man << (spec.manW-52) } // left shift
        else { roundToEven( 52 - spec.manW, man ) }
      val moreThan2AfterRound = bit(spec.manW+1, manRound)
      val manRoundNorm = if (moreThan2AfterRound!=0) { manRound >> 1 } else manRound
      val exD = slice(52,11,lv).toInt
      val ex  = exD - 0x3FF + moreThan2AfterRound + spec.exBias
      //println(f"$sgn%d $ex%x")
      if (ex<=0) {
        // subnormal should be considered later
        zero(spec)
      } else if (ex>=maskI(spec.exW)) {
        inf(spec,sgn)
      } else {
        new RealGeneric(spec, sgn, ex, manRoundNorm)
      }
    }
  }

  def fromFloat(spec : RealSpec, x : Float ) : RealGeneric = fromDouble(spec, x.toDouble)
  def fromInt(spec : RealSpec, x : Int ) : RealGeneric = fromDouble(spec, x.toDouble)

}
