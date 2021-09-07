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

  def W : Int = { exW+manW+(if (disableSign) 0 else 1) }

  def exMax : Int = { maskI(exW)-1-exBias }
  def exMin : Int = { 1-exBias }

  //def toString : String = s"(${exW},${manW}) exBias=${exBias}" +
  //s"disableSign=$disableSign disableNaN=$disableNaN disableSubNormal=$disableSubnormal"
  def toStringShort : String = s"(${exW},${manW})"

  override def equals(that: Any): Boolean =
    that match {
      case that: RealSpec => (that.exW==this.exW) &&
        (that.exBias==this.exBias) && (that.manW==this.manW)
      case _ => false
    }
}

object RealSpec {
  val Float64Spec = new RealSpec(11, 0x3FF, 52, false, false, true)
  val Float32Spec = new RealSpec(8, 0x7F, 23, false, false, true)
  val Float16Spec = new RealSpec(5, 0xF, 10, false, false, true)
  val BFloat16Spec = new RealSpec(8, 0x7F, 7, false, false, true)
}

//
// Again we do not use type generics for shift operators here...
//
// class RealGeneric
//   this uses SafeLong for representation
//
class RealGeneric ( val spec : RealSpec, val value: SafeLong  ) {

  def this( spec : RealSpec, s : Int, e : Int, m : SafeLong ) = {
    this(spec, (((SafeLong(s&1)<<spec.exW)+SafeLong(e&maskI(spec.exW)))<<spec.manW)+(m&maskSL(spec.manW)))
  }

  def this( spec : RealSpec, d : Double ) = {
    this(spec, doubleToRealGeneric(spec,d) )
  }

  def man   = (value & maskSL(spec.manW))
  def manW1 = ( (value & maskSL(spec.manW)) + (SafeLong(1)<<spec.manW) ) // with leading 1
  def ex  = slice(spec.manW,spec.exW,value).toInt
  def sgn = bit(spec.manW+spec.exW,value).toInt
  def exNorm = ex - spec.exBias

  def fracIsZero : Boolean = { man == 0 }

  def isNaN : Boolean = { (ex == maskI(spec.exW)) && !fracIsZero && !spec.disableNaN }
  def isInfinite : Boolean = { (ex == maskI(spec.exW)) && ( fracIsZero || spec.disableNaN ) }
  def isZero : Boolean = { (ex == 0) && (spec.disableSubnormal || fracIsZero) }
  def isSubnormal : Boolean = { (ex == 0) && !fracIsZero }

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
    val prod= manW1 * that.manW1
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
    if (this.isInfinite && that.isInfinite) {
      if (this.sgn == that.sgn) {
        if(resSpec.disableNaN) {
          return inf(resSpec, this.sgn)
        } else return nan(resSpec)
      } else return inf(resSpec, this.sgn)
    }
    if (this.isInfinite) return inf(resSpec, sgn)
    if (that.isInfinite) return inf(resSpec, that.sgn)

    // Normal values
    val thisEx=ex-this.spec.exBias
    val thatEx=that.ex-that.spec.exBias
    val (x,y,exMax) = if (thisEx<thatEx) (that, this, thatEx) else (this, that, thisEx)
    val diffEx = abs(thisEx-thatEx)
    val xshift = x.manW1 << diffEx
    val sum =  if (this.sgn!=that.sgn) xshift-y.manW1 else xshift+y.manW1
    //println(s"Sim: ${sum.toLong.toBinaryString}")
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

  // FMA
  //
  // - any of the arguments is NaN => NaN
  // - if z is Inf and others are finite => Inf
  //   - Note: while FMA operation, `x * y` has infinite precision, so overflow
  //           never happens. -2^1000 * 2^1000 + Inf == Inf, not NaN.
  // - otherwise => same as `x * y + z`
  //   - 0 * Inf = NaN
  //   - Inf - Inf = NaN
  //
  // TODO: currently subnormal is not supported
  //
  def fmadd( resSpec : RealSpec, roundSpec : RoundSpec, y : RealGeneric, z : RealGeneric ) : RealGeneric = {
    if (this.isNaN || y.isNaN || z.isNaN) {
      if (resSpec.disableNaN) {
        if (this.isNaN) {
          return inf(resSpec, sgn)
        } else if (y.isNaN) {
          return inf(resSpec, y.sgn)
        } else {
          return inf(resSpec, z.sgn)
        }
      } else {
        return nan(resSpec)
      }
    }
    val xySgn = this.sgn ^ y.sgn
//     println("xySgn = ", xySgn)

    if (z.isInfinite) {
      if (!this.isInfinite && !y.isInfinite) { // finite + Inf = Inf
        return inf(resSpec, z.sgn)
      } else if (xySgn == z.sgn) { // Inf + Inf = Inf
        return inf(resSpec, z.sgn)
      } else {                     // Inf - Inf = NaN
        if (resSpec.disableNaN) {
          return inf(resSpec, z.sgn) // sign?
        } else {
          return nan(resSpec)
        }
      }
    }

    if (this.isInfinite || y.isInfinite) { // here z is finite
      return inf(resSpec, xySgn)
    }

    if (this.isZero || y.isZero) {
      if (z.isZero) {
        if (roundSpec == RoundSpec.truncate) {
          if (xySgn == 0 && z.sgn == 0) { // +0 + +0 = +0
            zero(resSpec)
          } else {
            new RealGeneric( resSpec, 1, 0, 0 ) // -0
          }
        } else {
          if (xySgn == 1 && z.sgn == 1) {
            return z      // -0 + -0 = -0
          } else {
            zero(resSpec) // otherwise +0
          }
        }
      } else return z // if z is not zero, then return z.
    }

    val xyEx0       = (this.ex - this.spec.exBias) + (y.ex - y.spec.exBias)
    val xyProd      = this.manW1 * y.manW1
    val xyWidth     = this.spec.manW + y.spec.manW
    val xyMoreThan2 = bit(xyWidth+1, xyProd)
    val xyManWidth  = xyWidth + xyMoreThan2
    val xyExNobias  = xyEx0   + xyMoreThan2

//     println(s"xManW1     = ${this.manW1.toLong.toBinaryString}")
//     println(s"yManW1     = ${   y.manW1.toLong.toBinaryString}")
//     println(s"xyProd     = ${xyProd.toLong.toBinaryString}")
//     println(s"xyMoreThan2= ${xyMoreThan2}")

    if(z.isZero) {

      val xyRoundbits = xyManWidth - resSpec.manW
      val xyProdRound = roundBySpec(roundSpec, xyRoundbits, xyProd)
      val moreThan2AfterRound = bit(resSpec.manW+1, xyProdRound)
      val resMan = xyProdRound >> moreThan2AfterRound
      val resEx  = xyExNobias + moreThan2AfterRound + resSpec.exBias
      val resSgn = xySgn

      if (resEx <= 0) {
        return zero(resSpec)
      } else if (resEx >= maskI(resSpec.exW)) { // result is inf
        return inf(resSpec, resSgn)
      } else {
        return new RealGeneric(resSpec, resSgn, resEx, resMan)
      }
    }

    val zExNobias = z.ex - z.spec.exBias
    val diffEx    = abs(xyExNobias - zExNobias)
    val diffManW  = xyManWidth - z.spec.manW // > 0, normally

//     println(s"xyExNobias = ${xyExNobias}")
//     println(s" zExNobias = ${zExNobias}")

    val (sum, ex) = if (xyExNobias > zExNobias && diffEx > diffManW) {
      val xyshift = xyProd << (diffEx - diffManW)
//       println(s"xyshift = ${xyshift.toLong.toBinaryString}")
      if (xySgn == z.sgn) {
        (xyshift + z.manW1, zExNobias)
      } else {
        (xyshift - z.manW1, zExNobias)
      }
    } else {
      val zshift = z.manW1 << (diffManW + zExNobias - xyExNobias)
//       println(s"zshift     = ${zshift.toLong.toBinaryString}")
      if (xySgn == z.sgn) {
        (xyProd + zshift, xyExNobias)
      } else {
        (xyProd - zshift, xyExNobias)
      }
    }
//     println(s"sum = ${sum.toLong.toBinaryString}")
//     println(s"ex  = ${ex}")

    val (sumPos, resSgn) = if (sum < 0) { (-sum, xySgn^1) } else { (sum, xySgn) }

    if (sumPos == 0) {
      return zero(resSpec)
    }

    // normalize and round

    val leading1  = sumPos.bitLength-1
    val roundbits = leading1 - resSpec.manW
    val (resMan, exInc) = if (roundbits > 0) {
      val sumRound = roundBySpec(roundSpec, roundbits, sumPos)
      val moreThan2AfterRound = bit(resSpec.manW+1, sumRound)
      (sumRound >> moreThan2AfterRound, moreThan2AfterRound)
    } else {
      (sumPos << (-roundbits), 0)
    }

//     println(s"resMan = ${resMan.toLong.toBinaryString}")

//     println(s"ex         = ${ex}")
//     println(s"leading1   = ${leading1}")
//     println(s"xyManWidth = ${xyManWidth}")
//     println(s"exInc      = ${exInc}")
//     println(s"exBias     = ${resSpec.exBias}")
    val resEx = ex + (leading1 - xyManWidth) + exInc + resSpec.exBias
//     println(s"resEx  = ${resEx}")

    if (resEx <= 0) {
      zero(resSpec) // subnormal not supported
    } else if (resEx >= maskI(resSpec.exW)) { // result is inf
      inf(resSpec, resSgn)
    } else {
      new RealGeneric(resSpec, resSgn, resEx, resMan)
    }
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
      if (spec.manW>52) roundToEven(spec.manW-52,manW1)
      else ( manW1 << (52-spec.manW) )
    val manD = manRound.toDouble
    val z = scalb(manD, ex-spec.exBias-52)
    if (sgn!=0) -z else z
  }

  def convert( newSpec : RealSpec, roundSpec : RoundSpec ) : RealGeneric = {
    val manRound =
      if (newSpec.manW>=spec.manW) { man << (newSpec.manW - spec.manW) }
      else { roundBySpec( roundSpec, spec.manW - newSpec.manW, man ) }
    val moreThan2AfterRound = bit(newSpec.manW, manRound)
    val manRoundNorm = if (moreThan2AfterRound!=0) SafeLong(0) else manRound

    val exNew = if (isInfinite || isNaN) {
      maskI(newSpec.exW)
    } else {
      min(maskI(newSpec.exW), ex - spec.exBias + moreThan2AfterRound + newSpec.exBias)
    }

    if (isInfinite || isNaN) {
      // avoid exponent check
      // TODO: consider behavior of converter in case of NaN/Inf, especially mantissa of NaN
      new RealGeneric(newSpec, packToSafeLongWOCheck(newSpec, sgn, exNew, manRoundNorm))
    } else if (exNew == 0) { // to keep sign of zero
      new RealGeneric(newSpec, SafeLong((sgn&1L) << (newSpec.W-1)))
    } else {
      new RealGeneric(newSpec, packToSafeLong(newSpec, sgn, exNew, manRoundNorm))
    }
  }

  def scalbn( n : Int ) : RealGeneric = {
    if (this.isZero || this.isNaN || this.isInfinite) return this
    val exNew = this.ex + n
    if (exNew>=maskI(this.spec.exW)) return inf(this.spec, this.sgn)
    if (exNew<0) return zero(this.spec)
    new RealGeneric( this.spec, this.sgn, exNew, this.man )
  }

  def setSign( n : Int ) : RealGeneric = {
    val sgn = if (n!=0) 1 else 0
    new RealGeneric( this.spec, sgn, this.ex, this.man )
  }

  def negate( ) : RealGeneric = {
    this.setSign( 1-this.ex )
  }

}

object RealGeneric {
  def zero(spec : RealSpec) : RealGeneric = { new RealGeneric(spec, SafeLong(0)) }
  def nan(spec : RealSpec) : RealGeneric = { new RealGeneric(spec, 0, maskI(spec.exW), SafeLong(1)<<(spec.manW-1)) }
  def inf(spec : RealSpec, sgn: Int) : RealGeneric = { new RealGeneric(spec, sgn, maskI(spec.exW), SafeLong(0)) }

  // With exponent check
  private def packToSafeLong( spec: RealSpec, sgn : Int, ex: Int, man : SafeLong) : SafeLong = {
    if (ex<=0) {
      // subnormal should be considered later
      // keep sign of zero
      SafeLong((sgn&1L) << (spec.W-1))
    } else if (ex>=maskI(spec.exW)) {// Inf
      ( maskI(spec.exW+(sgn&1))<<spec.manW )
    } else {
      (man & maskSL(spec.manW))+ ( SafeLong( ex + ((sgn&1)<<spec.exW) ) << spec.manW )
    }
  }
  private def packToSafeLongWOCheck( spec: RealSpec, sgn : Int, ex : Int, man : SafeLong) : SafeLong = {
    (man & maskSL(spec.manW)) + ( SafeLong( ex + ((sgn&1)<<spec.exW) ) << spec.manW )
  }

  private def doubleToRealGeneric(spec : RealSpec, x : Double ) : SafeLong = {
    if (x == 0.0) SafeLong(0)
    else if (x.isNaN) (maskSL(spec.exW+1)<<(spec.manW-1))
    else if (x.isPosInfinity) ( maskSL(spec.exW)<<spec.manW )
    else if (x.isNegInfinity) ( maskSL(spec.exW+1)<<spec.manW )
    else {
      val (sgn,exD,manD) = unpackDouble(x)
      val man = SafeLong(manD)
      val manRound =
        if (spec.manW>=52) { man << (spec.manW-52) } // left shift
        else { roundToEven( 52 - spec.manW, man ) }
      val moreThan2AfterRound = bit(spec.manW, manRound)
      val manRoundNorm = if (moreThan2AfterRound!=0) SafeLong(0) else manRound
      val ex  = exD - 0x3FF + moreThan2AfterRound + spec.exBias
      //val manW = spec.manW; println(f"$x%f $sgn%d $ex%x $manD%x $manW%d")
      packToSafeLong(spec, sgn, ex, manRoundNorm)
    }
  }

  def fromDouble(spec : RealSpec, x : Double ) : RealGeneric = {
    val value = doubleToRealGeneric(spec,x)
    //val lv = value.toLong; println(f"$x%f $lv%x")
    new RealGeneric(spec, value)

  }

  def fromFloat(spec : RealSpec, x : Float ) : RealGeneric = fromDouble(spec, x.toDouble)
  def fromInt(spec : RealSpec, x : Int ) : RealGeneric = fromDouble(spec, x.toDouble)

}
