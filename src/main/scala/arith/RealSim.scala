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

    if (this.isZero || that.isZero) {
      return zero(resSpec)
    } else if (this.isNaN || that.isNaN) {
      if (resSpec.disableNaN) {
        return inf(resSpec, resSgn)
      } else {
        return nan(resSpec)
      }
    } else if (this.isInfinite || that.isInfinite) {
      return inf(resSpec, resSgn)
    }

    // Normal values

    val prod= this.manW1 * that.manW1
    val bp = this.spec.manW + that.spec.manW
    val moreThan2 = bit(bp+1,prod)
    val roundbits = bp-resSpec.manW+moreThan2
    val prodRound = if(roundbits >= 0) {
      roundBySpec(roundSpec, roundbits, prod)
    } else {
      prod << (-roundbits)
    }
    val moreThan2AfterRound = bit(resSpec.manW+1, prodRound)
    val resMan = prodRound >> moreThan2AfterRound

    val ex0 = (this.ex-this.spec.exBias) + (that.ex-that.spec.exBias)
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
        if (this.isNaN) {
          return inf(resSpec, sgn)
        } else {
          return inf(resSpec, that.sgn)
        }
      } else {
        return nan(resSpec)
      }
    }
    if (this.isInfinite && that.isInfinite) {
      if (this.sgn != that.sgn) { // inf - inf = nan, inf + inf = inf
        if(resSpec.disableNaN) {
          return inf(resSpec, this.sgn)
        } else {
          return nan(resSpec)
        }
      } else {
        return inf(resSpec, this.sgn)
      }
    }
    if (this.isInfinite) {return inf(resSpec, this.sgn)}
    if (that.isInfinite) {return inf(resSpec, that.sgn)}

    if (this.isZero && that.isZero) {
      val resSgn = if(roundSpec == RoundSpec.truncate) {
        (this.sgn, that.sgn) match {
          case (0, 0) => 0
          case _      => 1
        }
      } else {
        (this.sgn, that.sgn) match {
          case (1, 1) => 1
          case _      => 0
        }
      }
      return RealGeneric.zero(resSpec, resSgn)
    }
    if(this.isZero) {
      // convert that.spec to resSpec
      val resSgn = that.sgn
      val resEx0 = that.ex - that.spec.exBias + resSpec.exBias
      val (resMan, resExInc) = if(resSpec.manW >= that.spec.manW) {
        (that.man << (resSpec.manW - that.spec.manW), 0)
      } else {
        val thatRounded = roundBySpec(roundSpec, that.spec.manW - resSpec.manW, that.man)
        val moreThan2AfterRound = bit(resSpec.manW+1, thatRounded)
        (slice(0, resSpec.manW, thatRounded), moreThan2AfterRound)
      }
      val resEx = resEx0 + resExInc
      if(resEx >= maskI(resSpec.exW)) {
        return RealGeneric.inf(resSpec, resSgn)
      }
      return new RealGeneric(resSpec, resSgn, resEx, resMan)
    }
    if(that.isZero) {
      val resSgn = this.sgn
      val resEx0 = this.ex - this.spec.exBias + resSpec.exBias
      val (resMan, resExInc) = if(resSpec.manW >= this.spec.manW) {
        (this.man << (resSpec.manW - this.spec.manW), 0)
      } else {
        val thisRounded = roundBySpec(roundSpec, this.spec.manW - resSpec.manW, this.man)
        val moreThan2AfterRound = bit(resSpec.manW+1, thisRounded)
        (slice(0, resSpec.manW, thisRounded), moreThan2AfterRound)
      }
      val resEx = resEx0 + resExInc
      if(resEx >= maskI(resSpec.exW)) {
        return RealGeneric.inf(resSpec, resSgn)
      }
      return new RealGeneric(resSpec, resSgn, resEx, resMan)
    }

    // Normal values
    val thisEx = this.ex - this.spec.exBias
    val thatEx = that.ex - that.spec.exBias
    val (x,y,exMax) = if (thisEx<thatEx) {
      (that, this, thatEx)
    } else {
      (this, that, thisEx)
    }

    // resolve difference in mantissa width
    val (xmanW1, ymanW1) = if(x.spec.manW == y.spec.manW) {
      (x.manW1, y.manW1)
    } else if (x.spec.manW > y.spec.manW) {
      (x.manW1, y.manW1 << (x.spec.manW - y.spec.manW))
    } else {
      (x.manW1 << (y.spec.manW - x.spec.manW), y.manW1)
    }

    // align x with y
    val diffEx = abs(thisEx-thatEx)
    val xshift = xmanW1 << diffEx
    val sum =  if (this.sgn!=that.sgn) {
      xshift - ymanW1
    } else {
      xshift + ymanW1
    }
    //println(s"Sim: ${sum.toLong.toBinaryString}")

    // In case of diffEx==0, this can be negative
    val (sumPos, resSgn) = if (sum<0) {(-sum, x.sgn^1)} else {(sum,x.sgn)}

    if (sumPos==0) {return zero(resSpec)}

    // Normalization
    val leading1 = sumPos.bitLength-1
    //val ml=sumPos.toLong; println(f"$ml%x $leading1%d $diffEx%d")

    val roundbits = leading1 - resSpec.manW
    val (resMan, exInc) =
      if (roundbits>0) {
        val sumRound = roundBySpec(roundSpec, roundbits, sumPos)
        val moreThan2AfterRound = bit(resSpec.manW+1, sumRound)
        (sumRound >> moreThan2AfterRound, moreThan2AfterRound)
      } else {
        ( sumPos << (-roundbits), 0)
      }

    // no exponent change if leading1 = x.manW+diffEx
    val resEx = exMax + leading1 - max(x.spec.manW, y.spec.manW) - diffEx + exInc + resSpec.exBias

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
    val xyMoreThan2 = bit((xyWidth+2)-1, xyProd)
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
        return new RealGeneric(resSpec, resSgn, 0, 0) // underflow to zero but keep xySgn.
      } else if (resEx >= maskI(resSpec.exW)) { // result is inf
        return inf(resSpec, resSgn)
      } else {
        return new RealGeneric(resSpec, resSgn, resEx, resMan)
      }
    }

    val zExNobias = z.ex - z.spec.exBias
    val diffEx    = abs(xyExNobias - zExNobias)
    val diffManW  = xyManWidth - z.spec.manW // > 0, normally

//     println(s"--------------------------")
//     println(s"xyExNobias = ${xyExNobias}") // 78
//     println(s" zExNobias = ${zExNobias}")  //-60
//     println(s"diffManW   = ${diffManW}")   // 24
//     println(s"xyprod     = ${xyProd.toLong.toBinaryString}")
//     println(s"z          = ${z.manW1.toLong.toBinaryString}")

    val (sum, ex, sumManWidth) = if (xyExNobias > zExNobias && diffEx > diffManW) {
      // xyEx > zEx, and the difference exceeds the difference in mantissa width
      // normally, width of xyprod is larger than the width of z mantissa. so,
      // to align xyprod and z, we shift z.
      //   xx.xxxxxxxx*2^3     xx.xxxxxxxx
      //         z.zzz*2^2 ->     z.zzz <<
      //
      // But if the difference is larger than the difference in mantissa width,
      // we need to shift xyprod
      //
      //   xx.xxxxxxxx*2^10     xx.xxxxxxxx0000 <<
      //         z.zzz*2^2  ->            z.zzz

      val xyshift = xyProd << (diffEx - diffManW)
//       println(s"xyshift = ${xyshift.toLong.toBinaryString}")
      if (xySgn == z.sgn) {
        (xyshift + z.manW1, xyExNobias, xyManWidth + diffEx - diffManW)
      } else {
        (xyshift - z.manW1, xyExNobias, xyManWidth + diffEx - diffManW)
      }
    } else {
      val zshift = z.manW1 << (diffManW + zExNobias - xyExNobias)
//       println(s"zshift     = ${zshift.toLong.toBinaryString}")
      if (xySgn == z.sgn) {
        (xyProd + zshift, xyExNobias, xyManWidth)
      } else {
        (xyProd - zshift, xyExNobias, xyManWidth)
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

//     println(s"resMan = ${resMan.toLong.toBinaryString}(${resMan}))")
//     println(s"ex         = ${ex}")
//     println(s"leading1   = ${leading1}")
//     println(s"xyManWidth = ${xyManWidth}")
//     println(s"exInc      = ${exInc}")
//     println(s"exBias     = ${resSpec.exBias}")

    val roundExInc = leading1 - sumManWidth
    val resEx = ex + exInc + resSpec.exBias + roundExInc
//     println(s"roundExInc     = ${roundExInc}")
//     println(s"resEx  = ${resEx}")

    if (resEx <= 0) {
      zero(resSpec) // subnormal not supported
    } else if (resEx >= maskI(resSpec.exW)) { // result is inf
      inf(resSpec, resSgn)
    } else {
      new RealGeneric(resSpec, resSgn, resEx, resMan)
    }
  }

  def add3op( resSpec : RealSpec, roundSpec : RoundSpec, y : RealGeneric, z : RealGeneric ) : RealGeneric = {

    val (xnan,  ynan,  znan)  = (this.isNaN, y.isNaN, z.isNaN)
    val (xinf,  yinf,  zinf)  = (this.isInfinite, y.isInfinite, z.isInfinite)
    val (xzero, yzero, zzero) = (this.isZero, y.isZero, z.isZero)

    val (xsgn,   ysgn,   zsgn)   = (this.sgn, y.sgn, z.sgn)
    val (xex,    yex,    zex)    = (this.ex,  y.ex,  z.ex)
    val (xmanW1, ymanW1, zmanW1) = (this.manW1, y.manW1, z.manW1)

    val (xSpec, ySpec, zSpec) = (this.spec, y.spec, z.spec)

    val wnan = xnan || ynan || znan ||
              (xinf && yinf && (xsgn != ysgn)) ||
              (yinf && zinf && (ysgn != zsgn)) ||
              (zinf && xinf && (zsgn != xsgn))

    val winf0 = (xinf || yinf || zinf) && !wnan
    val winf0sgn = ! ((xinf && (xsgn != 1)) || (yinf && (ysgn != 1)) || (zinf && (zsgn != 1)))

    val wzero0 = xzero && yzero && zzero
    val wzero0sgn = if(roundSpec == RoundSpec.truncate) {
      if(xsgn == 1 || ysgn == 1 || zsgn == 1) {1} else {0}
    } else {
      if(xsgn == 1 && ysgn == 1 && zsgn == 1) {1} else {0}
    }

    if (wnan) {
      if (resSpec.disableNaN) {
        return inf(resSpec, 0)
      } else {
        return nan(resSpec)
      }
    }
    if (winf0) {
      return inf(resSpec, (if(winf0sgn) {1} else {0}))
    }
    if (wzero0) {
      return zero(resSpec, wzero0sgn)
    }

    val fracW = Seq(xSpec.manW, ySpec.manW, zSpec.manW).max
    val maxEx = Seq(xex - xSpec.exBias, yex - ySpec.exBias, zex - zSpec.exBias).max

    // 1 for leading1, 2 for guard and rounding bits (stickies follows later)
    val xmanPad = (1+fracW+fracW+2 - xSpec.manW)
    val ymanPad = (1+fracW+fracW+2 - ySpec.manW)
    val zmanPad = (1+fracW+fracW+2 - zSpec.manW)

    val xmanW1padded = if(xzero) {SafeLong(0)} else {xmanW1 << xmanPad}
    val ymanW1padded = if(yzero) {SafeLong(0)} else {ymanW1 << ymanPad}
    val zmanW1padded = if(zzero) {SafeLong(0)} else {zmanW1 << zmanPad}

    val xShift = maxEx - (xex - xSpec.exBias)
    val yShift = maxEx - (yex - ySpec.exBias)
    val zShift = maxEx - (zex - zSpec.exBias)

    val xSticky = if((xmanW1padded & maskSL(xShift - xmanPad)) != 0) {1} else {0}
    val ySticky = if((ymanW1padded & maskSL(yShift - ymanPad)) != 0) {1} else {0}
    val zSticky = if((zmanW1padded & maskSL(zShift - zmanPad)) != 0) {1} else {0}

    val xmanW1Aligned =  ((xmanW1padded >> xShift) << 1) + xSticky
    val ymanW1Aligned = (((ymanW1padded >> yShift) << 1) + ySticky) * (if(xsgn != ysgn) {-1} else {1})
    val zmanW1Aligned = (((zmanW1padded >> zShift) << 1) + zSticky) * (if(xsgn != zsgn) {-1} else {1})

    // note that the circuit uses 3:2 CSA and early normalization.
    val sumAligned  = xmanW1Aligned + ymanW1Aligned + zmanW1Aligned
    val sumNegative = sumAligned <  0
    val sumZero     = sumAligned == 0
    val sumPos      = if(sumNegative) {-sumAligned} else {sumAligned}

    val resSgn = if(sumNegative) {xsgn ^ 1} else {xsgn}

    val sumLeadingOne = sumPos.bitLength - 1
    val sumRoundBits  = sumLeadingOne - resSpec.manW
    val (resMan, exInc) = if(sumRoundBits > 0){
      val rounded = roundBySpec(roundSpec, sumRoundBits, sumPos)
      val moreThan2AfterRound = bit(resSpec.manW+1, rounded)
      (rounded >> moreThan2AfterRound, moreThan2AfterRound)
    } else {
      (sumPos << (-sumRoundBits), 0)
    }

    // +3 means +2+sticky
    val resEx = maxEx + (sumLeadingOne - (1+fracW+fracW+3)) + exInc + resSpec.exBias
    if(resEx <= 0 || sumZero) {
      zero(resSpec)
    } else if(resEx >= maskI(resSpec.exW)) {
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

  def toFloat : Float = {
    if (spec == RealSpec.Float32Spec) {
      val longBits = value & maskSL(spec.W)
      java.lang.Float.intBitsToFloat(longBits.toInt)
    } else {
      //println(sgn, ex, man)
      if (isNaN) {return Float.NaN}
      if (isInfinite) {
        if (sgn!=0) return Float.NegativeInfinity
        else        return Float.PositiveInfinity
      }
      if (isZero) {return 0.0f}

      val manRound = if (spec.manW>23) {
        roundToEven(spec.manW-23, manW1)
      } else if (spec.manW < 23) {
        manW1 << (23-spec.manW)
      } else {
        manW1
      }

      val manD = manRound.toFloat
      val z = scalb(manD, ex-spec.exBias-23)
      if (sgn!=0) {-z} else {z}
    }
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
    this.setSign( 1-this.sgn )
  }

}

object RealGeneric {
  def zero(spec : RealSpec, sgn: Int = 0) : RealGeneric = { new RealGeneric(spec, sgn, 0, SafeLong(0)) }
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
    if (x == 0.0) {
      val (sgn,exD,manD) = unpackDouble(x)
      if(sgn == 1) {
        SafeLong(1 << (spec.W-1))
      } else {
        SafeLong(0)
      }
    }
    else if (x.isNaN)         { maskSL(spec.exW+1) << (spec.manW-1) }
    else if (x.isPosInfinity) { maskSL(spec.exW  ) <<  spec.manW }
    else if (x.isNegInfinity) { maskSL(spec.exW+1) <<  spec.manW }
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
