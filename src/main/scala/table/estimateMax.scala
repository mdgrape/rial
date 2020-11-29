//
// Estimate Min / Max 
//   Copyright (C) 2019 Makoto Taiji, RIKEN BDR
//
package rial.table

import scala.math._
import java.lang.Math.scalb

// Estimate Minimum or Maximum
//   Since here we assume function is so smooth in the range,
//   2nd derivative is always positive or negative in the range.
//   In this case, the max and min value can be obtained by:
//     * deriv(xMin) and deriv(xMax) have the same sign
//       -> at xMin or xMax
//     * deriv(xMin) and deriv(xMax) have the different sign
//       -> at xMin or xMax or where deriv(x)==0
object estimateMinMax {
  def solveL( f : Long => Long, xMin : Long, xMax : Long) : Long = {
    // via bisection,
    // assume f(xMin) * f(xMax) < 0
    // assume only one solution
    val l = f(xMin)
    val r = f(xMax)
    val xMid = (xMin+xMax)/2L
    //println(xMin, xMax, xMid)
    val c = f(xMid)
    if      (l==0L) xMin
    else if (r==0L) xMax
    else if (c==0L) xMid
    else if (xMax<=xMin+1) xMin
    else if ( ((c<0L) && (l>0L)) || ((c>0L) && (l<0L)) ) solveL (f, xMin, xMid)
    else               solveL (f, xMid, xMax)
  }

  def getMinMaxLong( f : Long => Long, deriv : Long => Long,
    xMin : Long, xMax : Long) = {
    val ld = deriv(xMin)
    val rd = deriv(xMax)
    val lf = f(xMin)
    val rf = f(xMax)
    val max0 = max(lf,rf)
    val min0 = min(lf,rf)
    if ( ((ld<0L) && (rd>0L)) || ((ld>0L) && (rd<0L)) ) {
      val xc = solveL(deriv, xMin, xMax)
      val cf = f(xc)
      val min1 = min(min0,cf)
      val max1 = max(max0,cf)
      if (xc<xMax) {
        val cf2 = f(xc+1)
        ( min(min1, cf2), max(max1, cf2) )
      } else {
        ( min1, max1 )
      }
    } else {
      (min0, max0)
    }
  }

  def solveD( f : Double => Double, xMin : Double, xMax : Double, eps : Double = 1e-16, xeps : Double = 1e-16) : Double = {
    // via bisection,
    // assume f(xMin) * f(xMax) < 0
    // assume only one solution
    val l = f(xMin)
    val r = f(xMax)
    val xMid = (xMin+xMax)/2.0
    //println(xMin, xMax, xMid)
    val c = f(xMid)
    if      (abs(c)<eps) xMid
    else if (xMax<=xMin+xeps) xMid
    else if ( c*l < 0.0 ) solveD (f, xMin, xMid, eps, xeps)
    else                  solveD (f, xMid, xMax, eps, xeps)
  }

  def getMinMaxDouble( f : Double => Double, deriv : Double => Double,
    xMin : Double, xMax : Double, eps : Double = 1e-16, xeps : Double = 1e-16) = {
    val ld = deriv(xMin)
    val rd = deriv(xMax)
    val lf = f(xMin)
    val rf = f(xMax)
    val max0 = max(lf,rf)
    val min0 = min(lf,rf)
    if ( ld * rd < 0.0 ) {
      val xc = solveD(deriv, xMin, xMax, eps, xeps)
      val cf = f(xc)
      val min1 = min(min0,cf)
      val max1 = max(max0,cf)
      ( min1, max1 )
    } else {
      (min0, max0)
    }
  }
  
}
