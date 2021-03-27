//% @file exponential.scala
//
// Simulators for Exponential functions
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import chisel3._
import chisel3.util._

import rial.table._
import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric


// log2SimGeneric : Log base of 2 function,
//            sign is ignored (assume always positive)
//            x = 2^e * 1.m
//            z = log_2 x = e + log(1.m)
//               never overflow for the same input/output precision

class Log2Sim( val spec : RealSpec, val nOrder : Int, val aW : Int, extraBits : Int ) {
  val adrW  = if (spec.manW<aW) {
    sys.error(f"WARNING: Address width $aW > mantissa width ${spec.manW}")
    spec.manW
  } else {
    aW
  }
  val calcW = spec.manW + extraBits
  val maxZEx = log2Down(max(spec.exMax+1, -spec.exMin))+1
  val t : FuncTableInt = {
    val tableD = new FuncTableDouble( x => log(1.0+x)/log(2.0), nOrder )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, calcW )
  }

  def eval( x: RealGeneric ) : RealGeneric = {
    if (x.spec != spec) {
      sys.error(f"ERROR: input spec != output spec") 
      return RealGeneric.nan(spec)
    }
    val expW = spec.exW
    val manW = spec.manW
    val exBias = spec.exBias

    val man    = x.man
    val exRaw  = x.ex
    val sgn    = x.sgn // Ignored
    val ex     = exRaw-exBias

    if (x.isNaN)      return RealGeneric.nan(spec)
    if (x.isInfinite) return RealGeneric.inf(spec,0)
    if (x.isZero)     return RealGeneric.inf(spec,1) // -Inf

    val dxbp = manW-adrW-1
    val d    = slice(0, dxbp+1, man) - (SafeLong(1)<<dxbp)
    //println(f"d=${d.toLong}%x")
    val adr     = slice(dxbp+1, adrW, man).toInt
    // From here Long is used instead of SafeLong; should be fixed
    val logFrac = t.interval(adr).eval(d.toLong, dxbp)
    val log2 = (ex.toLong<<calcW) + logFrac
    //println(f"ex=${ex}%x logFrac=${logFrac}%x log2=${log2}%x")
    val log2sgn = if (log2 < 0) 1 else 0
    val log2abs = abs(log2)
    val zEx0   = log2Down(log2abs) // Leading 1 pos, -1 for 0
    val (zEx, zMan) =
      if (zEx0 < 0) {
        (0, 0L)
      } else {
        val z0w = maxZEx+calcW // Including leading 1
        val zShift = log2abs << (z0w-(zEx0+1))
        // b(z0w-1) b(z0w-2) ... b(z0w-(manW+1))   b(z0w-manW-2) b(z0w-manW-3) ... b0
        // 1        < --- mantissa (manW) ----->   round         <---- stickey  ---->
        val zman0   = slice(z0w-(manW+1), manW, zShift) // no leading 1
        //println(f"log2abs=${log2abs}%x zShift=${zShift}%x zman0=${zman0}%x")
        val lsb     = zman0 & 1
        val round   = bit(z0w-manW-2, zShift)
        val stickey = if ( slice(0, z0w-manW-2, zShift) == 0 ) 0 else 1
        val inc = round & ( stickey | lsb )
        val zmanR   = zman0 + inc
        val zEx1 = exBias + zEx0 - calcW
        if (bit(manW,zmanR)!=0) {
          ( zEx1+1, zmanR & maskL(manW) )
        } else {
          ( zEx1, zmanR )
        }
      }
    //println(f"zEx=${zEx}%x zMan=${zMan}%x")
    new RealGeneric(x.spec, log2sgn, zEx, SafeLong(zMan))
  }
}

class Log2F32Sim extends Log2Sim( RealSpec.Float32Spec, 2, 8, 2 ) { }
class Log2BF16Sim extends Log2Sim( RealSpec.BFloat16Spec, 0, 7, 3 ) { }

  /*
  def expSimGeneric( t : FuncTableInt, extraBits: Int, x: RealGeneric ) : RealGeneric = {
    val expW = x.spec.exW
    val manW = x.spec.manW
    val exBias = x.spec.exBias

    val man    = x.manW1
    val xe     = x.ex
    val sgn    = x.sgn

    if (x.isNaN || x.isInfinite || x.isZero) {
      pow2simGeneric(t, extraBits, x )
    } else {
      val log2e = 1.0d/java.lang.Math.log(2.0d)
      val log2eMan = BigInt(round(java.lang.Math.scalb(log2e, manW)))
      val prod = man * log2eMan
      val moreThan2 = bit(manW*2+1, prod)
      val ym = slice(manW+moreThan2, manW, prod)
      val r  = bit(manW-1+moreThan2, prod)
      val stickey = if (slice(0, manW-1+moreThan2, prod) != 0) 1 else 0 
      val inc = r & (stickey | bit(0,ym))
      val ymRound = ym + inc
      val ye = if ((bit(manW,ymRound)|moreThan2) !=0 ) xe+1 else xe
      val yman = ymRound & maskSL(manW)
      val y = new RealGeneric(x.spec, sgn, ye, yman)
      //println(f"${prod}%h ${ym}%h ${x.value}%h ${y.value}%h")
      val z = pow2simGeneric(t, extraBits, y)
      //println(f"${z.value}%h ${x.toDouble}%f ${y.toDouble}%f ${z.toDouble}%f")
      z
    }
  }
   */
