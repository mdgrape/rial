//% @file atan2Sim.scala
//
// Simulators for atan2(y, x) function
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.table._
import rial.util._
import rial.util.ScalaUtil._

import rial.arith.RealSpec
import rial.arith.RealGeneric
import rial.arith.Rounding._
import rial.arith._

object ATan2Sim {

  // atan2(y, x) = atan(y/x)
  //
  // atan(x) = x - x^3/3 + x^5/5
  //         < x - x^3/2 + x^5/4
  //
  // linearThreshold : x * 2^-manW > x^3/3 (in FP32, ex <= -13)
  //
  // if x > 2^manW:
  //   atan(x) = pi/2 - atan(1/x) ~ pi/2 - 2^-manW ~ pi/2
  // if x > 2^(-linearThreshold):
  //   atan(x) = pi/2 - atan(1/x) ~ pi/2 - 1/x
  //
  // we need table for 2^-8 ~ 2^8 range for FP32
  //
  // atan2(y, x) =
  //    atan(|y|/|x|)      .. x>0, y>0, |x|>|y|
  //   -atan(|y|/|x|)      .. x>0, y<0, |x|>|y|
  //   -atan(|y|/|x|) + pi .. x<0, y>0, |x|>|y|
  //    atan(|y|/|x|) - pi .. x<0, y<0, |x|>|y|
  //
  //    pi/2-atan(|x|/|y|) .. x>0, y>0, |y|>|x|
  //   -pi/2+atan(|x|/|y|) .. x>0, y<0, |y|>|x|
  //    pi/2+atan(|x|/|y|) .. x<0, y>0, |y|>|x|
  //   -pi/2-atan(|x|/|y|) .. x<0, y<0, |y|>|x|
  //
  // = ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|
  //   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y|
  //   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y|
  //   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y|
  //           ^^^^^^^^^^^^^^^^^^^
  //           this part is always positive in all the cases
  //
  def atan2SimGeneric( t_rec : FuncTableInt, ts : Seq[FuncTableIntFixedWidth], y : RealGeneric, x : RealGeneric ) : RealGeneric = {
//     println("==================================================")
    assert(x.spec == y.spec)

    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val xsgn   = x.sgn
    val ysgn   = y.sgn

//     println(f"x = ${x.sgn}|${x.ex}(${xex})|${xman.toLong.toBinaryString}")
//     println(f"y = ${y.sgn}|${y.ex}(${yex})|${yman.toLong.toBinaryString}")
//     println(f"x = ${x.toDouble}, y = ${y.toDouble}, atan2(y, x) = ${atan2(y.toDouble, x.toDouble)}")

    if(x.isNaN || y.isNaN) {
      return RealGeneric.nan(x.spec)
    }
    if(x.isInfinite) {
      return if (y.isInfinite) {
        (xsgn, ysgn) match {
          case    (0, 0)  => new RealGeneric(x.spec,  Pi * 0.25)
          case    (0, 1)  => new RealGeneric(x.spec, -Pi * 0.25)
          case    (1, 0)  => new RealGeneric(x.spec,  Pi * 0.75)
          case _/* 1, 1 */=> new RealGeneric(x.spec, -Pi * 0.75)
        }
      } else {
        (xsgn, ysgn) match {
          case    (0, 0)  => new RealGeneric(x.spec,   0)
          case    (0, 1)  => new RealGeneric(x.spec,  -0)
          case    (1, 0)  => new RealGeneric(x.spec,  Pi)
          case _/* 1, 1 */=> new RealGeneric(x.spec, -Pi)
        }
      }
    }

    if(y.isInfinite) {
      return (xsgn, ysgn) match {
          case    (0, 0)  => new RealGeneric(x.spec,  Pi * 0.5)
          case    (0, 1)  => new RealGeneric(x.spec, -Pi * 0.5)
          case    (1, 0)  => new RealGeneric(x.spec,  Pi * 0.5)
          case _/* 1, 1 */=> new RealGeneric(x.spec, -Pi * 0.5)
      }
    }

    if(x.isZero) {
      if(y.isZero) {
        return RealGeneric.nan(x.spec)
      }
      return (xsgn, ysgn) match {
          case    (0, 0)  => new RealGeneric(x.spec,  Pi * 0.5)
          case    (0, 1)  => new RealGeneric(x.spec, -Pi * 0.5)
          case    (1, 0)  => new RealGeneric(x.spec,  Pi * 0.5)
          case _/* 1, 1 */=> new RealGeneric(x.spec, -Pi * 0.5)
      }
    }
    if(y.isZero) {
      return (xsgn, ysgn) match {
        case    (0, 0)  => new RealGeneric(x.spec,   0)
        case    (0, 1)  => new RealGeneric(x.spec,  -0)
        case    (1, 0)  => new RealGeneric(x.spec,  Pi)
        case _/* 1, 1 */=> new RealGeneric(x.spec, -Pi)
      }
    }

    // ==================================================================
    // calc y / x and x / y

    val bpRec        = t_rec.bp
    val extraBitsRec = bpRec - manW
    val calcWRec     = bpRec

    // calc y / x if |x|>|y|.

    val swapped = slice(0, x.spec.W-1, x.value) < slice(0, x.spec.W-1, y.value)
    val xexS  = if(swapped) {y.ex } else {x.ex}
    val yexS  = if(swapped) {x.ex } else {y.ex}
    val xmanS = if(swapped) {y.man} else {x.man}
    val ymanS = if(swapped) {x.man} else {y.man}

    val y_over_x_ex0 = yexS - xexS
    val x_y_same_man = xmanS == ymanS

//     println(f"y_over_x_sgn = ${y_over_x_sgn}")
//     println(f"y_over_x_ex0 = ${y_over_x_ex0}")

    // 1+manW+1 width (x2)
    val one_over_x_manW1 = if (t_rec.nOrder == 0) {
      val adr  = xmanS
      val res0 = t_rec.interval(adr.toInt).eval(0L, 0) // in (0.5, 1]
      if(xmanS == 0) {
        1 << calcWRec // just 1
      } else {
        res0 // 2^-1 * 1.m [0.5, 1)
      }
    } else {
      val dxbp = manW - t_rec.adrW - 1
      val d    = (SafeLong(1) << dxbp) - slice(0, manW - t_rec.adrW, xmanS)-1
      val adr  = (1 << t_rec.adrW) - 1 - slice(manW - t_rec.adrW, t_rec.adrW, xmanS).toInt
      val res0 = t_rec.interval(adr.toInt).eval(d.toLong, dxbp)

      if(xmanS == 0) {
        1 << calcWRec // just 1
      } else {
        res0 // 2^-1 * 1.m [0.5, 1)
      }
    }
//     println(f"1/x ref = ${(1<<manW).toDouble/((1<<manW) + x.man.toInt)}")
//     println(f"1/x sim = ${one_over_x_manW1.toDouble / (1<<calcWRec)}")

    val y_over_x_man0 = if(x_y_same_man) { 1L<<(manW+calcWRec) } else {
      (one_over_x_manW1 * ((1<<manW) + ymanS)).toLong
    }
    val y_over_x_man0_clz = (calcWRec) + (1+manW) - y_over_x_man0.toLong.toBinaryString.length
    assert(y_over_x_man0_clz <= 1)

    val y_over_x_manW1     = roundBySpec(RoundSpec.round, calcWRec - y_over_x_man0_clz, y_over_x_man0)
    val y_over_x_MoreThan2 = bit(manW+1, y_over_x_manW1)
    val y_over_x_man       = if (y_over_x_MoreThan2 == 1) {
      ((y_over_x_manW1 >> 1) + bit(0, y_over_x_manW1)) - (1 << manW)
    } else {
      y_over_x_manW1 - (1 << manW)
    }
    val y_over_x_ex     = y_over_x_ex0 - y_over_x_man0_clz + y_over_x_MoreThan2.toInt
    val y_over_x_exBias = y_over_x_ex + exBias
//     println(f"y_over_x_ex  = ${y_over_x_ex}")
//     println(f"y_over_x = ref(${y.toDouble/x.toDouble}) = ${new RealGeneric(x.spec, y_over_x_sgn, y_over_x_ex.toInt+exBias, y_over_x_man).toDouble}")

    // ==================================================================
    // calc atan(y/x) = pi/2 - atan(x/y)

    val linearThreshold = -math.round(math.ceil(manW / 2.0 + 1.0))
//     println(f"constant threshold = ${manW}")
//     println(f"linear threshold   = ${linearThreshold}")

    val adrW      = ts(0).adrW
    val nOrder    = ts(0).nOrder
    val bp        = ts(0).bp
    val extraBits = bp - manW
    val calcW     = manW + extraBits

    assert(y_over_x_ex < 0 || (y_over_x_ex == 0 && y_over_x_man == 0))

    val (atanyx_ex, atanyx_manW1) = if (y_over_x_ex < linearThreshold) {
//       println(f"y/x < 2^${linearThreshold}. return y/x")

      // if y/x is enough small (ex < linearThreshold), atan(y/x) ~ y/x
      (y_over_x_ex, SafeLong(y_over_x_manW1.toLong))

    } else {
//       println(f"y/x <= 1. return atan(y/x) by using tables")

      // otherwise, use table for ex is in [linearThreshold, -1].
      // in FP32, ex is in [-12, -1].

      val exadr = (-y_over_x_ex - 1).toInt
      val t = ts(exadr)

      assert(adrW   == t.adrW)
      assert(nOrder == t.nOrder)
      assert(bp     == t.bp)

      if(adrW >= manW && nOrder != 0) {
        println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
      }
      val order = if(adrW >= manW) { 0 } else { nOrder }

      val (zEx, zMan0) = if (order == 0) {
        val adr   = y_over_x_man.toInt
        val res0  = t.interval(adr).eval(0L, 0)

        val scaling = (-y_over_x_ex - 1).toInt
        val shift = calcW+1 - res0.toLong.toBinaryString.length
        val res   = res0 << shift

        (-shift - scaling, res)
      } else {
        // TODO: here we use negative coefficient in x^2 term.

        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, y_over_x_man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, y_over_x_man).toInt
        val res0 = t.interval(adr).eval(d.toLong, dxbp)

        // the result is multiplied by 2^(-ex-1).
        // if x < 1, atan(x) = x - x^3/3 + x^5/5 + O(x^7) is less than x. So
        // in case of x is in [2^-8, 2^-7), atan(x) < 2^-7. So we can multiply
        // atan(x) by 2^7 = 2^(-(-8)-1) to contain more bits. By doing this, we
        // can keep the same precision in all 2^-8 ~ 2^0 regions.
        //   We need to correct the exponent term by this scaling term.

        val scaling = (-y_over_x_ex - 1).toInt
        val shift   = calcW+1 - res0.toLong.toBinaryString.length
        val res     = res0 << shift // normalize

        (-shift - scaling, res)
      }
      val zMan = if (extraBits > 0) {(zMan0 >> extraBits) + bit(extraBits-1, zMan0)} else {zMan0}

      (zEx, SafeLong(zMan.toLong))
    }
//     println(f"sim: atanyx      = ${atanyx_sgn}|${atanyx_ex}|${(atanyx_manW1 - (1<<manW)).toLong.toBinaryString}")
//     println(f"sim: atan(|y/x|) = ${new RealGeneric(x.spec, atanyx_sgn, atanyx_ex + exBias, atanyx_manW1 - (1<<manW)).toDouble}")
//     println(f"ref: atan(|y/x|) = ${atan(y.toDouble / x.toDouble)}")

    // ==================================================================
    // atan2(y, x)
    // = ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|
    //   ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y|
    //   ysgn * (pi/2-atan(|x|/|y|)) .. x>0, |x|<|y|
    //   ysgn * (pi/2+atan(|x|/|y|)) .. x<0, |x|<|y|

    val zsgn = ysgn

    if(swapped) {
      val halfpi = new RealGeneric(x.spec, Pi * 0.5)
      val halfpi_manW1sft3 = (halfpi.man + (1 << manW)) << 3
      val atan_shift       = -atanyx_ex
      val atan_aligned     = (atanyx_manW1 << 3) >> atan_shift

      if(xsgn == 0) {
        val zman0 = halfpi_manW1sft3 - atan_aligned
        // z is in [0.78 ~ 1.57)

        val zman0LessThan1 = if(bit(manW+3, zman0) == 0) { 1 } else { 0 }
        val zmanRound = roundBySpec(RoundSpec.round, 3-zman0LessThan1, zman0)
        val zman = zmanRound - (1<<manW)
        val zex  = exBias - zman0LessThan1

        return new RealGeneric(x.spec, zsgn, zex, zman)
      } else {
        val zman0 = halfpi_manW1sft3 + atan_aligned
        // z is in [1.57 ~ 2.35), ex is 0 or 1
        val zman0MoreThan2 = bit(manW+1+3, zman0)
        val zmanRound = roundBySpec(RoundSpec.round, 3+zman0MoreThan2, zman0)
        val zman = zmanRound - (1<<manW)
        val zex  = zman0MoreThan2 + exBias

        return new RealGeneric(x.spec, zsgn, zex, zman)
      }
    } else {
      if(xsgn == 0) {
        // ysgn *   atan(|y|/|x|)      .. x>0, |x|>|y|
        val atanyx_man = atanyx_manW1 - (1<<manW)
        return new RealGeneric(x.spec, zsgn, atanyx_ex + exBias, atanyx_man)
      } else {
        // ysgn * (-atan(|y|/|x|)+pi)  .. x<0, |x|>|y|

        val pi = new RealGeneric(x.spec, Pi)

        // pi.ex = 1
        val pi_manW1   = (pi.man + (1 << manW)) << 3 // pi = 2^1 * 1.57...
        val atan_shift = -atanyx_ex
        val atan_man   = (atanyx_manW1 << 2) >> atan_shift
        val zman0      = pi_manW1 - atan_man // pi ~ pi/2.
        val zman0LessThan2 = if(bit(manW+3, zman0) == 0) { 1 } else { 0 }
        val zman = roundBySpec(RoundSpec.round, 3-zman0LessThan2, zman0)

        return new RealGeneric(x.spec, zsgn, 1-zman0LessThan2 + exBias, zman - (1<<manW))
      }
    }
  }

  def reciprocalTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    val tableD = if (order==0) {
      new FuncTableDouble( x => 1.0/(1.0+x), order )
    } else {
      val eps=pow(2.0,-manW)
      new FuncTableDouble( x => 1.0/(2.0 - (x+eps)), order )
    }
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  def atanTableGeneration( order : Int, adrW : Int, manW : Int, fracW : Int ) = {
    val linearThreshold = -math.round(math.ceil(manW / 2.0 + 1.0)).toInt

    val nOrder = if (adrW >= manW) { 0 } else { order }

    val maxCalcWidth = (-1 to linearThreshold by -1).map(ex => {
        val tableD = new FuncTableDouble( x => scalb(atan(scalb(1.0 + x, ex)), -ex-1), nOrder)
        tableD.addRange(0.0, 1.0, 1<<adrW)
        val tableI = new FuncTableInt( tableD, fracW ) // convert float table into int
        tableI.calcWidth
      }).reduce( (lhs, rhs) => {
        lhs.zip(rhs).map( x => max(x._1, x._2))
      })

    (-1 to linearThreshold by -1).map( ex => {
      val tableD = new FuncTableDouble( x => scalb(atan(scalb(1.0 + x, ex)), -ex-1), nOrder)
      tableD.addRange(0.0, 1.0, 1<<adrW)
      new FuncTableIntFixedWidth( tableD, fracW, maxCalcWidth )
    })
  }

  val atan2F32ReciprocalTableI = ATan2Sim.reciprocalTableGeneration( 2, 8, 23, 23+2 )
  val atan2F32ATanTableI       = ATan2Sim.atanTableGeneration( 2, 8, 23, 23+2 )
  val atan2F32Sim              = atan2SimGeneric(atan2F32ReciprocalTableI, atan2F32ATanTableI, _, _ )

  val atan2BF16ReciprocalTableI = ATan2Sim.reciprocalTableGeneration( 0, 7, 7, 7 )
  val atan2BF16ATanTableI       = ATan2Sim.atanTableGeneration( 0, 7, 7, 7 )
  val atan2BF16Sim              = atan2SimGeneric(atan2BF16ReciprocalTableI, atan2BF16ATanTableI, _, _ )
}
