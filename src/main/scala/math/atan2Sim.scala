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
  // atan(1/x) = pi/2 - atan(x)
  //
  // if x > 2^manW:
  //   atan(x) = pi/2 - atan(1/x) ~ pi/2 - 2^-manW ~ pi/2
  // if x > 2^(-linearThreshold):
  //   atan(x) = pi/2 - atan(1/x) ~ pi/2 - 1/x
  //
  // we need table for 2^-8 ~ 2^8 range for FP32
  //
  def atan2SimGeneric( t_rec : FuncTableInt, ts : Seq[FuncTableIntFixedWidth], y : RealGeneric, x : RealGeneric ) : RealGeneric = {
    assert(x.spec == y.spec)

    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val xsgn    = x.sgn
    val xexRaw  = x.ex
    val xex     = xexRaw-exBias
    val xman    = x.man

    val ysgn    = y.sgn
    val yexRaw  = y.ex
    val yex     = yexRaw-exBias
    val yman    = y.man

    val adrW      = ts(0).adrW
    val nOrder    = ts(0).nOrder
    val bp        = ts(0).bp
    val extraBits = bp - manW
    val calcW     = manW + extraBits

    println(f"x = ${x.sgn}|${x.ex}(${xex})|${xman.toLong.toBinaryString}")
    println(f"y = ${y.sgn}|${y.ex}(${yex})|${yman.toLong.toBinaryString}")
    println(f"x = ${x.toDouble}, y = ${y.toDouble}, atan2(y, x) = ${atan2(y.toDouble, x.toDouble)}")

    val pi_1over4 = new RealGeneric(x.spec, Pi * 0.25)
    val pi_over_2 = new RealGeneric(x.spec, Pi * 0.50)
    val pi_3over4 = new RealGeneric(x.spec, Pi * 0.75)
    val pi        = new RealGeneric(x.spec, Pi)

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

    val y_over_x_sgn = xsgn ^ ysgn
    val y_over_x_ex0 = yex - xex

    println(f"y_over_x_sgn = ${y_over_x_sgn}")
    println(f"y_over_x_ex0 = ${y_over_x_ex0}")

    // 1+manW+1 width (x2)
    val one_over_x_manW1 = if (t_rec.nOrder == 0) {
      val adr  = x.man
      val res0 = t_rec.interval(adr.toInt).eval(0L, 0) // in (0.5, 1]
      if(bit(manW, res0) == 1) {
        1 << (manW+1) // just 1
      } else {
        res0 // 2^-1 * 1.m [0.5, 1)
      }
    } else {
      val dxbp = manW - t_rec.adrW - 1
      val d    = (SafeLong(1) << dxbp) - slice(0, manW - t_rec.adrW, x.man)-1
      val adr  = (1 << t_rec.adrW) - 1 - slice(manW - t_rec.adrW, t_rec.adrW, x.man).toInt
      val res0 = t_rec.interval(adr.toInt).eval(d.toLong, dxbp)

      println(f"xman = ${((1<<manW) + xman).toDouble / (1<<manW)}")
      println(f"1/xm = ${(1<<manW) / ((1<<manW) + xman).toDouble}")
      println(f"res0 = ${res0.toDouble / (1<<t_rec.bp)}")

      if(x.man == 0) {
        1 << (manW+1) // just 1
      } else {
        res0 // 2^-1 * 1.m [0.5, 1)
      }
    }
    // 1.m / 1.m is in (0.5, 2) range. required shift is less than 2.
    // here we already multiplied the denominator(x) by 2, so the resulting
    // range becomes (0.25, 1). FIXME ?
    val y_over_x_man0 = one_over_x_manW1 * ((1<<manW) + y.man)
    val y_over_x_man0_clz = (1+manW+1) + (1+manW) - y_over_x_man0.toLong.toBinaryString.length
    assert(y_over_x_man0_clz <= 2)

    val y_over_x_manW1= roundBySpec(RoundSpec.roundToEven, (1+manW+1) - y_over_x_man0_clz, y_over_x_man0)
    val y_over_x_man  = y_over_x_manW1 - (1 << manW)
    val y_over_x_ex   = y_over_x_ex0 - y_over_x_man0_clz
    println(f"y_over_x_ex  = ${y_over_x_ex}")

    println(f"y_over_x = ${y.toDouble/x.toDouble} = ${new RealGeneric(x.spec, y_over_x_sgn, y_over_x_ex+exBias, y_over_x_man).toDouble}")

    // ------------------------------------------------------------------
    // x / y

    val x_over_y_sgn = xsgn ^ ysgn
    val x_over_y_ex0 = xex - yex

    val one_over_y_manW1 = if (t_rec.nOrder == 0) {
      val adr  = y.man
      val res0 = t_rec.interval(adr.toInt).eval(0L, 0) // in (0.5, 1]
      if(bit(manW, res0) == 1) {
        1 << (manW+1) // just 1
      } else {
        res0 // 2^-1 * 1.m [0.5, 1)
      }
    } else {
      val dxbp = manW - t_rec.adrW - 1
      val d    = (SafeLong(1) << dxbp) - slice(0, manW - t_rec.adrW, y.man)-1
      val adr  = (1 << t_rec.adrW) - 1 - slice(manW - t_rec.adrW, t_rec.adrW, y.man).toInt
      val res0 = t_rec.interval(adr.toInt).eval(d.toLong, dxbp)

      println(f"yman = ${((1<<manW) + yman).toDouble / (1<<manW)}")
      println(f"1/ym = ${(1<<manW) / ((1<<manW) + yman).toDouble}")
      println(f"res0 = ${res0.toDouble / (1<<t_rec.bp)}")

      if(y.man == 0) {
        1 << (manW+1) // just 1
      } else {
        res0 // 2^-1 * 1.m [0.5, 1)
      }
    }
    val x_over_y_man0 = one_over_y_manW1 * ((1<<manW) + x.man)
    val x_over_y_man0_clz = (1+manW+1) + (1+manW) - x_over_y_man0.toLong.toBinaryString.length
    assert(x_over_y_man0_clz <= 2)

    // with leading 1
    val x_over_y_manW1= roundBySpec(RoundSpec.roundToEven, (1+manW+1) - x_over_y_man0_clz, x_over_y_man0)
    val x_over_y_man  = x_over_y_manW1 - (1 << manW)
    val x_over_y_ex   = x_over_y_ex0 - x_over_y_man0_clz
    println(f"x_over_y_ex0 = ${x_over_y_ex0}")
    println(f"x_over_y_clz = ${x_over_y_man0_clz}")
    println(f"x_over_y_ex  = ${x_over_y_ex}")
    val x_over_y_exBias = if(x_over_y_ex + exBias < 0) {
      0
    } else if ((1<<exW) <= x_over_y_ex + exBias) {
      maskI(exW)
    } else {
      x_over_y_ex + exBias
    }
    println(f"x_over_y = ${x.toDouble/y.toDouble} = ${new RealGeneric(x.spec, x_over_y_sgn, x_over_y_exBias, x_over_y_man).toDouble}")

    // ==================================================================
    // calc atan(y/x) = pi/2 - atan(x/y)

    val linearThreshold = -math.round(math.ceil(manW / 2.0 + 1.0))

    val (atanyx_sgn, atanyx_ex, atanyx_manW1) = if (manW < y_over_x_ex) {
      println("2^23 < y_over_x. return pi/2")

      // if y/x is enough large (ex > manW), atan(y/x) ~ pi/2
      (y_over_x_sgn, pi_over_2.ex - exBias, SafeLong((pi_over_2.man.toLong + (1L<<manW)) << extraBits))

    } else if (-linearThreshold-1 < y_over_x_ex) { // here we can use bitwise not
      println(f"2^${-linearThreshold-1} < y_over_x. return pi/2 - x/y")

      // if y/x is reasonably large,
      // atan(y/x) = pi/2 - atan(1/(y/x)) ~ pi / 2 - 1/(y/x).
      // if y/x.ex = -linearThreshold, x/y.ex = linearThreshold-1.
      // e.g. if x is in [2^+15, 2^+16), i.e. x.ex = 15,
      //  then 1/x is in (2^-16, 2^-15], i.e. 1/x.ex = -16.
      //
      assert(x_over_y_ex < linearThreshold)

      val halfpi_man     = (pi_over_2.man + (1<<manW)) << 3
      val xy_man_shift   = -x_over_y_ex - 3
      val xy_man_aligned = (x_over_y_man + (1<<manW)) >> xy_man_shift

      // pi/2 > 1.57_(10) = 1.1001_(2). and here x_over_y_ex < linearThreshold,
      // and linearThreshold ~ (manW-1)/3. Unless manW is too small, here we
      // don't need to consider normalization. the result is already normalized.

      val res0 = halfpi_man - xy_man_aligned
      val res  = roundBySpec(RoundSpec.roundToEven, 3, res0)

      (y_over_x_sgn, 0, SafeLong(res.toLong << extraBits))

    } else if (y_over_x_ex < linearThreshold) {
      println(f"y/x < 2^${linearThreshold}. return y/x")

      // if y/x is enough small (ex < linearThreshold), atan(y/x) ~ y/x
      (y_over_x_sgn, y_over_x_ex, SafeLong(y_over_x_manW1.toLong << extraBits))

    } else if (y_over_x_ex < 0) {
      println(f"y/x < 1. return atan(y/x) by using tables")

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

      val (zEx, zMan) = if (order == 0) {
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

        println(f"res0 = ${res0.toDouble / (1<<dxbp)}")
        println(f"res  = ${res0.toLong.toBinaryString}")

        (-shift - scaling, res)
      }

      (y_over_x_sgn, zEx, SafeLong(zMan.toLong))

    } else {
      println(f"1 < y/x. return pi/2 - atan(x/y) by using tables")

      assert(0 <= y_over_x_ex && y_over_x_ex <= -linearThreshold - 1)
      assert(linearThreshold <= x_over_y_ex && x_over_y_ex <= -1)

      // use pi/2 - atan(x/y) for ex is positive.
      //
      // The table used i the case of ex = [-12, -1] can also be used in case of
      // ex = [0, 11].

      // -----------------------------------------------------------------------
      // calc atan(x/y) first, by using the same table used above.

      val exadr     = (-x_over_y_ex - 1).toInt
      val t         = ts(exadr)

      assert(adrW   == t.adrW)
      assert(nOrder == t.nOrder)
      assert(bp     == t.bp)

      if(adrW >= manW && nOrder != 0) {
        println("WARNING: table address width >= mantissa width, but polynomial order is not zero. Polynomial order is set to zero.")
      }
      val order = if(adrW >= manW) { 0 } else { nOrder }

      val (rex, rman) = if (order == 0) {
        val adr   = x_over_y_man.toInt
        val res0  = t.interval(adr).eval(0L, 0)

        val scaling = (-x_over_y_ex - 1).toInt
        val shift = calcW+1 - res0.toLong.toBinaryString.length
        val res   = res0 << shift

        (-shift - scaling, SafeLong(res))
      } else {
        // TODO: here we use negative coefficient in x^2 term.

        val dxbp = manW - adrW - 1
        val d    = slice(0, dxbp+1, x_over_y_man) - (SafeLong(1)<<dxbp)
        val adr  = slice(dxbp+1, adrW, x_over_y_man).toInt
        val res0 = t.interval(adr).eval(d.toLong, dxbp)

        // the result is multiplied by 2^(-ex-1).
        // if x < 1, atan(x) = x - x^3/3 + x^5/5 + O(x^7) is less than x. So
        // in case of x is in [2^-8, 2^-7), atan(x) < 2^-7. So we can multiply
        // atan(x) by 2^7 = 2^(-(-8)-1) to contain more bits. By doing this, we
        // can keep the same precision in all 2^-8 ~ 2^0 regions.
        //   We need to correct the exponent term by this scaling term.

        val scaling = (-x_over_y_ex - 1).toInt
        val shift   = calcW+1 - res0.toLong.toBinaryString.length
        val res     = res0 << shift // normalize

        (-shift - scaling, SafeLong(res))
      }

      // -----------------------------------------------------------------------
      // calculate atan(y/x) = pi/2 - atan(x/y). we already have atan(x/y).

      // here, always |atan(x/y)| < pi/2. we can always anchor halfpi.
      //
      // if x in [2^-1, 1), i.e. x.ex = -1, 1/x is in (1, 2], i.e. 1/x.ex = 0.
      val halfpi_man     = (pi_over_2.man + (1 << manW)) << 3
      val xy_man_shift   = -rex
      val xy_man0        = if(extraBits >= 3) {
        rman >> (extraBits - 3)
      } else {
        rman << (3 - extraBits)
      }
      val xy_man_aligned = xy_man0 >> xy_man_shift

      val zman0  = halfpi_man - xy_man_aligned
      val zshift = if(bit(manW+3, zman0) == 0) {1} else {0}

      val zex   = 0 - zshift
      val zman  = if (extraBits + zshift <= 3) {
        roundBySpec(RoundSpec.roundToEven, 3 - extraBits - zshift, zman0)
      } else {
        zman0 << (extraBits + zshift - 3)
      }
      (y_over_x_sgn, zex, SafeLong(zman.toLong))
    }

    println(f"sim: atan(|y/x|) = ${new RealGeneric(x.spec, atanyx_sgn, atanyx_ex + exBias, atanyx_manW1 - (1<<manW)).toDouble}")
    println(f"ref: atan(|y/x|) = ${atan(y.toDouble / x.toDouble)}")

    // ==================================================================
    // correct x < 0 case.
    // - if 0 < x,          then atan2(y, x) = atan(y/x).
    // - if x < 0 && 0 < y, then atan2(y, x) = atan(y/x) + pi, atan2(y,x) < 0
    // - if x < 0 && y < 0, then atan2(y, x) = atan(y/x) - pi, atan2(y,x) > 0
    //
    // -|atan(y/x)| + pi
    //  |atan(y/x)| - pi = pi - |atan(y/x)|

    if(xsgn == 0) {
      val atanyx_man = roundBySpec(RoundSpec.roundToEven, extraBits, atanyx_manW1) - (1<<manW)
      return new RealGeneric(x.spec, atanyx_sgn, atanyx_ex + exBias, atanyx_man)
    } else {
      println("x < 0. correct the result by +/- pi")
      // x < 0 && 0 < y, atan(y/x) < 0, atan2(y,x) =   pi - |atan(y/x)|
      // x < 0 && y < 0, atan(y/x) > 0, atan2(y,x) = -(pi - |atan(y/x)|)
      // since |atan(y/x)| < pi/2, the sign depends only on the sign of y.
      //
      // sim: atan(|y/x|) = 1.601934782229364E-4
      // ref: atan(|y/x|) = 1.621186825437904E-4
      // x < 0. correct the result by +/- pi
      // atanyx_ex  = -13
      // atan_shift = 13
      // atan_manW1 =  10101001111111100110011100
      // atan_man   =  101010011111111
      // pi_man     = 110010010000111111011011000
      // zman0      = 110010010000010100111011001
      // zman0LT2   = 0
      // test: x    = -0.08609642833471298
      // test: y    = -1.3957839655631687E-5
      // test: y/x  = 1.621186839640834E-4
      // test: ref  = -3.1414305349072493
      // test: sim  = -3.140944242477417
      // test: test(1|128(1)|10010010000010100111011(49053b)) != ref(1|128(1)|10010010000110100110011(490d33))
      //
      //

      val zsgn = ysgn

      println(f"atanyx_ex = ${atanyx_ex}")
      assert(atanyx_ex <= 0)

      // pi.ex = 1
      val pi_manW1   = (pi.man + (1 << manW)) << 3 // pi = 2^1 * 1.57...
      val atan_shift = -atanyx_ex
      val atan_man0  = if(extraBits >= 2) {
        atanyx_manW1 >> (extraBits - 2)
      } else {
        atanyx_manW1 << (2 - extraBits)
      }
      val atan_man   = atan_man0 >> atan_shift
      val zman0      = pi_manW1 - atan_man // pi ~ pi/2.
      val zman0LessThan2 = if(bit(manW+3, zman0) == 0) { 1 } else { 0 }
      val zman = roundBySpec(RoundSpec.roundToEven, 3-zman0LessThan2, zman0)

      println(f"atan_shift = ${atan_shift}")
      println(f"atan_manW1 = ${atanyx_manW1.toLong.toBinaryString}")
      println(f"atan_man   = ${atan_man.toLong.toBinaryString}")
      println(f"pi_man     = ${pi_manW1.toLong.toBinaryString}")
      println(f"zman0      = ${zman0.toLong.toBinaryString}")
      println(f"zman0LT2   = ${zman0LessThan2}")

      return new RealGeneric(x.spec, zsgn, 1-zman0LessThan2 + exBias, zman - (1<<manW))
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
