//% @file atan2Sim.scala
//
// Simulators for atan2(y, x) stage 1, calculating min(x,y)/max(x,y)
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.math

import scala.language.reflectiveCalls
import scala.math._
import java.lang.Math.scalb

import chisel3.util.log2Up

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


//
// atan2-1 encodes...
// - ysgn into sgn.
// - state (x>0 && |x|>|y| etc, 2bit) into the MSB of the exponent and the LSB of the mantissa
//
// the MSB of the exponent will not be used if it is not a special value (then |x/y| < 1).
// the LSB of the mantissa will be zeroed in phase2.
//

object ATan2Phase1Sim {

  // (z, ATan2Status, ATan2SpecialValue, ysgn)
  def atan2Phase1SimGeneric(
    t_rec : FuncTableInt, y : RealGeneric, x : RealGeneric
  ): RealGeneric = {

//     println("==================================================")
    assert(x.spec == y.spec)

    val spec = x.spec

    val exW    = spec.exW
    val manW   = spec.manW
    val exBias = spec.exBias

    val fracW     = t_rec.bp
    val extraBits = fracW - manW

    val xsgn   = x.sgn
    val ysgn   = y.sgn

    val xnan  = x.isNaN
    val xinf  = x.isInfinite
    val xzero = x.isZero

    val ynan  = y.isNaN
    val yinf  = y.isInfinite
    val yzero = y.isZero

//     println(f"x = ${x.sgn}|${x.ex}(${x.ex-exBias})|${x.man.toLong.toBinaryString}")
//     println(f"y = ${y.sgn}|${y.ex}(${y.ex-exBias})|${y.man.toLong.toBinaryString}")
//     println(f"ATan2Phase1Sim: x = ${x.toDouble}, y = ${y.toDouble}, atan2(y, x) = ${atan2(y.toDouble, x.toDouble)}")

    val xySameMan = x.man == y.man

    val yIsLarger = slice(0, spec.W-1, x.value) < slice(0, spec.W-1, y.value)

    val minxy = if (yIsLarger) { x } else { y }
    val maxxy = if (yIsLarger) { y } else { x }

    val minex = minxy.ex
    val maxex = maxxy.ex

    val diffexDec = if(yIsLarger) {
      if(y.man > x.man) {1} else {0}
    } else {
      if(x.man > y.man) {1} else {0}
    }

    val zex0 = if (maxex == maskI(exW) && minex == maskI(exW)) {
      exBias
    } else if((maxex == maskI(exW) && minex != maskI(exW)) || minex == 0) {
      0
    } else {
      minex - maxex - diffexDec + exBias
    }

    val zEx  = if(zex0 < 0) {0} else {zex0}
    val zeroed = (zEx == 0)

    // println(f"x = ${xsgn}|${x.ex}|${x.man.toLong.toBinaryString}")
    // println(f"y = ${ysgn}|${y.ex}|${y.man.toLong.toBinaryString}")
    // println(f"minex = ${minex}")
    // println(f"minex = ${maxex}")
    // println(f"exdec = ${diffexDec}")
    // println(f"zex0  = ${zex0}")
    // println(f"zex   = ${zEx}")

    val tooLargeX = (zeroed && !yIsLarger) || ( xinf && !yinf) || (!xzero &&  yzero)
    val tooLargeY = (zeroed &&  yIsLarger) || (!xinf &&  yinf) || ( xzero && !yzero)

//     println(f"ATan2Phase1Sim: yIsLarger = ${yIsLarger}")

    // ------------------------------------------------------------------------

    val atan2Status = (if(yIsLarger) {2+x.sgn} else {x.sgn})
    val atan2Status0 = (atan2Status >> 0) & 1
    val atan2Status1 = (atan2Status >> 1) & 1

    // ------------------------------------------------------------------------
    // special values (nan boxing)

    val xpos = x.sgn == 0
    val xneg = x.sgn == 1
    val znan      =  (xnan ||  ynan) || (xzero && yzero)
    val zzero     = tooLargeX && xpos
    val zpi       = tooLargeX && xneg
    val zhalfpi   = tooLargeY
    val z1piover4 =  (xinf && yinf && xpos) || (x.ex == y.ex && x.man == y.man && xpos)
    val z3piover4 =  (xinf && yinf && xneg) || (x.ex == y.ex && x.man == y.man && xneg)

    if(znan) {
      return new RealGeneric( spec, ysgn, maskI(spec.exW), (SafeLong(1) << (spec.manW-1)) )
    } else if(zzero) {
      return new RealGeneric( spec, ysgn, maskI(spec.exW), 1 )
    } else if(zpi) {
      return new RealGeneric( spec, ysgn, maskI(spec.exW), 2 )
    } else if(zhalfpi) {
      return new RealGeneric( spec, ysgn, maskI(spec.exW), 3 )
    } else if(z1piover4) {
      return new RealGeneric( spec, ysgn, maskI(spec.exW), 4 )
    } else if(z3piover4) {
      return new RealGeneric( spec, ysgn, maskI(spec.exW), 5 )
    }

    // ------------------------------------------------------------------------
    // reciprocal table

    val recMan = if (t_rec.nOrder==0) {
      val adr = maxxy.man.toInt

      val res = if (maxxy.man == 0) {
        0L
      } else {
        t_rec.interval(adr).eval(0L, 0).toLong
      }
//       println(f"sim: atan2-1: rec adr = ${adr.toBinaryString}")
//       println(f"sim: atan2-1: rec res = ${res.toBinaryString}")
      res
    } else {
      val adrW = t_rec.adrW
      val dxbp = manW-adrW-1
      val d    = (SafeLong(1)<<dxbp) - slice(0, manW-adrW, maxxy.man)-1
      val adr  = maskI(adrW)-slice(manW-adrW, adrW, maxxy.man).toInt
      if (maxxy.man == 0) {
        0L
      } else {
        t_rec.interval(adr).eval(d.toLong, dxbp).toLong
      }
    }

//     println(f"sim: zres = ${recMan.toLong.toBinaryString}")
//     println(f"sim: zEx0 = ${zEx0.toLong.toBinaryString}")

    // ------------------------------------------------------------------------
    // atan2 stage1 postprocess (minxy * rec(maxxy))

    val denomW1 = (SafeLong(1)<<fracW) + recMan
    val numerW1 = (SafeLong(1)<<manW) + minxy.man

//     println(f"sim: denomW1 = ${denomW1.toLong.toBinaryString}")
//     println(f"sim: numerW1 = ${numerW1.toLong.toBinaryString}")

    val zProd = denomW1 * numerW1
    val bp    = fracW + manW

    val zProdMoreThan2 = bit(fracW+1 + manW+1 - 1, zProd)
    val roundBits = bp - manW + zProdMoreThan2
    val zProdRounded = roundBySpec(RoundSpec.roundToEven, roundBits, zProd)
    val zProdMoreThan2AfterRound = bit(manW+1, zProdRounded)
    val zMan0 = zProdRounded >> zProdMoreThan2AfterRound

    val zSgn = 0
    val zMan = if(zEx == 0 || xySameMan) {
      SafeLong(0)
    } else if(maxxy.man == 0) {
      minxy.man
    } else {
      zMan0
    }

    // -----------------------------------------------------------------------
    // encode: atan2Status(2bits -> sgn and ex.msb), ysgn (lsb of mantissa)

    assert( zEx <= spec.exBias )

    val zSgnEnc = ysgn
    val zExEnc  = zEx | (atan2Status1 << (spec.exW-1))
    val zManEnc = (zMan & ((SafeLong(1) << spec.manW) - 2)) + atan2Status0

    val z = new RealGeneric(x.spec, zSgnEnc, zExEnc, zManEnc)

//     println(f"sim: atan2-1: z = (${zSgn}|${zEx}(${zEx-exBias})|${zMan.toBinaryString}) (${z.toDouble}) should be (${min(x.toDouble, y.toDouble) / max(x.toDouble, y.toDouble)})")
    z
  }
}
