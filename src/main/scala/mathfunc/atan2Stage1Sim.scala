//% @file atan2Sim.scala
//
// Simulators for atan2(y, x) stage 1, calculating min(x,y)/max(x,y)
// Copyright (C) Toru Niina RIKEN BDR 2021
//
package rial.mathfunc

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

object ATan2Stage1Sim {

  def atan2Stage1SimGeneric( t_rec : FuncTableInt, y : RealGeneric, x : RealGeneric ):
    (RealGeneric, Int, Int, Int) = { // (z, ATan2Status, ATan2SpecialValue, ysgn)

//     println("==================================================")
    assert(x.spec == y.spec)

    val exW    = x.spec.exW
    val manW   = x.spec.manW
    val exBias = x.spec.exBias

    val fracW     = t_rec.bp
    val extraBits = fracW - manW

    val xsgn   = x.sgn
    val ysgn   = y.sgn

//     println(f"x = ${x.sgn}|${x.ex}(${xex})|${xman.toLong.toBinaryString}")
//     println(f"y = ${y.sgn}|${y.ex}(${yex})|${yman.toLong.toBinaryString}")
//     println(f"ATan2Stage1Sim: x = ${x.toDouble}, y = ${y.toDouble}, atan2(y, x) = ${atan2(y.toDouble, x.toDouble)}")

    val yIsLarger = slice(0, x.spec.W-1, x.value) < slice(0, x.spec.W-1, y.value)

//     println(f"ATan2Stage1Sim: yIsLarger = ${yIsLarger}")

    val xnan  = x.isNaN
    val xinf  = x.isInfinite
    val xzero = x.isZero

    val ynan  = y.isNaN
    val yinf  = y.isInfinite
    val yzero = y.isZero

    val xpos = x.sgn == 0
    val xneg = x.sgn == 1
    val znan      =  (xnan ||  ynan) || ( xzero &&  yzero)
    val zzero     = ((xinf && !yinf) || (!xzero &&  yzero)) && xpos
    val zpi       = ((xinf && !yinf) || (!xzero &&  yzero)) && xneg
    val zhalfpi   = (!xinf &&  yinf) || ( xzero && !yzero)
    val z1piover4 =  (xinf &&  yinf) &&  xpos
    val z3piover4 =  (xinf &&  yinf) &&  xneg

    val atan2Status       = (if(yIsLarger) {2+x.sgn} else {x.sgn})
    val atan2SpecialValue = if (znan)      { 1 }
                       else if (zzero)     { 2 }
                       else if (zpi)       { 3 }
                       else if (zhalfpi)   { 4 }
                       else if (z1piover4) { 5 }
                       else if (z3piover4) { 6 }
                       else                { 0 }

    val minxy = if (yIsLarger) { x } else { y }
    val maxxy = if (yIsLarger) { y } else { x }

    // ------------------------------------------------------------------------
    // reciprocal table

    val recMan = if (t_rec.nOrder==0) {
      val adr = maxxy.man.toInt
      if (maxxy.man==0) {
        0
      } else {
        t_rec.interval(adr).eval(0L, 0)
      }
    } else {
      val adrW = t_rec.adrW
      val dxbp = manW-adrW-1
      val d    = (SafeLong(1)<<dxbp) - slice(0, manW-adrW, maxxy.man)-1
      val adr  = maskI(adrW)-slice(manW-adrW, adrW, maxxy.man).toInt
      if (maxxy.man==0) {
        0
      } else {
        t_rec.interval(adr).eval(d.toLong, dxbp)
      }
    }

    val xySameMan = x.man == y.man
    val yOverXEx  = minxy.ex + exBias - 1 - maxxy.ex + (if(maxxy.man == 0) {1} else {0})
//     println(f"sim: yOverXEx = ${yOverXEx.toLong.toBinaryString}")

    // ------------------------------------------------------------------------
    // atan2 stage1 postprocess (minxy * rec(maxxy))

    val denomW1 = (1<<fracW) + recMan
    val numerW1 = (1<<manW) + minxy.man

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
    val zEx0 = yOverXEx + zProdMoreThan2 + zProdMoreThan2AfterRound + (if(xySameMan) {1} else {0})
//     val zEx0 = minxy.ex - maxxy.ex - 1 + exBias + zProdMoreThan2 + zProdMoreThan2AfterRound
    val zEx = if(zEx0 < 0) {0} else if((1<<exW) <= zEx0) {maskI(exW)} else {zEx0}

    val zMan = if(zEx == 0 || xySameMan) {0L} else {zMan0.toLong}

    (new RealGeneric(x.spec, zSgn, zEx, zMan), atan2Status, atan2SpecialValue, ysgn)
  }
}
