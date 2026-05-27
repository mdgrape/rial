//% @file GenRandomGaussianSim.scala
//
// Bit-exact software model of the GenRandomGaussian hardware generator. Given
// the same rndW-bit random input the HW consumes, it returns the identical
// output float bits. Pipeline stages do not change the functional result, so
// the model is a pure function of the input bits (and the table parameters).
//
package rial.rng

import rial.arith.{RealSpec, RealGeneric}
import rial.table._
import spire.math.SafeLong

object GenRandomGaussianSim {

  // Building a ProbitTable runs erfInv over every interval, so cache it per
  // (spec, nOrder, rndW).
  private val cache =
    scala.collection.mutable.Map.empty[(RealSpec, Int, Int),
      (ProbitTable, PolynomialBitWidth, PolynomialCoefficients)]

  private def tableFor(spec: RealSpec, nOrder: Int, rndW: Int) =
    cache.getOrElseUpdate((spec, nOrder, rndW), {
      val manW       = spec.manW
      val xex0Width  = rndW - 2 - manW
      val nIntervals = xex0Width + 1
      val table  = new ProbitTable(spec, nOrder, nIntervals, manW, manW)
      val bitW   = table.getBitWidth()
      val coeffs = table.generateTable(bitW)
      (table, bitW, coeffs)
    })

  def genGaussian(rnd: BigInt, rndW: Int, spec: RealSpec, nOrder: Int): RealGeneric = {
    val manW   = spec.manW
    val exW    = spec.exW
    val exBias = spec.exBias

    val (table, _, coeffs) = tableFor(spec, nOrder, rndW)
    val xex0Width  = rndW - 2 - manW
    val nIntervals = xex0Width + 1

    val zsgn = (rnd >> (rndW - 1)) & 1
    val page = (rnd >> (rndW - 2)) & 1
    val xex0 = (rnd >> manW) & ((BigInt(1) << xex0Width) - 1)
    val xman = rnd & ((BigInt(1) << manW) - 1)

    val xex = if (xex0 == 0) BigInt(xex0Width) else BigInt(xex0.lowestSetBit)

    val addrRaw = if (page == 0) BigInt(0) else xex + 1
    val addr    = (if (addrRaw >= nIntervals) BigInt(nIntervals - 1) else addrRaw).toInt

    // polynomial dx (unsigned; PolynomialEval.eval converts to signed internally)
    val dx =
      if (nOrder == 0) BigInt(0)
      else if (page == 0) ((BigInt(1) << (manW - 1)) | (xman >> 1)) >> xex.toInt
      else xman

    // polynomial output, clamped to outBitW like the HW PolynomialEval
    val resMax  = (BigInt(1) << manW) - 1
    val polyRaw = PolynomialEval.eval(addr, dx, coeffs)
    val polyOut = if (polyRaw < 0 || polyRaw > resMax) resMax else polyRaw

    val yManW1 = (BigInt(1) << manW) | polyOut

    val (mulEx, mulMan) =
      if (page == 0) {
        (BigInt(exBias - 2) - xex, xman)
      } else {
        val scale = table.scaleTable(addr)
        ((scale >> manW) & ((BigInt(1) << exW) - 1), scale & ((BigInt(1) << manW) - 1))
      }
    val mulManW1 = (BigInt(1) << manW) | mulMan

    val prodW         = (1 + manW) * 2
    val prod          = yManW1 * mulManW1
    val prodMoreThan2 = (prod >> (prodW - 1)) & 1
    val prodManW1 =
      if (prodMoreThan2 == 1) (prod >> (prodW - 1 - manW)) & ((BigInt(1) << (manW + 1)) - 1)
      else                    (prod >> (prodW - 2 - manW)) & ((BigInt(1) << (manW + 1)) - 1)

    val zman = prodManW1 & ((BigInt(1) << manW) - 1)
    val zex  = mulEx + prodMoreThan2

    val bits = (zsgn << (exW + manW)) | ((zex & ((BigInt(1) << exW) - 1)) << manW) | zman
    new RealGeneric(spec, SafeLong(bits))
  }
}
