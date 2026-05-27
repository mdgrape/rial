//% @file GenRandomGaussianSim.scala
//
// Bit-exact software model of the GenRandomGaussian hardware generator. Given
// the same rndW-bit random input the HW consumes, it returns the identical
// output float bits. Pipeline stages do not change the functional result, so
// the model is a pure function of the input bits (and the table parameters).
//
// Mirrors GenRandomGaussian (see GenGaussian.scala) exactly:
//   - page 0: central region.
//   - page 1, xex0 != 0: main octave, selected by the lowest set bit of xex0.
//   - page 1, xex0 == 0: deep tail, octave selected by xman's leading-zero
//     count (clz); the deepest `lumpOctaves` octaves merge into one capped
//     interval. dx is the normalized (leading-zeros + leading-1 removed) xman.
//   - the polynomial output clamps a negative result to 0 (clampNegToZero).
//
package rial.rng

import rial.arith.{RealSpec, RealGeneric}
import rial.table._
import spire.math.SafeLong

object GenRandomGaussianSim {

  // Building a ProbitTable runs erfc inversion over every interval, so cache it
  // per (spec, nOrder, rndW, lumpOctaves).
  private val cache =
    scala.collection.mutable.Map.empty[(RealSpec, Int, Int, Int),
      (ProbitTable, PolynomialBitWidth, PolynomialCoefficients)]

  private def tableFor(spec: RealSpec, nOrder: Int, rndW: Int, lumpOctaves: Int) =
    cache.getOrElseUpdate((spec, nOrder, rndW, lumpOctaves), {
      val table  = new ProbitTable(spec, nOrder, nIntervalsOf(rndW, spec, lumpOctaves),
                                   spec.manW, spec.manW, lumpOctaves)
      val bitW   = table.getBitWidth()
      val coeffs = table.generateTable(bitW)
      (table, bitW, coeffs)
    })

  private def mainOctavesOf(rndW: Int, spec: RealSpec): Int = rndW - 2 - spec.manW
  private def nIntervalsOf(rndW: Int, spec: RealSpec, lumpOctaves: Int): Int =
    1 + (mainOctavesOf(rndW, spec) + (spec.manW - lumpOctaves)) + 1

  def genGaussian(rnd: BigInt, rndW: Int, spec: RealSpec, nOrder: Int,
                  lumpOctaves: Int = 2): RealGeneric = {
    val manW   = spec.manW
    val exW    = spec.exW
    val exBias = spec.exBias

    val (table, _, coeffs) = tableFor(spec, nOrder, rndW, lumpOctaves)

    val mainOctaves = mainOctavesOf(rndW, spec)        // e.g. 39
    val deepIndiv   = manW - lumpOctaves               // e.g. 21
    val nIntervals  = nIntervalsOf(rndW, spec, lumpOctaves) // e.g. 62
    val xex0Width   = mainOctaves
    val manMask     = (BigInt(1) << manW) - 1

    val zsgn = (rnd >> (rndW - 1)) & 1
    val page = (rnd >> (rndW - 2)) & 1
    val xex0 = (rnd >> manW) & ((BigInt(1) << xex0Width) - 1)
    val xman = rnd & manMask

    // lowest set bit of xex0 (main range), and leading-zero count of xman (deep)
    val xex = if (xex0 == 0) xex0Width else xex0.lowestSetBit
    val clz = if (xman == 0) manW else manW - xman.bitLength

    val addr: Int =
      if (page == 0)            0
      else if (xex0 != 0)       xex + 1
      else if (clz < deepIndiv) mainOctaves + 1 + clz
      else                      nIntervals - 1 // lumped deep interval

    // polynomial input dx (unsigned; PolynomialEval.eval flips the MSB to signed)
    val dx =
      if (nOrder == 0)          BigInt(0)
      else if (page == 0)       ((BigInt(1) << (manW - 1)) | (xman >> 1)) >> xex
      else if (xex0 != 0)       xman
      else if (clz < deepIndiv) (xman << (clz + 1)) & manMask
      else                      (xman << (manW - lumpOctaves)) & manMask

    // polynomial output, clamped like PolynomialEval with clampNegToZero = true:
    // a negative undershoot saturates to 0, an overshoot to resMax.
    val resMax  = manMask
    val polyRaw = PolynomialEval.eval(addr, dx, coeffs)
    val polyOut = if (polyRaw < 0) BigInt(0) else if (polyRaw > resMax) resMax else polyRaw

    val yManW1 = (BigInt(1) << manW) | polyOut

    val (mulEx, mulMan) =
      if (page == 0) {
        (BigInt(exBias - 2) - xex, xman)
      } else {
        val scale = table.scaleTable(addr)
        ((scale >> manW) & ((BigInt(1) << exW) - 1), scale & manMask)
      }
    val mulManW1 = (BigInt(1) << manW) | mulMan

    val prodW         = (1 + manW) * 2
    val prod          = yManW1 * mulManW1
    val prodMoreThan2 = (prod >> (prodW - 1)) & 1
    val prodManW1 =
      if (prodMoreThan2 == 1) (prod >> (prodW - 1 - manW)) & ((BigInt(1) << (manW + 1)) - 1)
      else                    (prod >> (prodW - 2 - manW)) & ((BigInt(1) << (manW + 1)) - 1)

    val zman = prodManW1 & manMask
    val zex  = mulEx + prodMoreThan2

    val bits = (zsgn << (exW + manW)) | ((zex & ((BigInt(1) << exW) - 1)) << manW) | zman
    new RealGeneric(spec, SafeLong(bits))
  }
}
