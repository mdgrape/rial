//% @file GenRandomFloatSim.scala
//
// Bit-exact software models of the GenRandomFloat* hardware generators. Each
// takes the same random input bits the HW consumes and returns the identical
// output float bits.
//
package rial.rng

import rial.arith.{RealSpec, RealGeneric}
import spire.math.SafeLong

object GenRandomFloatSim {

  /** Bit-exact model of [[GenRandomFloat12]]: a float in [1, 2) whose mantissa
   *  is the manW-bit random input directly (sign 0, exponent = bias). */
  def genFloat12(rnd: BigInt, spec: RealSpec): RealGeneric = {
    val man  = rnd & ((BigInt(1) << spec.manW) - 1)
    val bits = (BigInt(spec.exBias) << spec.manW) | man
    new RealGeneric(spec, SafeLong(bits))
  }

  /** Bit-exact model of [[GenRandomFloat01]]: interpret the rndW-bit input as a
   *  fixed-point fraction rnd / 2^rndW in [0, 1); the leading-zero count sets the
   *  exponent and the bits below the leading 1 form the mantissa. */
  def genFloat01(rnd0: BigInt, rndW: Int, spec: RealSpec): RealGeneric = {
    val manW   = spec.manW
    val exBias = spec.exBias

    val rnd = rnd0 & ((BigInt(1) << rndW) - 1)

    if (rnd == 0) {
      new RealGeneric(spec, SafeLong(BigInt(0))) // sgn=0, ex=0, man=0
    } else {
      // PriorityEncoder(Reverse(rnd)): number of leading zeros within rndW bits.
      val zerosAtMSB = rndW - rnd.bitLength
      val shifted    = rnd << zerosAtMSB          // leading 1 now at bit rndW-1
      // manW1 = shifted[rndW-1 : rndW-1-manW]; the MSB is the leading 1.
      val manW1 = (shifted >> (rndW - 1 - manW)) & ((BigInt(1) << (manW + 1)) - 1)
      val man   = manW1 & ((BigInt(1) << manW) - 1)         // drop the leading 1
      val ex    = BigInt(exBias - 1 - zerosAtMSB)
      val bits  = (ex << manW) | man                        // sgn = 0
      new RealGeneric(spec, SafeLong(bits))
    }
  }
}
