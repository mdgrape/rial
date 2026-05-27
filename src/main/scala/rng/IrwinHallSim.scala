//% @file IrwinHallSim.scala
//
// Bit-exact software model of the IrwinHall hardware generator (IrwinHall.scala).
// Given the same 12 rndW-bit uniforms the HW consumes, it returns the identical
// output float bits. The tree-sum + FP-conversion pipeline does not change the
// functional result, so the model is a pure function of the inputs.
//
package rial.rng

import rial.arith.{RealSpec, RealGeneric}
import spire.math.SafeLong

object IrwinHallSim {

  /** Bit-exact model of [[IrwinHall]]: sum 12 rndW-bit uniforms (mean = 6*2^rndW,
   *  ~N(6,1) scaled) and normalize the centered sum into a float. */
  def irwinHall(rnds: Seq[BigInt], rndW: Int, spec: RealSpec): RealGeneric = {
    require(rnds.length == 12, s"IrwinHall needs 12 inputs, got ${rnds.length}")
    val manW   = spec.manW
    val exW    = spec.exW
    val exBias = spec.exBias

    val rndMask = (BigInt(1) << rndW) - 1
    val sum     = rnds.map(_ & rndMask).sum        // in [0, 12*2^rndW), fits in rndW+4 bits

    val sumW = rndW + 4
    val manWp = sumW - 1                            // mantissa width after dropping the leading 1
    require(manWp >= manW + 1, s"rndW=$rndW too small for manW=$manW (need rndW >= manW-2)")

    val six   = BigInt(6) << rndW                   // the mean, 6 in the top 4 bits
    val zZero = sum == six
    val zSgn  = if (sum >= six) BigInt(0) else BigInt(1)
    val zAbs  = if (sum >= six) sum - six else six - sum

    val bits =
      if (zZero) {
        zSgn << (exW + manW)                        // ex = 0, man = 0
      } else {
        // leading-zero count of zAbs within sumW bits (HW: PriorityEncoder(Reverse))
        val zShift    = sumW - zAbs.bitLength
        val zShifted0 = (zAbs << zShift) & ((BigInt(1) << sumW) - 1) // leading 1 at bit sumW-1
        val zShifted  = zShifted0 & ((BigInt(1) << manWp) - 1)       // drop the leading 1

        val topManW  = zShifted >> (manWp - manW)                   // top manW bits
        val roundBit = (zShifted >> (manWp - manW - 1)) & 1         // round-half-up
        val zRounded = topManW + roundBit

        val zMoreThan2 = (zRounded >> manW) & 1                     // rounding carry
        val zMan       = zRounded & ((BigInt(1) << manW) - 1)
        val zEx        = BigInt(exBias + 3) - zShift + zMoreThan2

        (zSgn << (exW + manW)) | ((zEx & ((BigInt(1) << exW) - 1)) << manW) | zMan
      }

    new RealGeneric(spec, SafeLong(bits))
  }
}
