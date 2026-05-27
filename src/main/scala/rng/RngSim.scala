//% @file RngSim.scala
//
// Stable software-simulation facade for the JNI binding (jni/, rialcpp), for
// the random-number generators. Mirrors rial.math.MathSim: the C++ side looks
// up exactly this object and these method names, so the binding does not depend
// on the internal generator/Sim object names.
//
// Every entry takes the raw random bits as an unsigned 64-bit value (`Long`
// bit-pattern) and returns the generated float, with the standard FP32 config.
//
package rial.rng

import rial.arith.RealSpec

object RngSim {
  private val spec = RealSpec.Float32Spec

  // reinterpret a (possibly negative) jlong as its unsigned 64-bit value
  private def u64(rnd: Long): BigInt = BigInt(rnd) & ((BigInt(1) << 64) - 1)

  /** uniform float in [1, 2) -- low manW bits of `rnd` become the mantissa. */
  def genFloat12(rnd: Long): Float =
    GenRandomFloatSim.genFloat12(u64(rnd), spec).toFloat

  /** uniform float in [0, 1) from a 32-bit random fraction. */
  def genFloat01(rnd: Long): Float =
    GenRandomFloatSim.genFloat01(u64(rnd), 32, spec).toFloat

  /** N(0,1) gaussian from a 64-bit random word (nOrder = 3). */
  def genGaussian(rnd: Long): Float =
    GenRandomGaussianSim.genGaussian(u64(rnd), 64, spec, 3).toFloat
}
