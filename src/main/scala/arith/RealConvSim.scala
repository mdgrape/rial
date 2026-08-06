//% @file RealConvSim.scala
//
// Conversion between the raw bit pattern of a RealSpec and float/double, for
// the JNI binding (jni/, rialcpp).
//
// The C++ side works on bit patterns (`std::uint64_t`) whenever the emulated
// format does not fit in a `float`, and needs a way to read those back as an
// ordinary number. The spec is passed per call, so this stays stateless: there
// is no table to keep and nothing to configure.
//
package rial.arith

import spire.math.SafeLong

object RealConvSim {

  private def maskOf(w: Int): BigInt = { (BigInt(1) << w) - 1 }

  private def realOf(exW: Int, exBias: Int, manW: Int, bits: Long): RealGeneric = {
    val spec = new RealSpec(exW, exBias, manW)
    new RealGeneric(spec, SafeLong(BigInt(bits) & maskOf(spec.W)))
  }

  private def bitsOf(r: RealGeneric): Long = {
    (r.value.toBigInt & maskOf(r.spec.W)).toLong
  }

  /** the value of the (exW, exBias, manW) number encoded in `bits`. */
  def toDouble(exW: Int, exBias: Int, manW: Int, bits: Long): Double = {
    realOf(exW, exBias, manW, bits).toDouble
  }
  def toFloat(exW: Int, exBias: Int, manW: Int, bits: Long): Float = {
    realOf(exW, exBias, manW, bits).toFloat
  }

  /** `x` rounded into the (exW, exBias, manW) format, as a bit pattern. */
  def fromDouble(exW: Int, exBias: Int, manW: Int, x: Double): Long = {
    bitsOf(RealGeneric.fromDouble(new RealSpec(exW, exBias, manW), x))
  }
  def fromFloat(exW: Int, exBias: Int, manW: Int, x: Float): Long = {
    fromDouble(exW, exBias, manW, x.toDouble)
  }
}
