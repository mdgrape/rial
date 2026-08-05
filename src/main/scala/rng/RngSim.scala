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

  /** Squares32 counter-based uniform RNG: 32-bit random word for (count, key),
   *  returned as the raw bit pattern in a jint. */
  def squares32(count: Long, key: Long): Int =
    SquaresSim.squares32(u64(count), u64(key)).toInt

  /** Approx-gaussian via Irwin-Hall: sum of 12 uniform 32-bit words `rnds`
   *  (low 32 bits of each used) converted to a float. */
  def irwinHall(rnds: Array[Long]): Float =
    IrwinHallSim.irwinHall(rnds.toIndexedSeq.map(u64), 32, spec).toFloat

  /** Threefry-4x32, 20 round */
  def threefry4x32(count: Array[Int], key: Array[Int]): Array[Int] = {
    require(count.length == 4, s"threefry4x32: count must have 4 words, got ${count.length}")
    require(key.length   == 4, s"threefry4x32: key must have 4 words, got ${key.length}")
    val sim = new ThreefrySim(20, key.toIndexedSeq.map(u32))
    sim.gen(count.toIndexedSeq.map(u32)).map(_.toInt).toArray
  }

  // reinterpret a jint as its unsigned 32-bit value
  private def u32(w: Int): Long = w.toLong & 0xFFFFFFFFL
}
