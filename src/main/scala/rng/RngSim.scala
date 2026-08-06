//% @file RngSim.scala
//
// Stable software-simulation facade for the JNI binding (jni/, rialcpp), for
// the random-number generators. Mirrors rial.math.MathSim: the C++ side looks
// up exactly these names, so the binding does not depend on the internal
// generator/Sim object names.
//
// It comes in two halves:
//
// - [[RngSim]] holds the generators that produce a float, so they depend on the
//   floating point format and on how many random bits they consume. Like
//   MathSim, one instance == one configuration.
// - [[BitRngSim]] holds the counter-based generators that produce raw random
//   bits. They have nothing to do with the floating point format, so they are
//   stateless and take their parameters (the Threefry round count) per call.
//
// The RngSim constructor takes only `Int`s so that the C++ side can build an
// instance through plain JNI (`NewObject` with a "(IIIIIIII)V" signature).
//
// Every float-returning entry takes the raw random bits as an unsigned 64-bit
// value (`Long` bit pattern).
//
package rial.rng

import rial.arith.RealSpec

/** Software simulation of the float-producing random-number generators.
 *
 *  @constructor set up the generator configuration.
 *  @param exW                 exponent width of the generated floats.
 *  @param exBias              exponent bias of the generated floats.
 *  @param manW                mantissa width of the generated floats.
 *  @param float01W            width of the random fraction [[genFloat01]] uses.
 *  @param gaussianW           width of the random word [[genGaussian]] uses.
 *  @param gaussianOrder       polynomial order of the probit approximation.
 *  @param gaussianLumpOctaves how many of the deepest tail octaves are merged
 *                             into one interval of the probit table.
 *  @param irwinHallW          width of each of the 12 words [[irwinHall]] sums.
 */
class RngSim(
  val exW:                 Int,
  val exBias:              Int,
  val manW:                Int,
  val float01W:            Int,
  val gaussianW:           Int,
  val gaussianOrder:       Int,
  val gaussianLumpOctaves: Int,
  val irwinHallW:          Int,
) {
  val spec = new RealSpec(exW, exBias, manW)

  /** uniform float in [1, 2) -- low manW bits of `rnd` become the mantissa. */
  def genFloat12(rnd: Long): Float =
    GenRandomFloatSim.genFloat12(RngSim.u64(rnd), spec).toFloat

  /** uniform float in [0, 1) from a `float01W`-bit random fraction. */
  def genFloat01(rnd: Long): Float =
    GenRandomFloatSim.genFloat01(RngSim.u64(rnd), float01W, spec).toFloat

  /** N(0,1) gaussian from a `gaussianW`-bit random word. */
  def genGaussian(rnd: Long): Float =
    GenRandomGaussianSim.genGaussian(RngSim.u64(rnd), gaussianW, spec,
      gaussianOrder, gaussianLumpOctaves).toFloat

  /** Approx-gaussian via Irwin-Hall: sum of 12 uniform `irwinHallW`-bit words
   *  `rnds` (low bits of each used) converted to a float. */
  def irwinHall(rnds: Array[Long]): Float =
    IrwinHallSim.irwinHall(rnds.toIndexedSeq.map(RngSim.u64), irwinHallW, spec).toFloat
}

/** Default configurations for [[RngSim]]. */
object RngSim {
  // reinterpret a (possibly negative) jlong as its unsigned 64-bit value
  private[rng] def u64(rnd: Long): BigInt = BigInt(rnd) & ((BigInt(1) << 64) - 1)

  /** The FP32 configuration the binding used before it became configurable. */
  def standardFP32: RngSim = new RngSim(
    exW                 = RealSpec.Float32Spec.exW,
    exBias              = RealSpec.Float32Spec.exBias,
    manW                = RealSpec.Float32Spec.manW,
    float01W            = 32,
    gaussianW           = 64,
    gaussianOrder       = 3,
    gaussianLumpOctaves = 2,
    irwinHallW          = 32,
  )
}

/** Software simulation of the counter-based raw bit generators.
 */
object BitRngSim {

  /** Squares32 counter-based uniform RNG: the 32-bit random word for
   *  (count, key), returned as the raw bit pattern in a jint. */
  def squares32(count: Long, key: Long): Int = {
    SquaresSim.squares32(RngSim.u64(count), RngSim.u64(key)).toInt
  }

  /** Threefry-4x32: the four 32-bit random words for (count, key).
   *
   *  @param rounds number of rounds (20 for the standard Threefry-4x32-20).
   */
  def threefry4x32(count: Array[Int], key: Array[Int], rounds: Int): Array[Int] = {
    require(count.length == 4, s"threefry4x32: count must have 4 words, got ${count.length}")
    require(key.length   == 4, s"threefry4x32: key must have 4 words, got ${key.length}")
    val sim = new ThreefrySim(rounds, key.toIndexedSeq.map(u32))
    sim.gen(count.toIndexedSeq.map(u32)).map(_.toInt).toArray
  }

  private def u32(w: Int): Long = { w.toLong & 0xFFFFFFFFL }
}
