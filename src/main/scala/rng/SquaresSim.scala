//% @file SquaresSim.scala
//
// Bit-exact software model of the Squares32 hardware generator (squares.scala).
// Given the same (count, key) the HW consumes, it returns the identical 32-bit
// random word. The 4 pipeline stages do not change the functional result, so
// the model is a pure function of (count, key).
//
//   Squares: A Fast Counter-Based RNG, B. Widynski, 2022 (arXiv:2004.06278)
//
package rial.rng

object SquaresSim {

  private val mask64 = (BigInt(1) << 64) - 1
  private val mask32 = (BigInt(1) << 32) - 1

  /** Bit-exact model of [[Squares32]]: the 32-bit random word for (count, key).
   *  Masking before/after each multiply is equivalent mod 2^64, so this matches
   *  the HW (which keeps every intermediate in 64 bits). */
  def squares32(count: BigInt, key: BigInt): BigInt = {
    var x = (count * key) & mask64
    val y = x
    val z = (y + key) & mask64

    // one round: x = rotate32(x*x + addend)
    def round(addend: BigInt): Unit = {
      x = (((x * x) & mask64) + addend) & mask64
      x = ((x >> 32) | ((x & mask32) << 32)) & mask64
    }
    round(y) // round 1
    round(z) // round 2
    round(y) // round 3
    round(z) // round 4

    // HW returns stage3.z(31,0) = high 32 bits of (x*x + z) = rotate(...)& mask32
    x & mask32
  }
}
