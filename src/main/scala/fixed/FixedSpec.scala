package rial.fixed

/** Fixed point spec.
 *
 * It contains the width of the fractional part.
 * In conversion modules, the fracW value can be overwritten by the second
 * (optional) input. If the second port is not defined (set None), the default
 * fracW is used. Otherwise, fracW will not be referenced.
 *
 * @param totalW   total width of the fixed point number.
 * @param fracW    the default length of the fractional part.
 * @param signed   spec of the output floating point number.
 * @param saturate Determines how many registers needed.
 */
case class FixedSpec (
  val totalW:    Int,
  val fracW:     Int,
  val signed:    Boolean = true, // 2's complement
  val saturate:  Boolean = true
) {
  val intW = totalW - fracW
  def W : Int = {
    totalW
  }
}
