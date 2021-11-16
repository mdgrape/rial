package rial.fixed

class FixedSpec (
  val totalW   : Int,
  val fracW    : Int,
  val signed   : Boolean = false, // 2's complement
  val saturate : Boolean = false
  ) {
  val intW = totalW - fracW
  def W : Int = {totalW}
}
