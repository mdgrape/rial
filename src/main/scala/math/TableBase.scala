package rial.math

import chisel3._
import chisel3.util._

import rial.table.{PolynomialBitWidth, PolynomialCoefficients}

private[rial] trait TableBase {
  def getBitWidth(): PolynomialBitWidth
  def generateTable(bitW: PolynomialBitWidth): PolynomialCoefficients
}

/** The coefficient tables of every function of a configuration, in one ROM.
 *
 * `io.adr` is a row of the merged ROM, which [[RomLayout]] lays out: the row
 * the table of a function starts at, with the address of that table below it.
 * The coefficients of that row leave on `io.cs`, `stage` cycles later, packed
 * as [[rial.table.PolynomialGenerator.packTable]] packs them; one row leaves
 * per cycle, so only one function may read the ROM at a time.
 *
 *  @param layout where the table of each function sits in the ROM
 *  @param stage  pipeline registers between the address and the row
 */
class TableCoeffRom(
  val layout: RomLayout,
  val stage:  Int,
) extends Module {

  val io = IO(new Bundle {
    val adr = Flipped(Valid(UInt(layout.addrW.W)))
    val cs  = Valid(UInt(layout.rowW.W))
  })

  val packed = VecInit(layout.packed.map(c => c.U(layout.rowW.W)))
  io.cs := Pipe(io.adr.valid, packed(io.adr.bits), stage)
}
