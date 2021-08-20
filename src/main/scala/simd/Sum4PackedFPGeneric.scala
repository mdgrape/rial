package rial.simd

import scala.collection.mutable.Queue
import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._
import rial.arith._
import rial.util._

/**
  * Compute horizontal sum of packed floating point numbers.
  *
  */
class Sum4PackedFPGeneric(
  xSpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig,
  val enableDebug : Boolean = false
) extends MultiIOModule with DebugControlSlave {

  val depth  = 2
  val nStage = stage.total * depth

  def getParam() = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage() = nStage

  val io = IO(new Bundle {
    val x    = Input (Vec(4, UInt(xSpec.W.W)))
    val mask = Input (Vec(4, Bool()))
    val z    = Output(UInt(zSpec.W.W))
  })

  val adder1 = Module(new AddFPGeneric(xSpec, xSpec, xSpec, roundSpec, stage));
  val adder2 = Module(new AddFPGeneric(xSpec, xSpec, xSpec, roundSpec, stage));
  val adder3 = Module(new AddFPGeneric(xSpec, xSpec, zSpec, roundSpec, stage));

  adder1.io.x := Mux(io.mask(0), 0.U, io.x(0));
  adder1.io.y := Mux(io.mask(1), 0.U, io.x(1));
  adder2.io.x := Mux(io.mask(2), 0.U, io.x(2));
  adder2.io.y := Mux(io.mask(3), 0.U, io.x(3));

  adder3.io.x := adder1.io.z
  adder3.io.y := adder2.io.z
  io.z := adder3.io.z
}
