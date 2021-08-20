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
class HSumPackedFPGeneric(
  len : Int, // length of SIMD pack
  xSpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig,
  val enableDebug : Boolean = false
) extends MultiIOModule with DebugControlSlave {

  val depth  = log2Up(len)
  val nStage = stage.total * depth

  def getParam() = { (xSpec, zSpec, roundSpec, nStage) }

  def getStage() = nStage

  val io = IO(new Bundle {
    val x    = Input (Vec(len, UInt(xSpec.W.W)))
    val mask = Input (Vec(len, Bool()))
    val z    = Output(UInt(zSpec.W.W))
  })

  // timing of conversion from xSpec to zSpec..
  val adders = (0 until len-2).map {i => Module(new AddFPGeneric(
    xSpec, xSpec, xSpec, roundSpec, stage
  ))}

  val in = Wire(Vec(len, UInt(xSpec.W.W)))
  for(i <- 0 until len) {
    in(i) := Mux(io.mask(i), 0.U, io.x(i))
  }

  var q = new Queue[UInt]
  for(i <- 0 until len) {
    q += in(i)
  }

  var idx = 0
  for(idx <- 0 until len-2) {
    adders(idx).io.x := q.dequeue
    adders(idx).io.y := q.dequeue
    q += adders(idx).io.z
  }
  assert(q.length == 2)

  // convert to zSpec
  val lastAdder = Module(new AddFPGeneric(xSpec, xSpec, zSpec, roundSpec, stage))

  lastAdder.io.x := q.dequeue
  lastAdder.io.y := q.dequeue
  io.z := lastAdder.io.z
}
