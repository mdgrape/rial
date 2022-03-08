package rial.simd

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._
import rial.arith._
import rial.util._

/**
  * Compute addition of packed floating point numbers.
  *
  */
class AddPackedFPGeneric(
  len : Int, // length of SIMD pack
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig,
  val enableDebug : Boolean = false
) extends Module with DebugControlSlave {

  val nStage = stage.total

  def getParam = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle {
    val x = Input (Vec(len, UInt(xSpec.W.W)))
    val y = Input (Vec(len, UInt(ySpec.W.W)))
    val z = Output(Vec(len, UInt(zSpec.W.W)))
  })

  val adders = (0 until len).map {i => Module(new AddFPGeneric(
    xSpec, ySpec, zSpec, roundSpec, stage
  ))}

  for(i <- 0 until len) {
    adders(i).io.x := io.x(i)
    adders(i).io.y := io.y(i)
    io.z(i) := adders(i).io.z
  }
}
