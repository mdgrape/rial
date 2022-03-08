//% @file MultPackedFPGeneric.scala

package rial.simd

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._
import rial.arith._
import rial.util._

class MultPackedFPGeneric(
  len : Int, // length of SIMD pack
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  def getParam = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x   = Input (Vec(len, UInt(xSpec.W.W)))
    val y   = Input (Vec(len, UInt(ySpec.W.W)))
    val z   = Output(Vec(len, UInt(zSpec.W.W)))
  })

  val mults = (0 until len).map {i => Module(new MultFPGeneric(
    xSpec, ySpec, zSpec, roundSpec, stage
  ))}

  for(i <- 0 until len) {
    mults(i).io.x := io.x(i)
    mults(i).io.y := io.y(i)
    io.z(i) := mults(i).io.z
  }
  //printf("x=%x y=%x z=%x\n", io.x, io.y, io.z)
}
