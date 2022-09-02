//% @file CrossMultPackedFPGeneric.scala

package rial.simd

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._
import rial.arith._
import rial.util._

/** Performs cross product of 3d vectors.
 *
 * {{{
 * z[0] = x[1] * y[2]
 * z[1] = x[2] * y[0]
 * z[2] = x[0] * y[1]
 * z[3] = x[3] * y[3]
 * }}}
 *
 * @constructor create a new CrossMultPackedFPGeneric module.
 * @param xSpec the spec of input floating point number
 * @param ySpec the spec of input floating point number
 * @param zSpec the spec of output floating point number
 * @param roundSpec the spec of rounding
 * @param stage the number of pipeline stages.
 *
 */
class CrossMultPackedFPGeneric(
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  def getParam = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x = Input (Vec(4, UInt(xSpec.W.W)))
    val y = Input (Vec(4, UInt(ySpec.W.W)))
    val z = Output(Vec(4, UInt(zSpec.W.W)))
  })

  val mults = (0 until 4).map {i => Module(new MultFPGeneric(
    xSpec, ySpec, zSpec, roundSpec, stage
  ))}

  mults(0).io.x := io.x(1)
  mults(1).io.x := io.x(2)
  mults(2).io.x := io.x(0)
  mults(3).io.x := io.x(3)

  mults(0).io.y := io.y(2)
  mults(1).io.y := io.y(0)
  mults(2).io.y := io.y(1)
  mults(3).io.y := io.y(3)

  io.z(0) := mults(0).io.z
  io.z(1) := mults(1).io.z
  io.z(2) := mults(2).io.z
  io.z(3) := mults(3).io.z
  //printf("x=%x y=%x z=%x\n", io.x, io.y, io.z)
}
