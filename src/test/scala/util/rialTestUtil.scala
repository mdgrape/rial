package rial.testUtil

import chisel3._
import chisel3.util._
import scala.math._
import scala.util.Random

import spire.math.SafeLong

import org.scalatest.TestSuite

import chisel3.simulator.HasSimulator
import chisel3.simulator.scalatest.ChiselSim
import svsim.verilator.Backend.CompilationSettings
import svsim.verilator.Backend.CompilationSettings.Parallelism

import rial.arith._
import rial.util.ScalaUtil._

object ScalaTestUtil {

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealWithin( from: Double, to: Double, spec: RealSpec, r: Random ) = {
    val rD : Double = (r.nextDouble())*(to - from) + from
    new RealGeneric(spec, rD)
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def generateRealWithinPos( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = r.nextDouble()*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFullPos ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W-1, r)))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }
  def errorLSB( x : RealGeneric, y : Float ) : Float = {
    val err = x.toFloat - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }
  def errorLSB( x : RealGeneric, y : RealGeneric ) : SafeLong = {
     x.value - y.value
  }


}

//
// ChiselSim with single-threaded verilation.
//
// Verilator 5.020 intermittently aborts with
//   %Error: Internal Error: attempted to destroy locked Thread Pool
// while verilating. The ChiselSim default backend settings pass `-j 0`, which
// makes both the verilation stage and the C++ build use every core; the bug is
// in the verilation thread pool only. Here verilation is pinned to a single
// thread (`--verilate-jobs 1`) while the C++ build keeps every core
// (`--build-jobs 0`), so suites stay fast. Everything else is taken from the
// default settings.
//
trait RialChiselSim extends ChiselSim {
  self : TestSuite =>

  implicit def rialSimulator : HasSimulator = {
    HasSimulator.simulators.verilator(
      verilatorSettings = CompilationSettings.default.withParallelism(
        Some(Parallelism.Different.default.withBuild(Some(0)).withVerilate(Some(1)))
      )
    )
  }
}
