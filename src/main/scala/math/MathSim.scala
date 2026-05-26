//% @file MathSim.scala
//
// Stable software-simulation facade for the JNI binding (jni/, rialcpp).
//
// The C++ side looks up exactly this object and these method names, so the
// binding does not depend on the internal per-function Sim/Table object names
// (which change as math/ is refactored). All entries use the `standard` FP32
// configuration. To expose more (BFloat16/FP8, custom nOrder/bitW, ...) add
// methods here rather than re-pointing the binding at internal objects.
//
package rial.math

object MathSim {
  def sqrt      (x: Float)          : Float = SqrtTable.sqrtDefaultFP32(x)
  def invsqrt   (x: Float)          : Float = InvSqrtTable.invsqrtDefaultFP32(x)
  def reciprocal(x: Float)          : Float = ReciprocalTable.reciprocalDefaultFP32(x)
  def sin       (x: Float)          : Float = SinCosTable.sinDefaultFP32(x)
  def cos       (x: Float)          : Float = SinCosTable.cosDefaultFP32(x)
  def acos      (x: Float)          : Float = ACosPhase2Sim.acosDefaultFP32(x)
  def atan2     (y: Float, x: Float): Float = ATan2Phase2Sim.atan2DefaultFP32(y, x)
  def exp       (x: Float)          : Float = ExpTable.expDefaultFP32(x)
  def log       (x: Float)          : Float = LogSim.logDefaultFP32(x)
}
