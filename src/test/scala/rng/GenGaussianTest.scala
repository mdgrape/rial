import org.scalatest._

import chisel3._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.arith._
import rial.math._
import rial.rng._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

class SinCos2Pi(
  rndW: Int,      // width of input fixedpoint
  spec: RealSpec, // output width
  polySpec: PolynomialSpec,
  stage: PipelineStageConfig
  ) extends Module {

  val io = IO(new Bundle {
    val isSin = Input(Bool())
    val x     = Input(UInt(rndW.W))
    val z     = Output(UInt(spec.W.W))
  })

  val cbit = BoxMullerSinCos2PiTableCoeff.getCBits(polySpec)

  val preProc  = Module(new BoxMullerSinCos2PiPreProc(rndW, spec, polySpec, stage))
  val polyEval = Module(new PolynomialEval(spec, polySpec, cbit, stage))
  val postProc = Module(new BoxMullerSinCos2PiPostProc(rndW, spec, polySpec, stage))

  preProc.io.isSin := io.isSin
  preProc.io.rnd   := io.x

  polyEval.io.coeffs := preProc.io.cs
  if(polySpec.order != 0) {
    polyEval.io.dx.get := preProc.io.dx.get
  }

  postProc.io.pre  := preProc.io.out
  postProc.io.zres := polyEval.io.result

  io.z := postProc.io.z
}

class SinCos2PiTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Fixed -> FP sin(2Pi*x)/cos(2Pi*x) for BoxMuller"

  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRandomUInt(width: Int, r: Random) = {
    val rnd = SafeLong(r.nextLong)
    rnd & maskSL(width)
  }

  private def runtest ( isSin: Boolean, rndW: Int, spec : RealSpec, polySpec: PolynomialSpec,
      stage: PipelineStageConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (Int, Random) => SafeLong),
      tolerance : Int
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"sin/cos(2pi*x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new SinCos2Pi(rndW, spec, polySpec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val nstage = stage.total

          val reference = (x: SafeLong) => {
            val xr = x.toDouble / (SafeLong(1) << rndW).toDouble
            val z  = if(isSin) {
              if(x == (SafeLong(0) << (rndW-2)) || x == (SafeLong(2) << (rndW-2))) {0.0} else {sin(2 * Pi * xr)}
            } else {
              if(x == (SafeLong(1) << (rndW-2)) || x == (SafeLong(3) << (rndW-2))) {0.0} else {cos(2 * Pi * xr)}
            }
            new RealGeneric(spec, z)
          }

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(rndW, r)
            val z0r= reference(xi)
            q += ((xi.toBigInt, z0r.value.toBigInt))

            c.io.isSin.poke(isSin.B)
            c.io.x.poke(xi.toBigInt.U(spec.W.W))
            val zi = c.io.z.peek().litValue.toBigInt
            c.clock.step(1)
            if (i > nstage) {
              val (xid,z0d) = q.dequeue()

              val xidsgn = bit(spec.W-1, xid).toInt
              val xidexp = slice(spec.manW, spec.exW, xid)
              val xidman = xid & maskSL(spec.manW)

              val zisgn = bit(spec.W-1, zi).toInt
              val ziexp = slice(spec.manW, spec.exW, zi)
              val ziman = zi & maskSL(spec.manW)

              val z0dsgn = bit(spec.W-1, z0d).toInt
              val z0dexp = slice(spec.manW, spec.exW, z0d)
              val z0dman = z0d & maskSL(spec.manW)

              val diff = (zi - z0d).abs

              val xr = xid.toDouble / (SafeLong(1) << rndW).toDouble
              val xc = if(xr < 0.25) {xr} else if(xr < 0.5) {0.5 - xr} else if(xr < 0.75) {xr - 0.5} else {1.0 - xr}

              assert(diff <= tolerance,
                     f"x = ${xid.toLong.toBinaryString}(${xr} -> ${xc}), z = ${sin(2 * Pi * xr)}, "+
                     f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble}) != " +
                     f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
            }
          }

          // special values: 0, 0.25, 0.5, 0.75, 1-eps
          for(xid <- Seq(SafeLong(0), SafeLong(1)<<(rndW-2), SafeLong(2)<<(rndW-2), SafeLong(3)<<(rndW-2), maskSL(rndW) )) {
            val z0r = reference(xid)
            val z0d = z0r.value.toBigInt

            c.io.isSin.poke(isSin.B)
            c.io.x.poke(xid.toBigInt.U(spec.W.W))
            for(j <- 0 until nstage) {
              c.clock.step(1)
            }
            val zi = c.io.z.peek().litValue.toBigInt

            val diff = (zi - z0d).abs

            val xr = xid.toDouble / (SafeLong(1) << rndW).toDouble
            val xc = if(xr < 0.25) {xr} else if(xr < 0.5) {0.5 - xr} else if(xr < 0.75) {xr - 0.5} else {1.0 - xr}
            val zr = if(isSin) {sin(2 * Pi * xr)} else {cos(2 * Pi * xr)}

            val xidsgn = bit(spec.W-1, xid).toInt
            val xidexp = slice(spec.manW, spec.exW, xid)
            val xidman = xid & maskSL(spec.manW)

            val zisgn = bit(spec.W-1, zi).toInt
            val ziexp = slice(spec.manW, spec.exW, zi)
            val ziman = zi & maskSL(spec.manW)

            val z0dsgn = bit(spec.W-1, z0d).toInt
            val z0dexp = slice(spec.manW, spec.exW, z0d)
            val z0dman = z0d & maskSL(spec.manW)

            if(diff > tolerance) {
              c.clock.step(1) // run printf
            }

            assert(diff <= tolerance,
                   f"x = ${xid.toLong.toBinaryString}(${xr} -> ${xc}), z(scala.math.FP64) = ${zr}, "+
                   f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble}) != " +
                   f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
          }
        }
      }
    }
  }
  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3
  val polySpecFP32  = new PolynomialSpec(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)

  runtest(true,  32, RealSpec.Float32Spec, polySpecFP32, PipelineStageConfig.none,
    n, r, "Test FP32 sin",generateRandomUInt, 3)
  runtest(false, 32, RealSpec.Float32Spec, polySpecFP32, PipelineStageConfig.none,
    n, r, "Test FP32 cos",generateRandomUInt, 3)
}


class Sqrt2LogXTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Fixed -> FP sqrt(-2log(x)) for Box-Muller"

  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRandomUInt(width: Int, r: Random) = {
    val rnd = SafeLong(r.nextLong)
    rnd & maskSL(width)
  }

  private def runtest ( rndW: Int, spec : RealSpec, polySpec: PolynomialSpec,
      stage: PipelineStageConfig,
      n : Int, r : Random, generatorStr : String,
      generator : ( (Int, Random) => SafeLong),
      tolerance : Double
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"sqrt(-2log(x)) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new Sqrt2LogX(rndW, spec, polySpec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val nstage = stage.total

          val reference = (x: SafeLong) => {
            val xr = x.toDouble / (SafeLong(1) << rndW).toDouble
            val z  = sqrt(-2.0 * log(1.0 - xr))
            new RealGeneric(spec, z)
          }

          val q  = new Queue[(BigInt,RealGeneric)]
          for(i <- 1 to n+nstage) {
            val xi  = generator(rndW, r)
            val z0r = reference(xi)
            q += ((xi.toBigInt, z0r))

            c.io.x.poke(xi.toBigInt.U(spec.W.W))
            val zi = c.io.z.peek().litValue.toBigInt
            c.clock.step(1)
            if (i > nstage) {
              val (xid,z0d) = q.dequeue()

              val xidsgn = bit(spec.W-1, xid).toInt
              val xidexp = slice(spec.manW, spec.exW, xid)
              val xidman = xid & maskSL(spec.manW)

              val zisgn = bit(spec.W-1, zi).toInt
              val ziexp = slice(spec.manW, spec.exW, zi)
              val ziman = zi & maskSL(spec.manW)
              val zid = new RealGeneric(spec, zisgn, ziexp.toInt, ziman)

              val z0dsgn = z0d.sgn
              val z0dexp = z0d.ex
              val z0dman = z0d.man

              val diff = (zid.toDouble - z0d.toDouble).abs

              val xr = xid.toDouble / (SafeLong(1) << rndW).toDouble

              assert(diff <= Seq(tolerance, z0d.toDouble * tolerance).max,
                     f"x = ${xid.toLong.toBinaryString}(${xr}), z = ${sqrt(-2.0*log(1.0 - xr))}, "+
                     f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble}) != " +
                     f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
            }
          }
        }
      }
    }
  }
  val rndW = 32
  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3 // should be >= 3
  val polySpecFP32  = new PolynomialSpec(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32, Some(rndW - adrWFP32))

  runtest(rndW, RealSpec.Float32Spec, polySpecFP32, PipelineStageConfig.none,
    n, r, "Test Fixed sqrt(-2log(x)) to FP32",generateRandomUInt, 15 * pow(2.0, -23))
}
