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

class BoxMullerSinCos2Pi(
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

  val preProc  = Module(new BoxMullerSinCos2PiPreProc(rndW, spec, polySpec, cbit, stage))
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

class BoxMullerSinCos2PiTest extends AnyFlatSpec
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
      test( new BoxMullerSinCos2Pi(rndW, spec, polySpec, stage)).
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

// ============================================================================

class BoxMullerLog(
  rndW: Int,      // width of input fixedpoint
  spec: RealSpec, // output width
  polySpec: PolynomialSpec,
  stage: PipelineStageConfig
  ) extends Module {

  val io = IO(new Bundle {
    val x     = Input(UInt(rndW.W))
    val z     = Output(UInt(spec.W.W))
  })

  val cbit = BoxMullerLogTableCoeff.getCBits(polySpec)

  val preProc  = Module(new BoxMullerLogPreProc(rndW, spec, polySpec, cbit, stage))
  val polyEval = Module(new PolynomialEval(spec, polySpec, cbit, stage))
  val postProc = Module(new BoxMullerLogPostProc(rndW, spec, polySpec, stage))

  preProc.io.x := io.x

  polyEval.io.coeffs := preProc.io.cs
  if(polySpec.order != 0) {
    polyEval.io.dx.get := preProc.io.dx.get
  }

  postProc.io.pre  := preProc.io.out
  postProc.io.zres := polyEval.io.result

  io.z := postProc.io.z
}

class BoxMullerLogTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Fixed -> FP -2log(1-x) for BoxMuller"

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
      tolerance : Int
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"-2log(1-x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new BoxMullerLog(rndW, spec, polySpec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val nstage = stage.total

          val reference = (x: SafeLong) => {
            val xr = ((SafeLong(1) << rndW) - x).toDouble / (SafeLong(1) << rndW).toDouble
            val z  = -2 * log(xr)
            new RealGeneric(spec, z)
          }

          var maxError    = 0.0
          var xatMaxError = 0.0
          var zatMaxError = 0.0
          val errs = collection.mutable.Map[Int, (Int, Int)]()

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(rndW, r)
            val z0r= reference(xi)
            q += ((xi.toBigInt, z0r.value.toBigInt))

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

              // relative tolerance
              val exBias = spec.exBias
              val exdiff = (ziexp - z0dexp.toInt).abs.toInt
              val diff   = if(z0dexp < exBias) {
                (zi - z0d).abs >> (exBias-z0dexp.toInt)
              } else {
                (zi - z0d).abs
              }

              val erri = (zi - z0d).toLong
              if(zi - z0d != 0) {
                val errkey = erri.abs.toInt
                if( ! errs.contains(errkey)) {
                  errs(errkey) = (0, 0)
                }
                if (erri >= 0) {
                  errs(errkey) = (errs(errkey)._1 + 1, errs(errkey)._2)
                } else {
                  errs(errkey) = (errs(errkey)._1, errs(errkey)._2 + 1)
                }
              }
              if (maxError < erri.abs) {
                maxError    = erri.abs
                xatMaxError = xid.toDouble / (SafeLong(1) << rndW).toDouble
                zatMaxError = new RealGeneric(spec, zisgn.toInt, ziexp.toInt, ziman).toDouble
              }

              val xr = ((SafeLong(1) << rndW) - xid).toDouble / (SafeLong(1) << rndW).toDouble

              assert((exdiff == 0) && (diff <= tolerance),
                     f"x = ${xid.toLong.toBinaryString}, z = ${-2*log(xr)}, "+
                     f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble}) != " +
                     f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
            }
          }
          println(f"${generatorStr} Summary")
          if(maxError != 0.0) {
            println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
                  + f"${zatMaxError} != ${-2*log(1.0 - xatMaxError)}, "
                  + f"diff = ${zatMaxError - (-2*log(1.0 - xatMaxError))}, x = ${xatMaxError}")
          }
          // if the resulting z is small, error becomes too large to print them all.
//           for(kv <- errs.toSeq.sortBy(_._1)) {
//             val (k, (errPos, errNeg)) = kv
//             println(f"N=$n%d : +/- $k errors positive $errPos%d / negative $errNeg%d")
//           }
          println( "---------------------------------------------------------------")

          // special values: 0, 1-eps
          for(xid <- Seq(SafeLong(0), maskSL(rndW))) {
            val z0r = reference(xid)
            val z0d = z0r.value.toBigInt

            c.io.x.poke(xid.toBigInt.U(spec.W.W))
            for(j <- 0 until nstage) {
              c.clock.step(1)
            }
            val zi = c.io.z.peek().litValue.toBigInt

            val diff = (zi - z0d).abs

            val xr = ((SafeLong(1) << rndW) - xid).toDouble / (SafeLong(1) << rndW).toDouble
            val zr = -2*log(xr)

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

            val x0r = xid.toDouble / (SafeLong(1) << rndW).toDouble
            assert(diff <= tolerance,
                   f"x = ${xid.toLong.toBinaryString}(${x0r} -> ${xr}), z(scala.math.FP64) = ${zr}, "+
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

  runtest(32, RealSpec.Float32Spec, polySpecFP32, PipelineStageConfig.none,
    n, r, "Test FP32 log",generateRandomUInt, 3)
}

// ============================================================================

class BoxMullerSqrt(
  rndW: Int,      // width of input fixedpoint
  spec: RealSpec, // output width
  polySpec: PolynomialSpec,
  stage: PipelineStageConfig
  ) extends Module {

  val io = IO(new Bundle {
    val x     = Input(UInt(spec.W.W))
    val z     = Output(UInt(spec.W.W))
  })

  val cbit = SqrtTableCoeff.getCBits(spec, polySpec)

  val preProc  = Module(new BoxMullerSqrtPreProc(rndW, spec, polySpec, cbit, stage))
  val polyEval = Module(new PolynomialEval(spec, polySpec, cbit, stage))
  val postProc = Module(new BoxMullerSqrtPostProc(rndW, spec, polySpec, stage))

  preProc.io.x := io.x

  polyEval.io.coeffs := preProc.io.cs
  if(polySpec.order != 0) {
    polyEval.io.dx.get := preProc.io.dx.get
  }

  postProc.io.pre  := preProc.io.out
  postProc.io.zres := polyEval.io.result

  io.z := postProc.io.z
}

class BoxMullerSqrtTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Fixed -> FP sqrt(x) for BoxMuller"

  var n = 10000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("10000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)


  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble())*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  private def runtest ( rndW: Int, spec : RealSpec, polySpec: PolynomialSpec,
    stage: PipelineStageConfig,
    n : Int, r : Random, generatorStr : String,
    generator : ( (RealSpec, Random) => RealGeneric),
    tolerance : Int
  ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"sqrt(x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new BoxMullerSqrt(rndW, spec, polySpec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val nstage = stage.total

          val reference = (x: RealGeneric) => {
            val z = sqrt(x.toDouble)
            new RealGeneric(spec, z)
          }

          var maxError    = 0.0
          var xatMaxError = 0.0
          var zatMaxError = 0.0
          val errs = collection.mutable.Map[Int, (Int, Int)]()

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec, r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt, z0r.value.toBigInt))

            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
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

              val xr = new RealGeneric(spec, xidsgn, xidexp.toInt, xidman)
              val zr = new RealGeneric(spec, zisgn,  ziexp.toInt,  ziman)

              val erri = (zi - z0d).toLong
              if(erri.abs != 0) {
                val errkey = erri.abs.toInt
                if( ! errs.contains(errkey)) {
                  errs(errkey) = (0, 0)
                }
                if (erri >= 0) {
                  errs(errkey) = (errs(errkey)._1 + 1, errs(errkey)._2)
                } else {
                  errs(errkey) = (errs(errkey)._1, errs(errkey)._2 + 1)
                }
              }
              if (maxError < erri.abs) {
                maxError    = erri.abs
                xatMaxError = xr.toDouble
                zatMaxError = zr.toDouble
              }

              assert(erri.abs <= tolerance,
                     f"x = ${xr.toDouble}, z = ${zr.toDouble}" +
                     f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble}) != " +
                     f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
            }
          }
          println(f"${generatorStr} Summary")
          if(maxError != 0.0) {
            println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
                  + f"${zatMaxError} != ${sqrt(xatMaxError)}, "
                  + f"diff = ${zatMaxError - sqrt(xatMaxError)}, x = ${xatMaxError}")
          }
          for(kv <- errs.toSeq.sortBy(_._1)) {
            val (k, (errPos, errNeg)) = kv
            println(f"N=$n%d : +/- $k errors positive $errPos%d / negative $errNeg%d")
          }
          println( "---------------------------------------------------------------")
        }
      }
    }
  }
  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3
  val polySpecFP32  = new PolynomialSpec(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)

  runtest(32, RealSpec.Float32Spec, polySpecFP32, PipelineStageConfig.none,
    n, r, "Test FP32 log",generateRealWithin(128.0,_,_), 3)
}
