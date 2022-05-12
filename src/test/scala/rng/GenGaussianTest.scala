import org.scalatest._

import chisel3._
import chiseltest._
import chiseltest.VerilatorBackendAnnotation

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import org.apache.commons.math3.special.Erf

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
    val en    = Input(Bool())
    val isSin = Input(Bool())
    val x     = Input(UInt(rndW.W))
    val z     = Output(UInt(spec.W.W))
  })

  val cbit = Seq(
    BoxMullerSinCos2PiTableCoeff.getCBits(polySpec),
    BoxMullerLogTableCoeff.getCBits(polySpec),
    SqrtTableCoeff.getCBits(spec, polySpec)
    ).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  val preProc  = Module(new BoxMullerSinCos2PiPreProc(rndW, spec, polySpec, cbit, stage))
  val polyEval = Module(new PolynomialEval(spec, polySpec, cbit, stage))
  val postProc = Module(new BoxMullerSinCos2PiPostProc(rndW, spec, polySpec, stage))

  preProc.io.en    := io.en
  preProc.io.isSin := io.isSin
  preProc.io.rnd   := io.x

  polyEval.io.coeffs := preProc.io.cs
  if(polySpec.order != 0) {
    polyEval.io.dx.get := preProc.io.dx.get
  }

  postProc.io.en    := io.en
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

          var maxError    = 0.0
          var xatMaxError = 0.0
          var zatMaxError = 0.0
          val errs = collection.mutable.Map[Int, (Int, Int)]()

          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(rndW, r)
            val z0r= reference(xi)
            q += ((xi.toBigInt, z0r.value.toBigInt))

            c.io.en.poke(true.B)
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
                zatMaxError = new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble
              }

              assert(diff <= tolerance,
                     f"x = ${xid.toLong.toBinaryString}(${xr} -> ${xc}), z = ${if(isSin){sin(2 * Pi * xr)} else {cos(2 * Pi * xr)}}, "+
                     f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble}) != " +
                     f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
            }
          }
          println(f"${generatorStr} Summary")
          if(maxError != 0.0) {
            println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
                  + f"${zatMaxError} != ${if(isSin){sin(2*Pi*xatMaxError)}else{cos(2*Pi*xatMaxError)}}, "
                  + f"diff = ${zatMaxError - (if(isSin){sin(2*Pi*xatMaxError)}else{cos(2*Pi*xatMaxError)})}, x = ${xatMaxError}")
          }
          for(kv <- errs.toSeq.sortBy(_._1)) {
            val (k, (errPos, errNeg)) = kv
            println(f"N=$n%d : +/- $k errors positive $errPos%d / negative $errNeg%d")
          }
          println( "---------------------------------------------------------------")

          maxError    = 0.0
          xatMaxError = 0.0
          zatMaxError = 0.0
          errs.clear()

          // special values: 0, 0.25, 0.5, 0.75, 1-eps
          for(x0 <- Seq(SafeLong(0), SafeLong(1)<<(rndW-2), SafeLong(2)<<(rndW-2), SafeLong(3)<<(rndW-2), maskSL(rndW) )) {
            for(i <- -255 to 255) {

              val xid = if(x0 == maskSL(rndW)) {x0 - i.abs} else {(x0 + i).abs}

              val z0r = reference(xid)
              val z0d = z0r.value.toBigInt

              c.io.en.poke(true.B)
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
                zatMaxError = new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble
              }

              assert(diff <= tolerance,
                     f"x = ${xid.toLong.toBinaryString}(${xr} -> ${xc}), z(scala.math.FP64) = ${zr}, "+
                     f"test(${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString})(${new RealGeneric(spec, zisgn, ziexp.toInt, ziman).toDouble}) != " +
                     f"ref(${z0dsgn}|${z0dexp}(${z0dexp - spec.exBias})|${z0dman.toLong.toBinaryString})(${new RealGeneric(spec, z0dsgn, z0dexp.toInt, z0dman).toDouble})")
            }
          }
          println(f"${generatorStr} Summary")
          if(maxError != 0.0) {
            println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
                  + f"${zatMaxError} != ${if(isSin){sin(2*Pi*xatMaxError)}else{cos(2*Pi*xatMaxError)}}, "
                  + f"diff = ${zatMaxError - (if(isSin){sin(2*Pi*xatMaxError)}else{cos(2*Pi*xatMaxError)})}, x = ${xatMaxError}")
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

  runtest(true,  32, RealSpec.Float32Spec, polySpecFP32, PipelineStageConfig.none,
    n, r, "Test FP32 sin",generateRandomUInt, 3)
  runtest(false, 32, RealSpec.Float32Spec, polySpecFP32, PipelineStageConfig.none,
    n, r, "Test FP32 cos",generateRandomUInt, 7)
}

// ============================================================================

class BoxMullerLog(
  rndW: Int,      // width of input fixedpoint
  spec: RealSpec, // output width
  polySpec: PolynomialSpec,
  stage: PipelineStageConfig
  ) extends Module {

  val io = IO(new Bundle {
    val en    = Input(Bool())
    val x     = Input(UInt(rndW.W))
    val z     = Output(UInt(spec.W.W))
  })

  val cbit = Seq(
    BoxMullerSinCos2PiTableCoeff.getCBits(polySpec),
    BoxMullerLogTableCoeff.getCBits(polySpec),
    SqrtTableCoeff.getCBits(spec, polySpec)
    ).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  val preProc  = Module(new BoxMullerLogPreProc(rndW, spec, polySpec, cbit, stage))
  val polyEval = Module(new PolynomialEval(spec, polySpec, cbit, stage))
  val postProc = Module(new BoxMullerLogPostProc(rndW, spec, polySpec, stage))

  preProc.io.x := io.x
  preProc.io.en := io.en

  polyEval.io.coeffs := preProc.io.cs
  if(polySpec.order != 0) {
    polyEval.io.dx.get := preProc.io.dx.get
  }

  postProc.io.en   := io.en
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

            c.io.en.poke(true.B)
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

            c.io.en.poke(true.B)
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
    val en    = Input(Bool())
    val x     = Input(UInt(spec.W.W))
    val z     = Output(UInt(spec.W.W))
  })

  val cbit = Seq(
    BoxMullerSinCos2PiTableCoeff.getCBits(polySpec),
    BoxMullerLogTableCoeff.getCBits(polySpec),
    SqrtTableCoeff.getCBits(spec, polySpec)
    ).reduce( (lhs, rhs) => { lhs.zip(rhs).map( x => max(x._1, x._2) ) } )

  val preProc  = Module(new BoxMullerSqrtPreProc(rndW, spec, polySpec, cbit, stage))
  val polyEval = Module(new PolynomialEval(spec, polySpec, cbit, stage))
  val postProc = Module(new BoxMullerSqrtPostProc(rndW, spec, polySpec, stage))

  preProc.io.en := io.en
  preProc.io.x := io.x

  polyEval.io.coeffs := preProc.io.cs
  if(polySpec.order != 0) {
    polyEval.io.dx.get := preProc.io.dx.get
  }

  postProc.io.en   := io.en
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

            c.io.en.poke(true.B)
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


//
// Test BoxMuller
//
class BoxMullerTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test normal distribution, mean = 0, stddev = 1"

  var n = 100000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("100000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRandomUInt(width: Int) = {
    val rnd = SafeLong(r.nextLong)
    rnd & maskSL(width)
  }

  private def runTest( spec: RealSpec, polySpec: PolynomialSpec, roundSpec: RoundSpec, tolerance: Int ) = {
    it should f"BoxMuller(x) spec ${spec.toStringShort} comparison" in {
      test( new BoxMuller(32, spec, polySpec, roundSpec)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val rndW = 32

          val reference = (xi: SafeLong, yi: SafeLong) => {
            val x  = ((SafeLong(1) << rndW) - xi).toDouble / (SafeLong(1) << rndW).toDouble
            val y  = yi.toDouble / (SafeLong(1) << rndW).toDouble
            val z1 = sqrt(-2.0 * log(x)) * sin(2.0 * Pi * y)
            val z2 = sqrt(-2.0 * log(x)) * cos(2.0 * Pi * y)
            val zr1 = new RealGeneric(spec, z1)
            val zr2 = new RealGeneric(spec, z2)
            (zr1, zr2)
          }

          var maxError    = 0.0
          var xatMaxError = 0.0
          var yatMaxError = 0.0
          var zatMaxError = 0.0
          var satMaxError = false
          val errs = collection.mutable.Map[Int, (Int, Int)]()

          //           pre   poly  post  multiply        | consume valid
          // --------------+ sqrt  cos                   | -------------
          // 0. rnd -> log +-----+ sqrt                  |    T      F
          // 1. rnd -> sin   log +-----+                 |    T      F
          // 2.     +> cos   sin   log | sqrt * sin -> z |    F      T
          // 3.        sqrt  cos   sin | sqrt * cos -> z |    F      T
          // --------------+ sqrt  cos +-----------------|---------------
          // 4. rnd -> log +-----+ sqrt                  |    T      F
          // 5. rnd -> sin   log +-----+                 |    T      F
          // 6.        cos   sin   log | sqrt * sin -> z |    F      T
          // 7.        sqrt  cos   sin | sqrt * cos -> z |    F      T
          // 8.              sqrt  cos +-----------------|---------------
          //

          for(i <- 1 to n) {
            val xi = generateRandomUInt(rndW)
            val yi = generateRandomUInt(rndW)
            val (z1r, z2r) = reference(xi, yi)

            c.io.x.poke(xi.toBigInt.U(spec.W.W))
            val v1 = c.io.valid.peek().litValue.toBigInt
            val c1 = c.io.consume.peek().litValue.toBigInt
            assert(v1 == 0)
            assert(c1 == 1)
            c.clock.step(1)

            c.io.x.poke(yi.toBigInt.U(spec.W.W))
            val v2 = c.io.valid.peek().litValue.toBigInt
            val c2 = c.io.consume.peek().litValue.toBigInt
            assert(v2 == 0)
            assert(c2 == 1)
            c.clock.step(1)

            c.io.x.poke(0.toBigInt.U(spec.W.W))
            val v3 = c.io.valid.peek().litValue.toBigInt
            val c3 = c.io.consume.peek().litValue.toBigInt
            assert(v3 == 1)
            assert(c3 == 0)
            c.clock.step(1)

            c.io.x.poke(0.toBigInt.U(spec.W.W))
            val v4 = c.io.valid.peek().litValue.toBigInt
            val c4 = c.io.consume.peek().litValue.toBigInt
            assert(v4 == 1)
            assert(c4 == 0)
            c.clock.step(1)

            c.clock.step(1)
            c.clock.step(1)

            val z1 = c.io.z.peek().litValue.toBigInt
            c.clock.step(1)

            c.io.x.poke(0.toBigInt.U(spec.W.W))
            val z2 = c.io.z.peek().litValue.toBigInt
            c.clock.step(1)

            val simx  = ((SafeLong(1) << rndW) - xi).toDouble / (SafeLong(1) << rndW).toDouble
            val simy  = yi.toDouble / (SafeLong(1) << rndW).toDouble

            val siny = new RealGeneric(spec, sin(2.0 * Pi * simy))
            val cosy = new RealGeneric(spec, cos(2.0 * Pi * simy))
            val logx = new RealGeneric(spec, -2.0 * log(simx))
            val sqrx = new RealGeneric(spec, sqrt(-2.0 * log(simx)))

            println(f"sinSim  = (${siny.sgn}|${siny.exNorm}|${siny.man.toLong.toBinaryString})")
            println(f"cosSim  = (${cosy.sgn}|${cosy.exNorm}|${cosy.man.toLong.toBinaryString})")
            println(f"logSim  = (${logx.sgn}|${logx.exNorm}|${logx.man.toLong.toBinaryString})")
            println(f"sqrtSim = (${sqrx.sgn}|${sqrx.exNorm}|${sqrx.man.toLong.toBinaryString})")

            val z1sgn = bit(spec.W-1, z1).toInt
            val z1exp = slice(spec.manW, spec.exW, z1)
            val z1man = z1 & maskSL(spec.manW)

            val z2sgn = bit(spec.W-1, z2).toInt
            val z2exp = slice(spec.manW, spec.exW, z2)
            val z2man = z2 & maskSL(spec.manW)

            val z1g = new RealGeneric(spec, z1sgn, z1exp.toInt, z1man)
            val z2g = new RealGeneric(spec, z2sgn, z2exp.toInt, z2man)

            // TODO: fails when x << 1.
            // - if x is small, 1 - x is close to 1.
            // - if 1-x is close to 1, no shift is needed to normalize it.
            // - the table retuns a small value if 1-x is close to 1.
            // - after normalizing a small value from a table, it loses its precision.

            assert((z1 - z1r.value).abs <= tolerance,
                   f"x = ${xi.toDouble / (SafeLong(1) << rndW).toDouble}, y = ${yi.toDouble / (SafeLong(1) << rndW).toDouble}, ztest = ${z1g.toDouble}, zref = ${z1r.toDouble}, " +
                   f"test(${z1g.sgn}|${z1g.ex}(${z1g.exNorm})|${z1g.man.toLong.toBinaryString})(${z1g.toDouble}) !=" +
                   f" ref(${z1r.sgn}|${z1r.ex}(${z1r.exNorm})|${z1r.man.toLong.toBinaryString})(${z1r.toDouble})")
            assert((z2 - z2r.value).abs <= tolerance,
                   f"x = ${xi.toDouble / (SafeLong(1) << rndW).toDouble}, y = ${yi.toDouble / (SafeLong(1) << rndW).toDouble}, ztest = ${z2g.toDouble}, zref = ${z2r.toDouble}, " +
                   f"test(${z2g.sgn}|${z2g.ex}(${z2g.exNorm})|${z2g.man.toLong.toBinaryString})(${z2g.toDouble}) !=" +
                   f" ref(${z2r.sgn}|${z2r.ex}(${z2r.exNorm})|${z2r.man.toLong.toBinaryString})(${z2r.toDouble})")
          }
          println(f"Summary")
          if(maxError != 0.0) {
            val zexpected = if(satMaxError) {
              sqrt(-2.0 * log(xatMaxError)) * sin(2*Pi*yatMaxError)
            } else {
              sqrt(-2.0 * log(xatMaxError)) * cos(2*Pi*yatMaxError)
            }
            println(f"N=$n%d : largest errors ${maxError.toInt}%d where the value is "
                  + f"${zatMaxError} != ${zexpected}, diff = ${zatMaxError - zexpected}, "
                  + f"x = ${xatMaxError}, y = ${yatMaxError}")
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

  private def runChiSquared( spec: RealSpec, polySpec: PolynomialSpec, roundSpec: RoundSpec ) = {
    it should f"BoxMuller(x) spec ${spec.toStringShort} Chi^2 test" in {
      test( new BoxMuller(32, spec, polySpec, roundSpec)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          var zs = scala.collection.mutable.ArraySeq.empty[Double]

          // +/- 3 sigma ~ 99.7%.
          // +/- 4 sigma ~ 99.994%.
          // +/- 5 sigma ~ 99.99994%.
          val xmax   =  5.0
          val xmin   = -5.0
          val xrange = xmax - xmin

          // Since we use FP that is represented in the binary numeral system,
          // we need to split the range by a power of 2. Otherwise, the min and
          // max numbers of a bin does not align to the FP precision and the
          // "true" number of FP numbers in a bin differs between each other.
          // This cause significant deviation from the uniform distribution.
          val nbins = 16
          val dx    = xrange / nbins

          val ndeg      = nbins - 1
          val threshold = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          for(i <- 1 to n * 2) {

            c.io.x.poke(generateRandomUInt(32).toBigInt.U(32.W))
            val valid = c.io.valid.peek().litValue == 1
            val zi    = c.io.z.peek().litValue.toBigInt

            if (valid) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"outlier z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              zs = zs :+ zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          val cdf = (x: Double) => {
            0.5 * (1.0 + Erf.erf(x / sqrt(2.0)))
          }

          val nz = zs.length

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx
            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

            val nref = (cdf(maxRange) - cdf(minRange)) * nz

            val term = (nsamples - nref) * (nsamples - nref) / nref
            chi2 += term
          }
          chi2 /= ndeg

//           if(threshold < chi2) {
            println(f"-----------------------------------------------")
            println(f"${nz} real values generated")
            for(i <- 0 until nbins) {
              val minRange = xmin +  i    * dx
              val maxRange = xmin + (i+1) * dx
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))

              val nref = (cdf(maxRange) - cdf(minRange)) * nz
              println("n in [%8.3f, %8.3f) = %8d should be %10.3f".format(minRange, maxRange, nsamples, nref))
            }
            println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")
//           }

          assert(chi2 < threshold)
        }
      }
    }
  }

  val nOrderFP32    = 2
  val adrWFP32      = 8
  val extraBitsFP32 = 3
  val polySpecFP32  = new PolynomialSpec(RealSpec.Float32Spec, nOrderFP32, adrWFP32, extraBitsFP32)

  runTest(RealSpec.Float32Spec, polySpecFP32, RoundSpec.roundToEven, 32)
  runChiSquared(RealSpec.Float32Spec, polySpecFP32, RoundSpec.roundToEven)
}
