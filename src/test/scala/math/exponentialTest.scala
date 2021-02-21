import org.scalatest._

import chisel3._
import chisel3.experimental.BundleLiterals._
import chiseltest._
//import org.scalatest.flatspec.AnyFlatSpec
//import org.scalatest.matchers.should.Matchers
//import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}
import org.scalatest.FlatSpec
import org.scalatest.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import spire.math.SafeLong
import spire.math.Numeric
import spire.implicits._

import rial.math._
import rial.arith._
import rial.table._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

//
// Testing Pow2F32 using ChiselTest
//

class Pow2F32Test extends FlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test pow2 F32"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

  def generateF32Within128 : Long = {
    val x  = ((r.nextDouble()-0.5)*128.0).toFloat
    val xi = java.lang.Float.floatToRawIntBits(x)
    xi.toLong&maskL(32)
  }

  def generateF32Within( p : Double, r : Random ) : Long = {
    val x  = ((r.nextDouble()-0.5)*p).toFloat
    val xi = java.lang.Float.floatToRawIntBits(x)
    xi.toLong&maskL(32)
  }

  def generateF32Full ( r : Random ) : Long = {
    val xi = r.nextInt
    xi.toLong&maskL(32)
  }

  private def runtest ( n : Int, stage : PipelineStageConfig ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"pipereg $pipeconfig" in {
      test(new Pow2F32(stage)) { c =>
        {
          var q  = new Queue[(Long,Long)]
          val nstage = c.getStage
          for (gen <- List( ("Test Within (-128,128)",generateF32Within(128.0,_)),
                            ("Test All range",generateF32Full _) ) ) {
            println(gen._1)
            for(i <- 1 to n+nstage) {
              val xi = gen._2(r)
              val xr = new RealGeneric( RealSpec.Float32Spec, SafeLong(xi) )
              val z0r= ExponentialSim.pow2F32Sim(xr)
              q += ((xi,z0r.value.toLong))
              c.io.x.poke(xi.U(32.W))
              val zi = c.io.z.peek.litValue.toLong
              if (i > nstage) {
                val (xid,z0d) = q.dequeue
                assert(zi == z0d, f"x=$xid%08x $zi%08x!=$z0d%08x")
              }
              c.clock.step(1)
            }
            q.clear
          }
        }
      }
    }
  }

  runtest(n, PipelineStageConfig.none())
  runtest(n, PipelineStageConfig.default(2))
}


class Pow2BF16Test extends FlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test pow2 BF16"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  private def runtest ( spec : RealSpec, stage : PipelineStageConfig,
    n : Int, r : Random,
    reference : (RealGeneric => RealGeneric),
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"pow(2,x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new Pow2Generic(spec,0,8,0, stage, false, false)) { c =>
        {
          val nstage = c.getStage
          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            val zi = c.io.z.peek.litValue.toBigInt
            if (i > nstage) {
              val (xid,z0d) = q.dequeue
              assert(zi == z0d, f"x=$xid%x $zi%x!=$z0d%x")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  def pow2TableGeneration( order : Int, adrW : Int, fracW : Int ) = {
    val tableD = new FuncTableDouble( x => pow(2.0,x)-1.0, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  val pow2BF16ExtraBits = 0
    val pow2BF16TableI = pow2TableGeneration( 0, 8, 7 )

  val pow2BF16sim = ExponentialSim.pow2simGeneric( pow2BF16TableI, 0, _ )

  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(),
    n, r, pow2BF16sim,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none(),
    n, r, pow2BF16sim,
    "Test All range",generateRealFull(_,_) )

}

class ExpF32Test extends FlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test exp F32"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(19660809)

  def generateRealWithin( p : Double, spec: RealSpec, r : Random ) = {
    val rD : Double = (r.nextDouble()-0.5)*p
    val x = new RealGeneric(spec, rD)
    new RealGeneric (spec, (x.value & (maskSL(spec.exW+1)<<spec.manW)) + SafeLong(BigInt(spec.manW, r)))
  }

  def generateRealFull ( spec: RealSpec, r : Random ) = {
    new RealGeneric (spec, SafeLong(BigInt(spec.W, r)))
  }

  def errorLSB( x : RealGeneric, y : Double ) : Double = {
    val err = x.toDouble - y
    java.lang.Math.scalb(err, -x.exNorm+x.spec.manW)
  }

  private def runtest ( spec : RealSpec, stage : PipelineStageConfig,
    n : Int, r : Random,
    reference : (RealGeneric => RealGeneric),
    generatorStr : String, generator : ( (RealSpec, Random) => RealGeneric) ) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"pow(2,x) pipereg $pipeconfig spec ${spec.toStringShort} $generatorStr " in {
      test( new ExponentialGeneric(spec,2,8,2,stage, false, false)) { c =>
        {
          val nstage = c.getStage
          val q  = new Queue[(BigInt,BigInt)]
          for(i <- 1 to n+nstage) {
            val xi = generator(spec,r)
            val z0r= reference(xi)
            q += ((xi.value.toBigInt,z0r.value.toBigInt))
            c.io.x.poke(xi.value.toBigInt.U(spec.W.W))
            val zi = c.io.z.peek.litValue.toBigInt
            if (i > nstage) {
              val (xid,z0d) = q.dequeue
              assert(zi == z0d, f"x=$xid%x $zi%x!=$z0d%x")
            }
            c.clock.step(1)
          }
        }
      }
    }
  }

  def pow2TableGeneration( order : Int, adrW : Int, fracW : Int ) = {
    val tableD = new FuncTableDouble( x => pow(2.0,x)-1.0, order )
    tableD.addRange(0.0, 1.0, 1<<adrW)
    new FuncTableInt( tableD, fracW )
  }

  val pow2F32ExtraBits = 0
  val pow2F32TableI = pow2TableGeneration( 2, 8, 23 )
  val pow2F32sim = ExponentialSim.expSimGeneric( pow2F32TableI, 2, _ )

  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
    n, r, pow2F32sim,
    "Test Within (-128,128)",generateRealWithin(128.0,_,_))
  runtest(RealSpec.Float32Spec, PipelineStageConfig.none(),
    n, r, pow2F32sim,
    "Test All range",generateRealFull(_,_) )
}

