package rial.tests

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
import rial.math._
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
    xi.toLong&mask(32)
  }

  def generateF32Within( p : Double, r : Random ) : Long = {
    val x  = ((r.nextDouble()-0.5)*p).toFloat
    val xi = java.lang.Float.floatToRawIntBits(x)
    xi.toLong&mask(32)
  }

  def generateF32Full ( r : Random ) : Long = {
    val xi = r.nextInt
    xi.toLong&mask(32)
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
              //println(xi)
              val z0i= ExponentialSim.pow2F32Sim(xi)
              q += ((xi,z0i))
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

