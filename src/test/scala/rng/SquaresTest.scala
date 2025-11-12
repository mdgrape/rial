package rial.tests

import org.scalatest._

import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

import chisel3._
import chisel3.experimental.BundleLiterals._
import chisel3.util._

import chiseltest._

import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers
import org.scalatest.{BeforeAndAfterAllConfigMap, ConfigMap}

import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig
import rial.util.DebugControlMaster
import rial.rng._

class SquaresTest extends AnyFlatSpec
  with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Squares32"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  def squares(cnt: BigInt, key: BigInt): BigInt = {
    val mask64 = (BigInt(1) << 64) - BigInt(1)
    val mask32 = (BigInt(1) << 32) - BigInt(1)

    var x = cnt * key
    val y = cnt * key
    val z = y + key

    x = (((x * x) & mask64) + y) & mask64
    x = ((x >>32) & mask32) | ((x & mask32) << 32)

    x = (((x * x) & mask64) + z) & mask64
    x = ((x >>32) & mask32) | ((x & mask32) << 32)

    x = (((x * x) & mask64) + y) & mask64
    x = ((x >>32) & mask32) | ((x & mask32) << 32)

    x = (((x * x) & mask64) + z) & mask64
    x = ((x >>32) & mask32) | ((x & mask32) << 32)

    return x & mask32
  }

  def runTest(key: BigInt) = {
    test(new Squares32) {
      c => {
        c.io.en.poke(true.B)

        c.io.input.valid.poke(true.B)
        c.io.input.key.poke(key.U)

        for (i <- 0 to n-1) {
          c.io.input.count.poke(i.U)

          if(i >= 4) {
            assert(c.io.output.valid.peek().litValue == 1)
            assert(c.io.output.rand.peek().litValue == squares(BigInt(i-4), key))
          } else {
            assert(c.io.output.valid.peek().litValue == 0)
          }

          c.clock.step(1)
        }
      }
    }
  }

  it should "Squares should be equal to software impl" in {
    runTest(BigInt("2467cb532b5ce8d1", 16))
  }
}

class Squares32GeneratorTest extends AnyFlatSpec
  with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test Squares32Generator"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  def squares(cnt: BigInt, key: BigInt): BigInt = {
    val mask64 = (BigInt(1) << 64) - BigInt(1)
    val mask32 = (BigInt(1) << 32) - BigInt(1)

    var x = cnt * key
    val y = cnt * key
    val z = y + key

    x = (((x * x) & mask64) + y) & mask64
    x = ((x >>32) & mask32) | ((x & mask32) << 32)

    x = (((x * x) & mask64) + z) & mask64
    x = ((x >>32) & mask32) | ((x & mask32) << 32)

    x = (((x * x) & mask64) + y) & mask64
    x = ((x >>32) & mask32) | ((x & mask32) << 32)

    x = (((x * x) & mask64) + z) & mask64
    x = ((x >>32) & mask32) | ((x & mask32) << 32)

    return x & mask32
  }

  def runTest(key: BigInt) = {
    test(new Squares32Generator) {
      c => {
        c.io.initialized.expect(false.B)

        c.io.init.en   .poke(true.B)
        c.io.init.key  .poke(key.U)
        c.io.init.count.poke(  0.U)

        c.io.rand.ready.poke(false.B)

        c.clock.step(1) // --------------------------------------------

        c.io.initialized.expect(true.B)

        c.io.init.en   .poke(false.B)
        c.io.init.key  .poke(0.U)
        c.io.init.count.poke(0.U)

        c.io.rand.ready.poke(true.B)
        for (i <- 0 to n-1) {
          c.io.rand.valid.expect( (i >= c.nStage).B )

          if(c.io.rand.valid.peek().litValue == 1) {
            c.io.rand.bits.expect(squares(BigInt(i-c.nStage), key).U)
          }
          c.clock.step(1)
        }
      }
    }
  }

  it should "Squares should be equal to software impl" in {
    runTest(BigInt("2467cb532b5ce8d1", 16))
  }
}
