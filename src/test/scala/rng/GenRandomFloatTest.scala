package rial.tests

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
import rial.rng._
import rial.util.ScalaUtil._
import rial.util.PipelineStageConfig

import scala.util.Random
import scala.math._
import scala.collection.mutable.Queue
import scala.language.reflectiveCalls

//
// Test GenRandomFloat12 using chi-squared test.
//

class Uniform12Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test uniform real distribution [1, 2)"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRandomUInt(width: Int) = {
    val rnd = SafeLong(r.nextLong)
    rnd & maskSL(width)
  }

  private def runtest ( spec : RealSpec, stage: PipelineStageConfig) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"uniform12(x) pipereg $pipeconfig spec ${spec.toStringShort}" in {
      test( new GenRandomFloat12(32, spec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val nstage = stage.total
          val nnum = c.nRndInt
          var zs = scala.collection.mutable.ArraySeq.empty[Double]

          val xmax   = 2.0
          val xmin   = 1.0
          val xrange = xmax - xmin

          // Since we use FP that is represented in the binary numeral system,
          // we need to split the range by a power of 2. Otherwise, the min and
          // max numbers of a bin does not align to the FP precision and the
          // "true" number of FP numbers in a bin differs between each other.
          // This cause significant deviation from the uniform distribution.
          val nbins = 16
          val dx    = xrange / nbins
          val nref  = n.toDouble / nbins

          val ndeg      = nbins - 1
          val threshold = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          for(i <- 1 to n+nstage) {
            for(j <- 0 until nnum) {
              c.io.rnds(j).poke(generateRandomUInt(32).toBigInt.U(32.W))
            }
            val zi = c.io.z.peek().litValue.toBigInt

            if (i > nstage) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              assert(xmin <= zd.toDouble && zd.toDouble < xmax)

              zs = zs :+ zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx

            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

            val term = (nsamples - nref) * (nsamples - nref) / nref
            chi2 += term
          }
          chi2 /= ndeg

          if(threshold < chi2) {
            println(f"-----------------------------------------------")
            for(i <- 0 until nbins) {
              val minRange = xmin +  i    * dx
              val maxRange = xmin + (i+1) * dx
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))
              println(f"n in [${minRange}%6f, ${maxRange}%6f) = ${nsamples} should be ${nref}")
            }
            println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")
          }

          assert(chi2 < threshold)
        }
      }
    }
  }

  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none)
  runtest(RealSpec.Float32Spec,  PipelineStageConfig.none)
  runtest(RealSpec.Float64Spec,  PipelineStageConfig.none)
}


//
// Test GenRandomFloat01 using chi-squared test.
//

class Uniform01Test extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test uniform real distribution [0, 1)"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRandomUInt(width: Int) = {
    val rnd = SafeLong(r.nextLong)
    rnd & maskSL(width)
  }

  private def runtest ( spec : RealSpec, stage: PipelineStageConfig) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"uniform01(x) pipereg $pipeconfig spec ${spec.toStringShort}" in {
      test( new GenRandomFloat01(32, spec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val nstage = stage.total
          val nnum = c.nRndInt
          var zs = scala.collection.mutable.ArraySeq.empty[Double]

          val xmax   = 1.0
          val xmin   = 0.0
          val xrange = xmax - xmin

          // Since we use FP that is represented in the binary numeral system,
          // we need to split the range by a power of 2. Otherwise, the min and
          // max numbers of a bin does not align to the FP precision and the
          // "true" number of FP numbers in a bin differs between each other.
          // This cause significant deviation from the uniform distribution.
          val nbins = 16
          val dx    = xrange / nbins
          val nref  = n.toDouble / nbins

          val ndeg      = nbins - 1
          val threshold = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          for(i <- 1 to n+nstage) {
            for(j <- 0 until nnum) {
              c.io.rnds(j).poke(generateRandomUInt(32).toBigInt.U(32.W))
            }
            val zi = c.io.z.peek().litValue.toBigInt

            if (i > nstage) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              assert(xmin <= zd.toDouble && zd.toDouble < xmax)

              zs = zs :+ zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx

            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

            val term = (nsamples - nref) * (nsamples - nref) / nref
            chi2 += term
          }
          chi2 /= ndeg

          if(threshold < chi2) {
            println(f"-----------------------------------------------")
            for(i <- 0 until nbins) {
              val minRange = xmin +  i    * dx
              val maxRange = xmin + (i+1) * dx
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))
              println(f"n in [${minRange}%6f, ${maxRange}%6f) = ${nsamples} should be ${nref}")
            }
            println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")
          }

          assert(chi2 < threshold)
        }
      }
    }
  }

  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none)
  runtest(RealSpec.Float32Spec,  PipelineStageConfig.none)
  runtest(RealSpec.Float64Spec,  PipelineStageConfig.none)
}


//
// Test GenRandomFloat01Full using chi-squared test.
//

class Uniform01FullTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test uniform real distribution [0, 1)"

  var n = 1000

  override def beforeAll(configMap: ConfigMap) = {
    n = configMap.getOptional[String]("n").getOrElse("1000").toInt
    println(s"ncycle=$n")
  }

  val r = new Random(123456789)

  def generateRandomUInt(width: Int) = {
    val rnd = SafeLong(r.nextLong)
    rnd & maskSL(width)
  }

  private def runtest ( spec : RealSpec, stage: PipelineStageConfig) = {
    val total = stage.total
    val pipeconfig = stage.getString
    it should f"uniform01Full(x) pipereg $pipeconfig spec ${spec.toStringShort}" in {
      test( new GenRandomFloat01Full(32, spec, stage)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          val nnum = c.nRndInt

          // check if it works with rnd === 0.U first
          var zerochecked = false
          for(i <- 1 until 2) {
            for(j <- 0 until nnum) {
              c.io.rnds(j).poke(0.U(32.W))
            }
            val z = c.io.z.peek().litValue.toBigInt
            assert(z == 0, "if all the rnds are zero, z should also be zero")
            zerochecked = true
            c.clock.step(1)
          }
          assert(zerochecked, "zero case should be checked")

          val nstage = stage.total
          var zs = scala.collection.mutable.ArraySeq.empty[Double]

          val xmax   = 1.0
          val xmin   = 0.0
          val xrange = xmax - xmin

          // Since we use FP that is represented in the binary numeral system,
          // we need to split the range by a power of 2. Otherwise, the min and
          // max numbers of a bin does not align to the FP precision and the
          // "true" number of FP numbers in a bin differs between each other.
          // This cause significant deviation from the uniform distribution.
          val nbins = 16
          val dx    = xrange / nbins
          val nref  = n.toDouble / nbins

          val ndeg      = nbins - 1
          val threshold = 1.8 // probability of getting this value is ~2.9%

          // generate random numbers

          for(i <- 1 to n+nstage) {
            for(j <- 0 until nnum) {
              c.io.rnds(j).poke(generateRandomUInt(32).toBigInt.U(32.W))
            }
            val zi = c.io.z.peek().litValue.toBigInt

            if (i > nstage) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              assert(xmin <= zd.toDouble && zd.toDouble < xmax)

              zs = zs :+ zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx

            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

            val term = (nsamples - nref) * (nsamples - nref) / nref
            chi2 += term
          }
          chi2 /= ndeg

          if(threshold < chi2) {
            println(f"-----------------------------------------------")
            for(i <- 0 until nbins) {
              val minRange = xmin +  i    * dx
              val maxRange = xmin + (i+1) * dx
              val nsamples = zs.count(a => (minRange <= a && a < maxRange))
              println(f"n in [${minRange}%6f, ${maxRange}%6f) = ${nsamples} should be ${nref}")
            }
            println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")
          }

          assert(chi2 < threshold)
        }
      }
    }
  }

  runtest(RealSpec.BFloat16Spec, PipelineStageConfig.none)
  runtest(RealSpec.Float32Spec,  PipelineStageConfig.none)
  runtest(RealSpec.Float64Spec,  PipelineStageConfig.none)
}


//
// Test GenRandomFloat01Full using chi-squared test.
//

class Uniform01FullPipeTest extends AnyFlatSpec
    with ChiselScalatestTester with Matchers with BeforeAndAfterAllConfigMap {

  behavior of "Test uniform real distribution [0, 1)"

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

  private def runtest ( spec : RealSpec) = {
    it should f"uniform01FullPipe(x) spec ${spec.toStringShort}" in {
      test( new GenRandomFloat01FullPipe(32, spec)).
        withAnnotations(Seq(VerilatorBackendAnnotation)) { c =>
        {
          // check if it works with rnd === 0.U first
          var zerochecked = false
          println(f"nRndsForEx  = ${c.nRndsForEx }")
          println(f"nRndsForMan = ${c.nRndsForMan}")
          for(i <- 1 to (c.nRndsForEx + c.nRndsForMan) * 2) {
            c.io.rnd.poke(0.U(32.W))
            c.clock.step(1)
            val valid = c.io.valid.peek().litValue == 1
            val z     = c.io.z.peek().litValue.toBigInt
            if(valid) {
              assert(z == 0, "if all the rnds are zero, z should also be zero")
              zerochecked = true
            }
          }
          assert(zerochecked, "zero case should be checked")

          // check if it generates the correct uniform distribution

          var zs = scala.collection.mutable.ArraySeq.empty[Double]

          val xmax   = 1.0
          val xmin   = 0.0
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

          for(i <- 1 to n) {

            c.io.rnd.poke(generateRandomUInt(32).toBigInt.U(32.W))
            val valid = c.io.valid.peek().litValue == 1
            val zi    = c.io.z.peek().litValue.toBigInt

            if (valid) {
              val zd = new RealGeneric(spec, zi)

              if(zd.toDouble < xmin || xmax <= zd.toDouble){
                val zisgn = bit(spec.W-1, zi).toInt
                val ziexp = slice(spec.manW, spec.exW, zi)
                val ziman = zi & maskSL(spec.manW)
                println(f"z = (${zisgn}|${ziexp}(${ziexp - spec.exBias})|${ziman.toLong.toBinaryString}) = ${zd.toDouble}")
              }
              assert(xmin <= zd.toDouble && zd.toDouble < xmax)

              zs = zs :+ zd.toDouble
            }
            c.clock.step(1)
          }

          // calculate chi^2

          val nz    = zs.length
          val nref  = nz.toDouble / nbins

          var chi2 = 0.0
          for(i <- 0 until nbins) {
            val minRange = xmin +  i    * dx
            val maxRange = xmin + (i+1) * dx

            val nsamples = zs.count(a => (minRange <= a && a < maxRange))

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
              println(f"n in [${minRange}%6f, ${maxRange}%6f) = ${nsamples} should be ${nref}")
            }
            println(f"chi^2 = ${chi2}, threshold(<3%%) = ${threshold}")
//           }

          val zmin = new RealGeneric(spec, zs.min)
          println(f"min z = ${zs.min}(0|${zmin.ex-spec.exBias}|${zmin.man.toLong.toBinaryString})")

          assert(chi2 < threshold)
        }
      }
    }
  }

  runtest(RealSpec.BFloat16Spec)
  runtest(RealSpec.Float32Spec)
  runtest(RealSpec.Float64Spec)
}
