//
// Table Generation for general function
//

package rial.table

import scala.math._
import java.io._
import java.lang.Math.scalb
import chisel3._
import chisel3.util._
import org.apache.commons.math3._
import rial.util.ScalaUtil._
import rial.table._
import rial.table.NonbondFunctions._
import scala.reflect.runtime.universe
import scala.tools.reflect.ToolBox

private[rial] object GenGenericMain extends App {
  val allargs = args.mkString(" ")

  def compile[A](code: String): (Map[String, Any]) => A = {
    val tb = universe.runtimeMirror(getClass.getClassLoader).mkToolBox()
    val tree = tb.parse(
      s"""
       |import scala.math._
       |import org.apache.commons.math3._
       |def wrapper(context: Map[String, Any]): Any = {
       |  $code
       |}
       |wrapper _
      """.stripMargin)
    val f = tb.compile(tree)
    val wrapper = f()
    wrapper.asInstanceOf[Map[String, Any] => A]
  }

  case class Config(
    name:String="",
    func:String="",
    dir:String="",
    xmin:Double=0.0,
    xmax:Double=0.0,
    order:Int=3,
    bp:Int=24,
    imode:Int=0,
    omode:Int=0,
    signmode:Int=0,
    ndivlog:Int=0,
    unified:Boolean= false,
    cOutput:Boolean= false,
    verilog:Boolean= false
    )

  val parser = new scopt.OptionParser[Config]("GenGeneric") {
    head("scopt", "3.x")
    arg[String]("<name>").required() action { (x, c) =>
      c.copy(name = x) } text("name is the name of the function")
    arg[String]("<function>").required() action { (x, c) =>
      c.copy(func = x) } text("func is the equation for the function")

    opt[String]('d', "dir") action { (x, c) =>
      c.copy(name = x) } text("directory to put generated files")

    opt[Double]('l', "lower").required() action { (x, c) => c.copy(xmin = x) } text("input lower bound")
    opt[Double]('r', "upper").required() action { (x, c) => c.copy(xmax = x) } text("input upper bound")
    opt[Int]('o', "order").required() action { (x, c) => c.copy(order = x) } text("polynomial order")
    opt[Int]('b', "bp").required() action { (x, c) => c.copy(bp = x) } text("binary point of output")
    opt[Int]('w', "adrwidth").required() action { (x, c) => c.copy(ndivlog = x) } text("Address width w of table; number of segments = 2^w")
    opt[Int]('s', "imode") action { (x, c) =>
      c.copy(imode = x) } text("Input mode: 0: fixed, 1: floating, 2: floating with an additional segment [0,xmin)")
    opt[Int]('t', "omode") action { (x, c) => c.copy(omode = x) } text("Output mode: 0: fixed, 1: floating")
    opt[Int]('z', "signmode") action { (x, c) => c.copy(signmode = x)
    } text("Sign mode: 0: always include sign bit, 1: 2's compliment, delete sign bit 2: absolute value\n"+
      "if coefficients have both signs, a sign bit will be included (correspond to signmode==0).")
    opt[Unit]('u', "unified") action { (_, c) => c.copy(unified = true) } text("Concatenated output in verilog")
    opt[Unit]('c', "coutput") action { (_, c) => c.copy(cOutput = true) } text("Output C header (<name>.h)")
    opt[Unit]('v', "verilog") action { (_, c) => c.copy(verilog = true) } text("Output Verilog  (<name>.v)")
  }

  parser.parse(args, Config()) map { config =>
    println(f"name          = ${config.name}")
    println(f"dir           = ${config.dir}")
    println(f"function      = ${config.func}")
    println(f"order         = ${config.order}")
    println(f"range         = [${config.xmin},${config.xmax})")
    println(f"address width = ${config.ndivlog}, number of segments=${1<<config.ndivlog})")
    val imodeStr = Array("FIXED", "FLOATING", "FLOATING WITH SPECIAL FIXED REGION AROUND 0", "UNKNOWN")
    val imodeNorm = if ((config.imode>=0)||(config.imode<=2)) config.imode else 3
    println(f"input mode    = ${config.imode} : ${imodeStr(imodeNorm)}")
    val omodeStr = Array( "FIXED", "FLOATING", "UNKNOWN" )
    val omodeNorm = if ((config.omode>=0)||(config.omode<=1)) config.omode else 2
    println(f"output mode   = ${config.omode} : ${omodeStr(omodeNorm)}")
    val signmodeStr = Array( "Always include sign bit", "2's compliment, delete sign bit", "absolute value", "UNKNOWN" )
    val signmodeNorm = if ((config.signmode>=0)||(config.signmode<=2)) config.signmode else 3
    println(f"sign mode     = ${config.signmode} : ${signmodeStr(signmodeNorm)}")

    val code =
      s"""
       | val x=context("x").asInstanceOf[Double]
       | ${config.func}
      """
    val fref = compile[Double](code)

    def f ( x: Double ) : Double = {
      fref(Map("x" -> x))
    }

    //println(f(1.41421356))
    val ndiv    = 1<<config.ndivlog;
    val tbl = new FuncTableDouble( f, config.order )

    if (config.imode>0) {
      println("Input floating mode")
      println(f"x range = [${config.xmin},${config.xmax})")
      // for example:
      //  (xmin, xmax) = (0.5, 2.0) -> (eMin, eMax) = (-1, 1)
      //  (xmin, xmax) = (0.5, 2.1) -> (eMin, eMax) = (-1, 2)
      // eMax is actually eMax+1
      var eMin = log2DownD(config.xmin)
      val eMax = log2UpD(config.xmax)
      var de = eMax-eMin;
      if (config.imode==2) {
        println("imode==2: [0, xmin) is included as a fixed-point range")
        de+=1
      }
      val wexp = log2Up(de)
      val wman = config.ndivlog - wexp
      println(f"exp_min=${eMin} exp_max=${eMax} exp_bits=${wexp} man_bits=${wman}")
      // if 2^wexp > eMax-eMin, we can extend range for large or small values.
      // Here, we always extend for smaller values than xmin.
      if ((1<<wexp)>de) {
        eMin = eMax - (1<<wexp)
        if (config.imode==2) eMin += 1
        println(f"new exp_min=${eMin}")
      }
      val xmin_new = scalb(1.0, eMin)
      val xmax_new = scalb(1.0, eMax);
      println(f"new range = [${xmin_new},${xmax_new})")
      println(f"new range sqrt = [${sqrt(xmin_new)},${sqrt(xmax_new)})")

      val nblock = 1<<wman
      println(f"man_bits=${wman} nblock=${nblock}")
      var n=0
      if (config.imode==2) {
        // Add interval from 0~xmin
        val w = scalb(1.0, eMin)
        tbl.addRange(0.0, w, nblock)
        println(f"Block $n : xmin=0.0 xmax=${w} nblock=${nblock} width=${w/nblock}")
        n+=1
      }
      (eMin to eMax-1).toSeq.foreach( e => {
        val x = scalb(1.0, e)
        println(f"Block $n : xmin=$x xmax=${x+x} nblock=${nblock} width=${x/nblock}")
        tbl.addRange(x, x+x, nblock)
        n+=1
      } )
    } else {
      println("Input fixedpoint mode")
      println(f"x range = [${config.xmin},${config.xmax})")
      tbl.addRange(config.xmin, config.xmax, ndiv)
    }

    if (config.omode != 0) {
      println("Floating output");
      //mdg_table_integer_conversion_float(nbit, &t);
    } else {
      println("Fixed output");
      val tblI = new FuncTableInt( tbl, config.bp )
      if (config.verilog) {
        tblI.verilogOut( config.dir+config.name+".v", config.name, config.signmode, config.unified, allargs)
      }
      if (config.cOutput) {
        tblI.COut( config.dir+config.name+".h", config.name, allargs)
      }
    }
    // Plot

  }
}
