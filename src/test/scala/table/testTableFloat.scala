// See README.md for license details.

package rial.table

import chisel3._
import scopt.OptionParser
import scala.math._
//import scalax.chart.api._

import org.jfree.chart.ChartFactory
//import org.jfree.chart.ChartUtilities
import org.jfree.chart.JFreeChart
import org.jfree.chart.axis._
import org.jfree.chart.renderer.xy._
import org.jfree.chart.plot._
import org.jfree.data.general._
import org.jfree.data._
import org.jfree.data.xy._
import org.jfree.chart.ChartFrame
import org.jfree.chart.ChartPanel

import javax.swing.WindowConstants
import javax.swing.JFrame
import java.awt.GridLayout
import javax.swing.JPanel

object testTableFloat extends App {
  //Console.println("Hello World: " + (args mkString ", "))
  /*
  case class Config(ncycle:Int=1000)

  val parser = new scopt.OptionParser[Config]("RSquare") {
    head("scopt", "3.x")
    opt[Int]('n', "ncycle") action { (x, c) =>
      c.copy(ncycle = x) } text("ncycle is the number of simulation cycles")
  }

  parser.parse(arg1, Config()) map { config =>
    val ncycle = config.ncycle
    println(f"ncycle=$ncycle")
    iotesters.Driver.execute(arg0, () => new RSquare(28, 25, 0, 0) ) {
      c => new RSquareUnitTester(c, ncycle)
    }
  }
   */
  println("Test interpolation table, floating output")

  val targetFunc = sin(_)

  val xMin = 0.0
  val xMax = 0.5*Pi
  val nDiv = 64
  val bp   = 24
  val bpDx = 20

  val t = new FuncTableDouble( targetFunc, 3 )
  t.addRange(xMin, xMax, nDiv )
  val ti = new FuncTableInt( t, bp )

  //val calcW = ti.checkWidthAll( 20 )
  val cbit = Array(bp+2)
  val tf = new FuncTableFloat( t, true, cbit )

  val n = 100
  val xList = (0 to n-1).toList.map( y => xMin + y*(xMax-xMin)/n )

  //x.foreach(y => println(t.err(y), ti.eval(y, 20)))

  //ti.verilogOut("test.v", "test", 2, true)
  //ti.COut("test.h", "test")

  ////////////////////////////////////////////////////////////////////////
  // Absolute Error 
  var errset = new XYSeriesCollection
  var errSeries0 = new XYSeries("Error")
  xList.foreach( x => errSeries0.add( x, abs(tf.eval(x, bpDx)-targetFunc(x)) ) )

  var estimatedErrSeries = new XYSeries("Estimated Truncation Error")
  t.interval.foreach( iv => estimatedErrSeries.add( iv.xMin, iv.getEstimatedError() ) )

  errset.addSeries(errSeries0)
  errset.addSeries(estimatedErrSeries)

  var xAxis = new NumberAxis("x")
  var eAxis = new LogAxis("Error")
  eAxis.setBase(2)
  eAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits())
  var erenderer = new StandardXYItemRenderer
  var ePlot  = new XYPlot(errset, xAxis, eAxis, erenderer)
  var eChart = new JFreeChart("Absolute Error", ePlot)
  var ePanel = new ChartPanel(eChart)

  ////////////////////////////////////////////////////////////////////////
  // Relative Error
  val maxRelErr = 1e-3
  var relErrset = new XYSeriesCollection
  var relErrSeries = new XYSeries("Relative Error")
  xList.foreach( x => {
    val y0 = targetFunc(x)
    val err = if (y0 != 0.0) {
      val e = abs(tf.eval(x, bpDx)-y0)/y0
      if (e<maxRelErr) e
      else maxRelErr
    } else {
      maxRelErr
    }
    relErrSeries.add( x, err )
  } )

  relErrset.addSeries(relErrSeries)

  var rAxis = new LogAxis("rError")
  rAxis.setBase(2)
  rAxis.setStandardTickUnits(NumberAxis.createIntegerTickUnits())
  var rRenderer = new StandardXYItemRenderer
  var rPlot  = new XYPlot(relErrset, xAxis, rAxis, rRenderer)
  var rPanel = new ChartPanel(new JFreeChart("Relative Error", rPlot) )

  ////////////////////////////////////////////////////////////////////////
  // Caluculated value
  var dataset = new XYSeriesCollection
  var xySeries0 = new XYSeries("Interpolation")
  xList.foreach( x => xySeries0.add( x, tf.eval(x, bpDx) ) )
  var xySeries1 = new XYSeries("Target")
  xList.foreach( x => xySeries1.add( x, targetFunc(x) ) )

  dataset.addSeries(xySeries0)
  dataset.addSeries(xySeries1)

  var yAxis = new NumberAxis("y")
  var drenderer = new StandardXYItemRenderer
  var dPlot  = new XYPlot(dataset, xAxis, yAxis, drenderer)
  var dChart = new JFreeChart("Function value", dPlot)
  var dPanel = new ChartPanel(dChart)

  //plot.setRenderer(1, renderer)
  //plot.setDataset(1,  dataset)
  //plot.setRangeAxis(1, yAxis)
  //plot.mapDatasetToRangeAxis(1, 1)
  //var chart = new JFreeChart(plot)
  //var frame = new ChartFrame("sin", chart)

  var frame = new JFrame("sin")
  frame.setLayout(new GridLayout(1, 0))
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setBounds(10, 10, 900, 300)
  frame.add(dPanel)
  frame.add(ePanel)
  frame.add(rPanel)
  frame.setVisible(true)

  //val z = xList.map( y => (y, ti.eval(y, bpDx)) )
  //val chart = XYLineChart(z)
  //chart.show()
}

