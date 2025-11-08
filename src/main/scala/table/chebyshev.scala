package rial.table

import scala.math._

private[rial] object Chebyshev {

  // Fit function by Chebyshev
  // size of c: norder+1
  def fit (f : Double => Double, xMin : Double, xMax: Double, nOrder : Int) = {
    val z = (0 to nOrder).toArray.map( x => Pi*(x+0.5)/(nOrder+1))
    val xz = z.map( x => (x, xMin + (cos(x)+1.0) * (xMax-xMin) *0.5) )
    val c = (0 to nOrder).toArray.map( j => xz.map( x => f(x._2)*cos(j*x._1) ).sum *2.0/(nOrder+1) )
    c.head*0.5 +: c.tail
  }

  // Fit function by Chebyshev
  // Divide xmin-xmax to ndiv region
  // size of c: ndiv * (norder+1)
  def fitDiv (f : Double => Double, xMin : Double, xMax: Double, nOrder : Int, nDiv : Int ) = {
    val w = (xMax-xMin)/nDiv
      (0 to nDiv-1).toArray.map( n => xMin+n*w ).map( l => fit(f, l, l+w, nOrder) )
  }

  // Convert Chebyshev coefficients to polynomial coefficients
  // T0 = 1, T1=x, Tn+1+Tn-1 = 2xTn
  // by Clenshaw's recurrence
  //    http://en.wikipedia.org/wiki/Clenshaw_algorithm
  def conv ( c : Array[Double] ) = {
    def step ( k : Double, c: Double, cb : ( Array[Double], Array[Double] ) ) : Array[Double] = {
      (c +: cb._1.map( x=> k*x ) ).zipAll(cb._2, 0.0, 0.0).map( x => x._1 - x._2 )
    }
    val initCB = ( Array(c.last), Array[Double]() )
    val cp = c.init.tail.foldRight(initCB)( (x, cb) => ( step(2.0, x, cb), cb._1 ) )
    step(1.0, c.head, cp)
  }

}

