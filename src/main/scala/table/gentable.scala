//
// Libraries for table generation
//   Copyright (C) 2019 Makoto Taiji, RIKEN BDR
//
package rial.table

import scala.math._
import java.io._
import java.lang.Math.scalb
import chisel3._
import chisel3.util._
import rial.util.ScalaUtil._

class FuncTableIntervalDouble (
  f : Double => Double,
  val nOrder: Int,
  val xMin  : Double,
  val w     : Double
) {
  val c  = if (nOrder==0) {
    Array(f(xMin))
  } else {
    Chebyshev.conv(Chebyshev.fit(f, xMin, xMin+w, nOrder))
  }
  val c1 = Chebyshev.conv(Chebyshev.fit(f, xMin, xMin+w, nOrder+1))

  def xMax = xMin+w
  def inRange( x : Double ) = (xMin <= x) && (x < xMin+w)

  def evalNorm(d : Double, skip : Int = 0) = {
    c.drop(skip).foldRight(0.0)( (cj, z) => z*d+cj )
  }

  def eval(x : Double, skip : Int = 0) = {
    val d = 2.0*(x - xMin)/w - 1.0
    c.drop(skip).foldRight(0.0)( (cj, z) => z*d+cj )
  }

  def derivCoeff( order : Int ) = {
    (1 to order).toSeq.zip(c.drop(nOrder+1-order)).map( x => ( x._1 * x._2 ) )
  }

  def evalDerivNorm(d : Double, order : Int) = {
    val dc = derivCoeff( order )
    dc.foldRight(0.0)( (cj, z) => z*d+cj )
  }

  def coeff(n : Int) = { c(n) }

  def getMinMax(order : Int) = {
    def deriv( x : Double ) = { 
      evalDerivNorm( x, order)
    }
    def func( x : Double ) = { 
      evalNorm( x, nOrder-order)
    }
    if (nOrder==0) {
      ( abs(c(0)), abs(c(0)) )
    } else if (order==0) {
      ( abs(c.last), abs(c.last) )
    } else {
      estimateMinMax.getMinMaxDouble( func, deriv, -1.0, 1.0 )
    } 
  }

  def getAbsMax(order : Int) = {
    def deriv( x : Double ) = { 
      evalDerivNorm( x, order)
    }
    def func( x : Double ) = { 
      evalNorm( x, nOrder-order)
    }
    if (nOrder==0) {
      abs(c(0))
    } else if (order==0) {
      abs(c.last)
    } else {
      val minmax = estimateMinMax.getMinMaxDouble( func, deriv, -1.0, 1.0 )
      max(abs(minmax._1),abs(minmax._2))
    } 
  }

  // Calculate relative scale to final value
  def getRelativeScale(order : Int) = {
    val l2absMax0 = log2UpD(getAbsMax(0))
    log2UpD(c(order))-l2absMax0
  }

  def getEstimatedError() = {
    c1(nOrder+1)
  }
}

class FuncTableDouble (f : Double => Double, val nOrder : Int) {
  var interval = new Array[FuncTableIntervalDouble](0)

  def nInterval = interval.length

  def adrW = log2Up(interval.length)

  def add(xMin : Double, w : Double) = {
    interval = interval :+ new FuncTableIntervalDouble(f, nOrder, xMin, w)
  }

  def addN(xMin : Double, w : Double, n : Int) = {
    (0 to n-1).toSeq.foreach( x => add(xMin+x*w, w) )
  }

  def addRange(xMin : Double, xMax : Double, n : Int) = {
    val w = (xMax-xMin)/n
    (0 to n-1).toSeq.foreach( x => add(xMin+x*w, w) )
  }

  def eval(x : Double) = {
    interval.find( t => t.inRange(x) ) match {
      case Some(iv) => iv.eval(x)
      case None     => { println(f"${this.getClass.getName}.eval: Range Error"); 0.0 }
    }
  }

  def err(x : Double) = {
    val z  = eval(x)
    val z0 = f(x)
    (z, z0, z-z0)
  }

  // Return sign information
  // 0 : All 0
  // 1 : All 0 or positive
  // 2 : All 0 or negative
  // 3 : include positive and negative
  def checkSign() = {
    (0 to nOrder).toSeq.
      map( i =>
        interval.foldLeft(0)(
          (s, iv) =>
          if      (iv.c(i)<0.0) s|2
          else if (iv.c(i)>0.0) s|1
          else s
        ) )
  }

  // Make a list of coefficients for specified order
  def coeffList( n : Int ) = {
    interval.map( iv => iv.c(n) )
  }

  def minMax() = {
    val cmin = (0 to nOrder).toSeq.map( i => coeffList(i).min )
    val cmax = (0 to nOrder).toSeq.map( i => coeffList(i).max )
    cmin.zip(cmax)
  }

  def getMinMaxAll() = {
    (nOrder to 0 by -1).toSeq.map(
      n => interval.foldLeft((1e16,0.0))( (w, iv) => {
        val mm = iv.getMinMax(n)
        ( min(mm._1, w._1), max(mm._2, w._2) )
      } )
    )
  }

  def getAbsMaxAll() = {
    (nOrder to 0 by -1).toSeq.map(
      n => interval.foldLeft(0.0)( (w, iv) => max(w, iv.getAbsMax(n)) )
    )
  }

}

class FuncTableIntervalInt (iv : FuncTableIntervalDouble, val floating: Boolean, val cbit : Seq[Int], val bp : Int = 0 ) {

  //c.foreach(x => println(f"${x._1}%h"))
  val xMin   = iv.xMin
  val w      = iv.w
  val nOrder = iv.nOrder

  def coeff() = cw.map( _._1 )
  def coeff(n: Int) = cw(n)._1

  def width() = cw.map( _._2 )
  def width(n: Int) = cw(n)._2

  def derivCoeff( order : Int ) = {
    val n = nOrder
    (1 to order).toSeq.zip(cw.drop(n+1-order)).map( x => ( x._1 * x._2._1, x._2._2 ) )
  }

  def horner( c: (Long, Int), z: (Long, Int), dx: Long, dxBp: Int, checkOverflow: Boolean ) = {
    val zw   = z._2
    val d    = c._2
    val dxw  = min(zw, dxBp)
    val dxl  = dx >> (dxBp - dxw)// MSBs?
    val prod = dxl * z._1
    val prod_sft = prod >> dxw
    val term = c._1 + prod_sft

    if (checkOverflow) {
      // Check width
      val zmax = (1L<<(d-1))-1
      val zmin = -zmax-1
      if ( (term > zmax) || (term < zmin) ) {
        println(f"INFO (${this.getClass.getName}) : Range error: dx=$dx%f, dx=$dx%h, z=${term}%h(${term.toDouble}%e), c=${c._1}%h, cmin=$zmin%h, zmax=$zmax%h")
        //println(f"${c._1} ${z._1} ${term} $zw $d $dxBp $dx%h $dxl")
      }
      //println(f"${c._1}%x ${z._1}%x ${term}%x $zw $d $dxBp $dx%x $dxl")
    }
//     println("sim----------------------------------------")
//     println(f"c     : ${c._1.toBinaryString}(${c._1})")
//     println(f"zw    : ${zw  .toBinaryString}(${zw  })")
//     println(f"dx    : ${dx.toBinaryString}(${dx})")
//     println(f"dxBp  : ${dxBp.toBinaryString}(${dxBp})")
//     println(f"dxw   : ${dxw                }         ")
//     println(f"dxl   : ${dxl                }         ")
//     println(f"z     : ${z                  }         ")
//     println(f"prod  : ${prod}")
//     println(f"prodsf: ${prod_sft}")
//     println(f"term  : ${term.toBinaryString}(${term})")
    (term, d)
  }

  def eval ( dx : Long, dxBp : Int, checkOverflow: Boolean = true, skip : Int = 0 ) = {
    // Note: in Scala and C, a result of division is rounded toward 0.
    // That's why I use double for this calculation..
    val res = cw.init.drop(skip).foldRight(cw.last)( (cj,z) => horner( cj, z, dx, dxBp, checkOverflow ) )
    res._1.toLong
  }

  def evalDeriv ( dx : Long, dxBp : Int, order : Int = 1 ) = {
    val dc = derivCoeff(order)
    val res = dc.init.foldRight(dc.last)( (cj,z) => horner( cj, z, dx, dxBp, false ) )
    res._1.toLong
  }

  def inRange( x : Double ) = (xMin <= x) && (x < xMin+w)

  def evalD ( x : Double,  dxBP : Int, checkOverflow: Boolean = true) = {
    val d = 2.0*(x-xMin)/w - 1.0
    if ((d < -1.0)||(d>=1.0)) {
      println(f"ERROR (${this.getClass.getName}) : Range error: x=$x%f, xMin=${xMin}%f, w=${w}%f dx=$d%f")
      0.0
    } else { 
      val dx : Long = scalb(d, dxBP).toLong
      val z = eval( dx, dxBP, checkOverflow ).toDouble
      if (floating) {
        scalb( z, exponent-(cbit(0)-2) )
      } else {
        scalb( z, -bp )
      }
    }
  }

  def checkWidth( dxBp : Int, order : Int) = {
    def deriv( x : Long ) = { 
      evalDeriv( x, dxBp, order)
    }
    def func( x : Long ) = { 
      eval( x, dxBp, true, nOrder-order)
    }
    var width  = cbit(nOrder-order)
    var maxV   = 1L<<(width-1)

    val xMax = 1L<<dxBp

    val minmax = estimateMinMax.getMinMaxLong( func, deriv, -xMax, xMax-1 )
    val mn = minmax._1
    val mx = minmax._2
    //if (minmax._1 * minmax._2)
    if ( (mn < -maxV) || (mx >= maxV) ) {
      println(f"${this.getClass.getName} : Overflow at ${iv.xMin}, min=${mn}(${mn.toDouble/maxV}), max=${mx.toDouble}(${mx/maxV}), limit=$maxV")
      while ( (mn < -maxV) || (mx >= maxV) ) {
        width += 1
        maxV *= 2
      }
    }
    width
  }

  def checkWidthAll( dxBp : Int ) = {
    (nOrder to 1 by -1).toSeq.map( n => checkWidth(dxBp, n) )
  }

  def scaleCoeff( c : Seq[Double], cbit : Seq[Int] ) = {
    // Tune coefficients to fit given width
    val exc = iv.c.map( x => log2UpD(abs(x)) )
    val dif = exc.zip(cbit).map( x => x._1 - (x._2-2) ) // cbit include sign bits
    val maxVal = dif.max
    val maxPos = dif.indexWhere(x=> x == maxVal)
    // Scale to fit all coefficients within cbit including sign
    (maxPos, maxVal)
  }

  val (cw, exponent) = if (floating) {
    val (maxPos, scaleBits) = scaleCoeff(iv.c, cbit)
    // When c0=1, exp=0 and maxVal = -(cbit(0)-2)
    ( iv.c.map( x => round(scalb(x, -scaleBits)).toLong ).zip(cbit),
      scaleBits + cbit(0) - 2 )
  } else {
    (iv.c.map( x => round(scalb(x, bp)).toLong ).zip(cbit), 0)
  }
  
}

trait FuncTable {
  val nOrder : Int
  val sign : Seq[Int]
  val cbit : Seq[Int]
  val calcWidth : Seq[Int]

  val interval : Array[FuncTableIntervalInt]

  def nInterval = interval.length
  def adrW = log2Up(interval.length)

  //   signMode=0 : always include sign bit
  //   signMode=1 : 2's complement and no sign bit (if possible)
  //   signMode=2 : absolute and no sign bit (if possible)
  def getVectorUnified( signMode : Int = 0 ) = {
    val nAdr = log2Up(interval.size);
    val w = cbit.zip(sign).map( x => if ((x._2 != 3) && (signMode!=0)) x._1-1 else x._1 )
//     println(f"getVectorUnified: cbit = ${cbit}")
//     println(f"getVectorUnified: w    = ${w}")
    val s = sign.map( x =>
      if ((x == 2) && (signMode!=0)) -1
      else if ((x == 1) && (signMode!=0)) 1
      else 0
    )
    val takeAbs = sign.map( x => (x!=3) && (signMode==2) )
    println(f"//  sign Mode : $signMode")
    println(f"//  sign  : "+s.mkString("", ", ", ""))
    println(f"//  width : "+w.mkString("", ", ", ""))
    println(f"//  absolute value : "+ takeAbs.mkString("", ", ", ""))
    val wTotal = w.sum
    // combined, Vector of Vector
    val cmask = interval.map( iv =>
      iv.coeff.zip(takeAbs).map( x => if (x._2) abs(x._1) else x._1).
        zip(w).map( x => toBinStringFill( x._1 & ((1L<<x._2)-1), x._2)).mkString("b","_","")
    )
    //cmask.foreach( x=> println(x) )
    val z = VecInit( cmask.toArray.map( x => x.U(wTotal.W) ) )
    (z,w)
  }

  def getCBitWidth(signMode : Int = 0 ) = {
    cbit.zip(sign).map( x => if ((x._2 != 3) && (signMode!=0)) x._1-1 else x._1 )
  }
  def getVectorWithWidth(w : Seq[Int], signMode : Int = 0 ) = {
    val nAdr = log2Up(interval.size);
    val s = sign.map( x =>
      if ((x == 2) && (signMode!=0)) -1
      else if ((x == 1) && (signMode!=0)) 1
      else 0
    )
    val takeAbs = sign.map( x => (x!=3) && (signMode==2) )
    val wTotal = w.sum
    // combined, Vector of Vector
    val cmask = interval.map( iv =>
      iv.coeff.zip(takeAbs).map( x => if (x._2) abs(x._1) else x._1).
        zip(w).map( x => toBinStringFill( x._1 & ((1L<<x._2)-1), x._2)).mkString("b","_","")
    )
    //cmask.foreach( x=> println(x) )
    val z = VecInit( cmask.toArray.map( x => x.U(wTotal.W) ) )
    z
  }

  def getVectorSeparate( signMode : Int = 0) = {
    val nAdr = log2Up(interval.size);
    val w = cbit.zip(sign).map( x => if ((x._2 != 3) && (signMode!=0)) x._1-1 else x._1 )
    val s = sign.map( x =>
      if ((x == 2) && (signMode!=0)) -1
      else if ((x == 1) && (signMode!=0)) 1
      else 0
    )
    val takeAbs = sign.map( x => (x!=3) && (signMode==2) )
    println(f"//  sign Mode : $signMode")
    println(f"//  sign  : "+s.mkString("", ", ", ""))
    println(f"//  width : "+w.mkString("", ", ", ""))
    println(f"//  absolute value : "+ takeAbs.mkString("", ", ", ""))
    // separate array
    (0 to nOrder).toSeq.map( i => {
      val mask = (1L<<w(i))-1
      val cmask = interval.map( z => {
        val c = z.coeff.apply(i)
        val ca = if (takeAbs(i)) abs(c) else c
        val cm = (ca & mask).U(w(i).W)
        cm
      } )
      val z = VecInit( cmask.toArray )
      (z, w(i))
    } )
  }
  
  //   signMode=0 : always include sign bit
  //   signMode=1 : 2's complement and no sign bit (if possible)
  //   signMode=2 : absolute and no sign bit (if possible)
  def verilogOut( filename : String, moduleName : String, signMode : Int = 0, unified : Boolean = false, cmdline : String = "") = {
    val p = new PrintWriter(filename)
    try {
      val nAdr = log2Up(interval.size);
      val w = cbit.zip(sign).map( x => if ((x._2 != 3) && (signMode!=0)) x._1-1 else x._1 )
      val s = sign.map( x =>
        if ((x == 2) && (signMode!=0)) -1
        else if ((x == 1) && (signMode!=0)) 1
        else 0
      )
      val takeAbs = sign.map( x => (x!=3) && (signMode==2) )

      if (unified) {
        val wTotal = w.sum
        p.println(f"// $moduleName");
        p.println(f"// $cmdline");
        p.println(f"//  sign Mode : $signMode")
        p.println(f"//  sign  : "+s.mkString("", ", ", ""))
        p.println(f"//  width : "+w.mkString("", ", ", ""))
        p.println(f"//  absolute value : "+ takeAbs.mkString("", ", ", ""))
        p.println(f"module ${moduleName} (input [${nAdr}-1:0] adr, output reg [$wTotal-1:0] c);")
        p.println(f"  always_comb begin");
        p.println(f"    case (adr)");
        for ( j <- 0 to interval.size-1 ) {
          val c  = interval(j).coeff
          val ca = c.zip(takeAbs).map( x => if (x._2) abs(x._1) else x._1)
          val cmask = ca.zip(w).map( x => x._1 & ((1L<<x._2)-1))
          val cmaskStr = cmask.zip(w).map( x => toBinStringFill(x._1, x._2) )
          p.println(f"      $nAdr'h$j%h: c = $wTotal'b"+cmaskStr.mkString("", "_", ";"))
          for ( i <- 0 to nOrder ) {
            val maxV = 1L<<(cbit(i)-1)
            if ( ( ca(i) >= maxV ) || ( ca(i) < -maxV ) ) {
              p.println(f"    !!!RANGE ERROR in $i !!!");
            }
          }
        }
        p.println(f"      default: c = $wTotal'bx;")
        p.println(f"    endcase;");
        p.println(f"  end");
        p.println(f"endmodule // test\n")
      } else {
        for ( i <- 0 to nOrder ) {
          val maxV = 1L<<(cbit(i)-1)
          val mask = (1L<<w(i))-1
          p.println(f"//  sign : ${s(i)}")
          p.println(f"//  absolute value : ${takeAbs(i)}")
          p.println(f"module ${moduleName}_$i (input [${nAdr}-1:0] adr, output reg [$w(i)-1:0] c);")
          p.println(f"  always_comb begin");
          p.println(f"    case (adr)");
          for ( j <- 0 to interval.size-1 ) {
            val c  = interval(j).coeff(i)
            val ca = if (takeAbs(i)) abs(c) else c
            if ( ( ca >= maxV ) || ( ca < -maxV ) ) {
              p.println(f"//    !!!RANGE ERROR!!!");
            }
            val cmask = ca & mask
            p.println(f"      $nAdr'h$j%h: c = $w'h$cmask%h;")
          }
          p.println(f"      default: c = $w'bx;")
          p.println(f"    endcase;");
          p.println(f"  end");
          p.println(f"endmodule // test_$i\n")
        }
      }
    } finally {
      p.close
    }
  }

  // C Generation
  def COut( filename : String, name : String, cmdline : String = "") = {
    val p = new PrintWriter(filename)
    val nAdr = log2Up(interval.size);

    try {
      p.println(f"// $name");
      p.println(f"// $cmdline");
      p.println(f"#define NORDER_$name $nOrder")
      p.println(f"#define ADRMAX_$name $nAdr")
      p.println(f"#define NAME_$name ${'"'}$name${'"'}")
      p.println(f"int ${name}_sign[${nOrder+1}] = "+sign.mkString("{", ", ", "};"))
      p.println(f"int ${name}_cbit[${nOrder+1}] = "+cbit.mkString("{", ", ", "};"))
      p.println(f"int ${name}_calcWidth[${nOrder+1}] = "+calcWidth.mkString("{", ", ", "};"))
      p.println(f"int64_t ${name}_table[$nAdr][${nOrder+1}] = {")
      for ( iv <- interval ) {
        p.print(f"  "+iv.coeff.mkString("{", ", ", "},"))
        p.println(f" // [${iv.xMin}, ${iv.xMin+iv.w}) d=${iv.w}")
      }
      p.println("};\n")
    } finally {
      p.close
    }
  }
}

class FuncTableInt (t: FuncTableDouble, val bp: Int,
    // if provided, use this as calcWidth and cbit, respectively.
    calcWidthSetting: Option[Seq[Int]] = None,
    cbitSetting: Option[Seq[Int]] = None
  ) extends FuncTable {

//   println(s"bp        : $bp")

  val nOrder = t.nOrder
  val sign = t.checkSign
  // Coefficient max
  val minMax = t.minMax
  val absMax = minMax.map( x => max( abs(x._1), abs(x._2) ) )
//   println("minMax    : "+minMax.mkString("", ", ",""))
//   println("absMax    : "+absMax.mkString("", ", ",""))

  val cbit   = cbitSetting match {
    case Some(cb) => cb
    case None     => absMax.map {
      // Here, we always include sign bit, so plus 1.
      x => round(scalb(x, bp)).toLong.toBinaryString.length()+1
    }
  }
//   println("Cbits     : "+cbit.mkString("", ", ",""))

  // Calculate max
  val minMaxCalc = t.getMinMaxAll
  val absMaxCalc = minMaxCalc.map( x => max( abs(x._1), abs(x._2) ) )
//   println("minMaxCalc: "+minMaxCalc.mkString("", ", ",""))
//   println("absMaxCalc: "+absMaxCalc.mkString("", ", ",""))

  val calcWidth = calcWidthSetting match {
    case Some(cW) => cW
    case None     => absMaxCalc.map {
      x => round(scalb(x, bp)).toLong.toBinaryString.length()+1
    }
  }
//   println("CalcWidth : "+calcWidth.mkString("", ", ",""))
//   println(f"ResRange  : ${minMaxCalc.head._1} ~ ${minMaxCalc.head._2}")

  val interval = t.interval.map( x => new FuncTableIntervalInt( x, false, calcWidth, bp ) )

  //interval.foreach(iv => println(iv.xMin, iv.w))

  def checkWidthAll( dxBp : Int ) = {
    (t.nOrder to 1 by -1).toSeq.map(
      n => interval.foldLeft(0)( (w, iv) => max(w, iv.checkWidth(dxBp, n)) )
    )
  }

  def eval(x : Double, dxBp : Int) = {
    interval.find( t => t.inRange(x) ) match {
      case Some(iv) => iv.evalD(x, dxBp, true)
      case None     => { println(f"${this.getClass.getName}.eval: Range Error"); 0.0 }
    }
  }
}

// Currently, always add sign bits
// cbitGiven must include sign bits
class FuncTableFloat (t: FuncTableDouble, val getWidth: Boolean, val cbitGiven : Seq[Int] ) extends FuncTable {
  val nOrder = t.nOrder
  val sign = t.checkSign
  println("Cbits     : "+cbitGiven.mkString("", ", ",""))
  val cbit = if (getWidth) {
    // Get max after calc
    val relScaleMax = (1 to nOrder).toSeq.map(n => t.interval.map(x => x.getRelativeScale(n)).max)
    cbitGiven(0) +: relScaleMax.map( n => n + cbitGiven(0) )
  } else {
    cbitGiven
  }
  println("Cbits     : "+cbit.mkString("", ", ",""))

  // Calculate max
  val minMaxCalc = t.getMinMaxAll
  val absMaxCalc = minMaxCalc.map( x => max( abs(x._1), abs(x._2) ) )
  println("CalcMax   : "+absMaxCalc.mkString("", ", ",""))
  val calcWidth = cbit
  //val calcWidth = absMaxCalc.map(x => bp + java.lang.Math.getExponent(x)+1+1)
  //println("CalcWidth : "+calcWidth.mkString("", ", ",""))
  //println(f"ResRange  : ${minMaxCalc.head._1} ~ ${minMaxCalc.head._2}")

  val interval = t.interval.map( x => new FuncTableIntervalInt(x, true, cbit) )

  //interval.foreach(iv => println(iv.xMin, iv.w))
  /*
  def checkWidthAll( dxBp : Int ) = {
    (t.nOrder to 1 by -1).toSeq.map(
      n => interval.foldLeft(0)( (w, iv) => max(w, iv.checkWidth(dxBp, n)) )
    )
  }
   */

  def eval(x : Double, dxBp : Int) = {
    interval.find( t => t.inRange(x) ) match {
      case Some(iv) => iv.evalD(x, dxBp, true)
      case None     => { println(f"${this.getClass.getName}.eval: Range Error"); 0.0 }
    }
  }

}

/*


void
mdg_table_integer_conversion_coulomb
(
 int one,     // value for "1"
 mdg_func_table *t
 ) {
  t->nprec=one;
  t->floating=0;
  //  int cbit_total=0;
  check_sign(t);
  for(int i=0;i<=t->norder;++i) {
    double cmin=t->cp[i];
    double cmax=t->cp[i];
    double cabsmax;
    for(int j=1;j<t->nentries;++j) {
      int k=j*(t->norder+1);
      if (t->cp[k+i]>cmax) cmax=t->cp[k+i];
      if (t->cp[k+i]<cmin) cmin=t->cp[k+i];
    }
    cabsmax = (fabs(cmax)>fabs(cmin)) ? fabs(cmax) : fabs(cmin);

    // ilogb returns 0 when [1,2)
    // c in [1/2,1)   -> Nprecision - Ntable bit
    // c in [1,2)     -> Nprecision - Ntable + 1bit
    // c in [1/4,1/2) -> Nprecision - Ntable - 1bit
    // Here, we always include sign bit, so plus 1.
    t->cbit[i] = ilogb(cabsmax*one) + 1 + 1;  

    if(MDG_DEBUG_FLAG)fprintf(stderr, "Min: %le Max: %le Absmax : %le Nbit : %d sign: %d\n",
	    cmin, cmax, cabsmax, t->cbit[i], t->sign[i]);
  }

  double *cp = t->cp;
  int64_t *cip = t->ci;
  for(int i=0;i<t->nentries;++i) {
    for(int j=0;j<=t->norder;++j) {
      *cip = lround(*cp * one );
      ++cip;
      ++cp;
    }
  }
}

void
mdg_table_integer_conversion_float
(
 int nprec,    // precision for fractions
 mdg_func_table *t
 ) {
  t->nprec=nprec;
  t->floating=1;
  //  int cbit_total=0;
  int n=t->norder+1;
  // check whether sign of each coefficient changes 
  int *sign = t->sign;
  for(int i=0;i<n;++i) sign[i]=2;
  for(int j=0;j<t->nentries;++j) {
    for(int i=0;i<n;++i) {
      if (sign[i]==2) {
	sign[i] =
	  ((t->cp[j*n+i]) > 0) ? 1 :
	  ( ((t->cp[j*n+i]) < 0) ? -1 : 0);
      } if ( (sign[i]==1)||(sign[i]==-1) ) {
	if (sign[i]*t->cp[j*n+i] < 0.0) {
	  sign[i] = 0;
	}
      }
    }
  }
  for(int i=0; i<n; ++i) t->cbit[i] = 0;

  double *cp = t->cp;
  int64_t * cip = t->ci;
  t->cbit[0]=nprec+2; // with sign
  for(int j=0;j<t->nentries;++j) {
    // In some case c0<cn, so search the max coefficient
    double cmax=0.0;
    int imax=-1;
    for(int i=0;i<=t->norder;++i) {
      if (fabs(cp[i])>cmax) {
	cmax = fabs(cp[i]);
	imax = i;
      }
    }
    if(MDG_DEBUG_FLAG)fprintf(stderr, "cmax %f at %d\n",cmax,imax);
    
    // note : if c0==-2^k, it will fit in k+1 bits, but here k+2 bits will be used.
    int c0exp = ilogb(cmax); // 2^e <= c0 < 2^(e+1)
    // scale by 2^(nprec-1-e) (without sign case)
    int shift = nprec-c0exp;
    t->exp[j] = -shift;
    
    for(int i=0;i<=t->norder;++i) {
      double cn = *cp;
      cn = scalbn(cn, shift);
      int cnbits = ilogb(fabs(cn));
      if (cn==-(1<<cnbits)) cnbits +=1;
      else cnbits +=2;
      if (cnbits > t->cbit[i]) t->cbit[i]=cnbits;
      *cip = llround(cn);
      ++cp;
      ++cip;
    }
  }

  for(int i=0;i<n;++i) {
    if(MDG_DEBUG_FLAG)fprintf(stderr, "order %d: Nbit %d sign: %d\n",
	    i, t->cbit[i], t->sign[i]);
  }

  t->hsign[0] = t->sign[0];
  for(int i=1;i<=t->norder;++i) {
    t->hsign[i] = t->sign[i] * t->sign[i-1];
  }

  // Process exponents
  int expmax = t->exp[0];
  int expmin = expmax;
  for(int j=1;j<t->nentries;++j) {
    if (t->exp[j] > expmax) expmax = t->exp[j];
    if (t->exp[j] < expmin) expmin = t->exp[j];
  }
  int dexpn = expmax-expmin+1; // include exponent for 0
  int expbits = 0;
  while (dexpn) {
    ++expbits;
    dexpn>>=1;
  }
  t->expbits = expbits;
  for(int j=0;j<t->nentries;++j)
    t->exp[j] += -expmin + 1; // expmin = 1
  t->exporg = expmin-1;

  if(MDG_DEBUG_FLAG)fprintf(stderr, "exponent Nbit %d: Max %d\n", expbits, expmax-expmin);
}


#ifndef GNUPLOT_EXE
#define GNUPLOT_EXE "/usr/local/bin/gnuplot"
//#define GNUPLOT_EXE "/usr/bin/gnuplot"
#endif

void
mdg_table_plot(int n, int ylog, double (*f)(double), mdg_func_table *t)
{
  char *datafilename = (char *) malloc(sizeof(char)*(strlen(t->name)+100));
  char *plotfilename = (char *) malloc(sizeof(char)*(strlen(t->name)+100));
  strcpy(datafilename, t->name);
  strcpy(plotfilename, t->name);
  strcat(datafilename, "_err.dat");
  strcat(plotfilename, "_err.gp");
  FILE *fp = fopen(datafilename, "w");
  double xmin, xmax;
  get_minmax(t, &xmin, &xmax);

  // Index 0: Theoretical Error, absolute and relative
  for(int i=0; i<t->nentries; ++i) {
    double err = fabs(t->c1[(t->norder+2)*(i+1)-1]);
    double errlsb = scalbn(err, t->nprec);
    if (t->floating) {
      errlsb = scalbn(err, -(t->exp[i]+t->exporg));
    }
    double relerr = 0.0;
    double x = t->xmin[i]+t->w[i]*0.5;
    double exact = (*f)(x);
    if (exact!=0.0) relerr = fabs(err/exact);
    fprintf(fp,"%f %le %le %le\n", x, err, errlsb, relerr);
  }
  fprintf(fp,"\n\n");

  // Index 1: Error based on calculations in double and int
  fprintf(fp,"#1:x 2:err(dbl) 3:err(int) 4:err(dbl/lsb) 5:err(int/lsb) 6:y(exact) 7:y(d) 8:y(i) 9:relerr(d) 10:relerr(i)\n");
  for(int i=0; i<t->nentries; ++i) {
    double x = t->xmin[i];
    double w = t->w[i]/n;
    double lsb = scalbn(1.0, t->exp[i]+t->exporg);
    for(int j=0; j<n; ++j) {
      double y_exact  = (*f)(x);
      double y_poly_d = mdg_table_eval(x,i,t);
      double y_poly_i = mdg_table_eval_i(x,18,t);
      double err_d = fabs(y_poly_d-y_exact);
      double err_i = fabs(y_poly_i-y_exact);
      double relerr_d = 0.0;
      double relerr_i = 0.0;
      if (y_exact!=0.0) {
	relerr_d = fabs(err_d/y_exact);
	relerr_i = fabs(err_i/y_exact);
      }
      fprintf(fp,"%f %le %le %le %le %le %le %le %le %le\n",x,
	      err_d, err_i, err_d/lsb, err_i/lsb, y_exact, y_poly_d, y_poly_i, relerr_d, relerr_i);
      x += w;
    }
  }
  fprintf(fp,"e\n");
  fclose(fp);

  //fp = popen(GNUPLOT_EXE, "w");
  //pclose(fp);
  fp = fopen(plotfilename, "w");

  fprintf(fp,"set xrange [%f:%f]\n", xmin, xmax);
  if (ylog) fprintf(fp,"set logscale y\n");
  fprintf(fp,"file='%s'\n", datafilename);
  fprintf(fp,"plot file index 0 using 1:2 with steps lt rgb 'black' t 'Theoretical',");
  //  fprintf(fp,"     file index 1 using 1:2 with lines lt rgb 'blue' ,");
  fprintf(fp,"     file index 1 using 1:3 with lines lt rgb 'red' t 'Actual'\n");
  fprintf(fp,"pause -1\n");
  fprintf(fp,"plot file index 0 using 1:3 with steps lt rgb 'black' t 'Theoretical/lsb',");
  //  fprintf(fp,"     file index 1 using 1:4 with lines lt rgb 'blue' ,");
  fprintf(fp,"     file index 1 using 1:5 with lines lt rgb 'red' t 'Actual/lsb'\n");
  fprintf(fp,"pause -1\n");
  fprintf(fp,"plot file index 0 using 1:4 with steps lt rgb 'black' t 'Theoretical relative',");
  fprintf(fp,"     file index 1 using 1:10 with lines lt rgb 'red' t 'Actual relative'\n");
  fprintf(fp,"pause -1\n");
  if (ylog) fprintf(fp,"set nologscale y\n");
  fprintf(fp,"plot file index 0 using 1:4 with steps lt rgb 'black' t 'Theoretical relative'\n");
  fprintf(fp,"pause -1\n");
  
  fclose(fp);
  
}

 */
