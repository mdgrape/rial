package rial.util

import chisel3._
import chisel3.util._
import scala.math._

import spire.math.SafeLong
import spire.implicits._

object ScalaUtil {

  def separateOptions(args: Array[String]) : (Array[String], Array[String]) = {
    val pos = args.indexWhere(_=="--")
    val arg0 = if (pos<0) args else args.take(pos)
    val nullArray : Array[String] = Array.empty
    val arg1 = if (pos<0) nullArray else args.drop(pos+1)
    (arg0, arg1)
  }

  def toBinStringFill( x: Long, w: Int ) = {
    var xx = x
    var ww = w
    var s  = ""
    while ( (xx!=0) || (ww>0) ) {
      s = f"${xx%2}" + s
      xx >>=1
      ww -= 1
    }
    s
  }

  def log2UpD( x: Double ) = {
    if (x % 1.0 == 0.0) {
      log2Up(x.toInt)
    } else {
      ceil(log(x)/log(2.0)).toInt
    }
  }

  def log2DownD( x: Double ) = {
    if (x % 1.0 == 0.0) {
      log2Down(x.toInt)
    } else {
      floor(log(x)/log(2.0)).toInt
    }
  }

  def maskL( n : Int ) = {
    if (n<=0) { 0 }
    else if (n<63) { (1L << n)-1 }
    else if (n==63) 0x7FFFFFFFFFFFFFFFL
    else 0xFFFFFFFFFFFFFFFFL
  }

  def maskI( n : Int ) : Int = {
    if (n<=0) { 0 }
    else if (n<31) { (1 << n)-1 }
    else if (n==31) 0x7FFFFFFF
    else 0xFFFFFFFF
  }

  def maskSL( n : Int ) : SafeLong = {
    if (n<=0) { SafeLong(0) }
    else (SafeLong(1) << n) - 1
  }

  def bit( n : Int, x : Long ) : Long = { ((x>>n) & 1) }
  def slice( n : Int, w : Int, x : Long ) : Long = { (x>>n) & maskL(w) }
  def bit( n : Int, x : Int ) : Int = { (x>>n) & 1 }
  def slice( n : Int, w : Int, x : Int ) : Int = { (x>>n) & maskI(w) }

  def bit( n : Int, x : SafeLong ) : Int = { ((x>>n) & 1).toInt }
  def slice( n : Int, w : Int, x : SafeLong ) : SafeLong = { (x>>n) & maskSL(w) }

  def normalizeFloat( x : Float ) : Float = {
    val x0 = if (x.isNaN) 0x7FC00000 else java.lang.Float.floatToRawIntBits(x)
    if (slice(23,8,x0)==0) 0.0f else java.lang.Float.intBitsToFloat(x0)
  }
}
