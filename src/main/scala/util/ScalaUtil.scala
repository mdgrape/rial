package rial.util

import chisel3._
import chisel3.util._
import scala.math._

import spire.math.SafeLong
import spire.implicits._
import java.nio.ByteBuffer

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

  def maskBI( n : Int ) : BigInt = {
    if (n<=0) { BigInt(0) }
    else (BigInt(1) << n) - 1
  }

  def bit( n : Int, x : Long ) : Long = { ((x>>n) & 1) }
  def slice( n : Int, w : Int, x : Long ) : Long = { (x>>n) & maskL(w) }

  def bit( n : Int, x : Int ) : Int = { (x>>n) & 1 }
  def slice( n : Int, w : Int, x : Int ) : Int = { (x>>n) & maskI(w) }

  def bit( n : Int, x : BigInt ) : Int = { ((x>>n) & 1).toInt }
  def slice( n : Int, w : Int, x : BigInt ) : BigInt = { (x>>n) & maskBI(w) }

  def bit( n : Int, x : SafeLong ) : Int = { ((x>>n) & 1).toInt }
  def slice( n : Int, w : Int, x : SafeLong ) : SafeLong = { (x>>n) & maskSL(w) }

  def normalizeFloat( x : Float ) : Float = {
    val x0 = if (x.isNaN) 0x7FC00000 else java.lang.Float.floatToRawIntBits(x)
    if (slice(23,8,x0)==0) 0.0f else java.lang.Float.intBitsToFloat(x0)
  }

  def castFromDouble( x : Double ) : Long = {
    java.lang.Double.doubleToRawLongBits(x)
  }

  def castToDouble( x : Long ) : Double = {
    java.lang.Double.longBitsToDouble(x)
  }

  def unpackDouble( x : Double ) : (Int, Int, Long) = {
    val lv = castFromDouble(x)
    ( bit(63,lv).toInt, slice(52,11,lv).toInt, slice(0,52,lv) )
  }

  def packDouble( sgn : Int, ex : Int, man : Long ) : Double = {
    val lv : Long = ((sgn&1).toLong<<63)+((ex&0x7FF).toLong<<52)+
    (man & 0xFFFFFFFFFFFFFL)
    castToDouble(lv)
  }

  def castFromFloat( x : Float ) : Int = {
    java.lang.Float.floatToRawIntBits(x)
  }

  def castToFloat( x : Int ) : Float = {
    java.lang.Float.intBitsToFloat(x)
  }

  def unpackFloat( x : Float ) : (Int, Int, Int) = {
    val iv = castFromFloat(x)
    ( bit(31,iv), slice(23,8,iv), slice(0,23,iv) )
  }

  def packFloat( sgn : Int, ex : Int, man : Int ) : Float = {
    val iv : Int = ((sgn&1)<<31)+((ex&0xFF)<<23)+(man & 0x7FFFFF)
    castToFloat(iv)
  }

}
