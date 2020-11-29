package crial

import com.sun.jna.Library
import com.sun.jna.Native
import com.sun.jna.NativeLibrary
import com.sun.jna.Pointer
import java.nio._

object crial {

  val lib=NativeLibrary.getInstance("./crial/libcrial.so")

  val mdg_inv_sqrt=lib.getFunction("mdg_inv_sqrt")
  def invSqrt(x:Int) : Int = {
    mdg_inv_sqrt.invokeInt(Array(x.asInstanceOf[java.lang.Integer]))
  }

  val mdg_inv_sqrt3=lib.getFunction("mdg_inv_sqrt3")
  def invSqrt3(x:Int) : Int = {
    mdg_inv_sqrt3.invokeInt(Array(x.asInstanceOf[java.lang.Integer]))
  }

  val mdg_sub32=lib.getFunction("mdg_sub32");
  def sub32(x:Int, y:Int) : (Int, Int, Boolean) = {
    val z  = IntBuffer.allocate(1)
    val zs = IntBuffer.allocate(1)
    val inrange = IntBuffer.allocate(1)

    mdg_sub32.invokeVoid(
      Array(
        x.asInstanceOf[java.lang.Integer],
        y.asInstanceOf[java.lang.Integer],
        z, zs, inrange
      ));
    (z.get, zs.get, inrange.get!=0)
  }

  val mdg_rsquare=lib.getFunction("mdg_rsquare");
  def rsquare(rij: Array[Int], rij_valid:Boolean) : (Int, Int, Boolean, Long, Boolean) = {
    val rij_valid_i = if (rij_valid) 1 else 0
    val r2  = IntBuffer.allocate(1)
    val r2p = IntBuffer.allocate(1)
    val r2_valid = IntBuffer.allocate(1)
    val r2_soft = LongBuffer.allocate(1)
    val r2_soft_gt1 = IntBuffer.allocate(1)

    mdg_rsquare.invokeVoid(
      Array(rij, rij_valid_i.asInstanceOf[java.lang.Integer],
        r2, r2p, r2_valid, r2_soft, r2_soft_gt1));
    (r2.get, r2p.get, r2_valid.get!=0, r2_soft.get, r2_soft_gt1.get!=0)
  }

  val mdg_threefry4x32_init=lib.getFunction("mdg_threefry4x32_init");
  def threefry4x32_init(ctr: Array[Int], key: Array[Int]) : Pointer = {
    mdg_threefry4x32_init.invokePointer( Array(ctr, key) )
  }

  val mdg_threefry4x32=lib.getFunction("mdg_threefry4x32");
  def threefry4x32( p: Pointer ) : Array[Int] = {
    val z  = IntBuffer.allocate(4)
    mdg_threefry4x32.invokeVoid( Array(p, z) )
    z.array()
  }

}
