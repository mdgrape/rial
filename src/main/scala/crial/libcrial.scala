package crial

import com.sun.jna.Library
import com.sun.jna.Native
import com.sun.jna.NativeLibrary
import com.sun.jna.Pointer
import java.nio._

/** C library interface for testing.
 */
object crial {

  // Resolve the platform-specific name (libcrial.so / libcrial.dylib / crial.dll)
  // from ./lib, instead of hardcoding the ".so" suffix.
  NativeLibrary.addSearchPath("crial", "./lib")
  val lib = NativeLibrary.getInstance("crial")

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
