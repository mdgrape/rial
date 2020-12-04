import spire.implicits._
import spire.math.Integral
import spire.algebra._

// private def mask[T]( pos : Int ) (implicit numeric: Numeric[T]) : T = { ...
//  private def bit[T]( x : T, pos : Int )(implicit ev: T <:< { def >>(x:Int) : T }) : Int = {

object genericsTest extends App {
  def test[T](x: T, y: T)(implicit ev: Integral[T]) : T = { x + y }

  println( test[Int](3, 5) )
  val x = test[Long](0xFFFFFFFFFFFL, 1L) 
  println(f"$x%x")

  /*
  def shift[T](x: T, y: T)(implicit ev: spire.math.Integral[T]) : T = { x << y }

  println( shift[Int](1, 5) )
  val z = shift[Long](1L, 40L) 
  println(f"$z%x")
   */

  /*
  //def test[T](x: T, y: T)(implicit ev: HasAnd[T]) : T = { x & y }
  def test[T](x: T, y: T)(implicit ev: T <:< { def &(x: T) : T }) : T = { x & y }

  println( test[Int](3, 7) )
  println( test[Long](33L, 77L) )
   */
}
