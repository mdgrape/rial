package rial.util

import chisel3._
import chisel3.util._
import scala.math._

class RealSpec (
  exW : Int,
  manW : Int,
  disableSign : Bool,
  disableNaN  : Bool,
  enableSubnormal : Bool )

//
// class Real
//   this uses BigInt for mantissa, so the speed may be slow
//
class Real ( sgn : Bool, ex : Int, man : BigInt, spec : RealSpec ) {

  def multiply( r : Real, resSpec : RealSpec ) : Real = {
    r
  }

}

