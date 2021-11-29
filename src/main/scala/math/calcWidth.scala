package rial.math

import scala.math._
import rial.math._
import rial.arith._

//
// cbit : acosTable       : Vector(27, 21, 20)
// cbit : atan2Table      : Vector(26, 16, 5)
// cbit : cosTable        : Vector(25, 16, 6)
// cbit : sinTable        : Vector(25, 16, 6)
// cbit : logTable        : Vector(26, 18, 8)
// cbit : expTable        : Vector(26, 18, 7)
// cbit : sqrtTable       : Vector(26, 17, 6)
// cbit : invsqrtTable    : Vector(26, 17, 8)
// cbit : reciprocalTable : Vector(26, 18, 9)
// cbit : maximum         : Vector(27, 21, 20)
// calcWidth : acosTable       : Vector(27, 22, 20)
// calcWidth : atan2Table      : Vector(26, 17, 5)
// calcWidth : cosTable        : Vector(26, 16, 6)
// calcWidth : sinTable        : Vector(26, 16, 6)
// calcWidth : logTable        : Vector(27, 18, 8)
// calcWidth : expTable        : Vector(27, 18, 7)
// calcWidth : sqrtTable       : Vector(27, 17, 6)
// calcWidth : invsqrtTable    : Vector(27, 17, 8)
// calcWidth : reciprocalTable : Vector(27, 18, 9)
// calcWidth : maximum         : Vector(27, 22, 20)
//

object Float32MathTables {
  def main(args: Array[String]) {

    val cbits = Vector(
      ("acosTable      ", ACosSim .acosF32TableI         (0).cbit),
      ("atan2Table     ", ATan2Sim.atan2F32ATanTableI    (0).cbit), // atan2 uses both atan and reciprocal, but reciprocalTable is shared
      ("cosTable       ", CosPiSim.cosPiF32TableI        (0).cbit),
      ("sinTable       ", SinPiSim.sinPiF32TableI        (0).cbit),
      ("logTable       ", new Log2Sim( RealSpec.Float32Spec, 2, 8, 2 ).t.cbit),
      ("expTable       ", ExponentialSim.pow2F32TableI      .cbit), // exp(x) is calculated from pow2
      ("sqrtTable      ", SqrtSim.sqrtF32TableI             .cbit),
      ("invsqrtTable   ", InvSqrtSim.invsqrtF32TableI       .cbit),
      ("reciprocalTable", ReciprocalSim.reciprocalF32TableI .cbit),
      )

    val calcWidths = Vector(
      ("acosTable      ", ACosSim .acosF32TableI         (0).calcWidth),
      ("atan2Table     ", ATan2Sim.atan2F32ATanTableI    (0).calcWidth), // atan2 uses both atan and reciprocal, but reciprocalTable is shared
      ("cosTable       ", CosPiSim.cosPiF32TableI        (0).calcWidth),
      ("sinTable       ", SinPiSim.sinPiF32TableI        (0).calcWidth),
      ("logTable       ", new Log2Sim( RealSpec.Float32Spec, 2, 8, 2 ).t.calcWidth),
      ("expTable       ", ExponentialSim.pow2F32TableI      .calcWidth), // exp(x) is calculated from pow2
      ("sqrtTable      ", SqrtSim.sqrtF32TableI             .calcWidth),
      ("invsqrtTable   ", InvSqrtSim.invsqrtF32TableI       .calcWidth),
      ("reciprocalTable", ReciprocalSim.reciprocalF32TableI .calcWidth),
      )

    for((n, w) <- cbits) {
      println(f"cbit : ${n} : ${w}")
    }
    val maxCbit = cbits.map(t => t._2).reduce( (lhs, rhs) => {
      lhs.zip(rhs).map(x => max(x._1, x._2))
    })
    println(f"cbit : maximum         : ${maxCbit}")

    for((n, w) <- calcWidths) {
      println(f"calcWidth : ${n} : ${w}")
    }
    val maxWidth = calcWidths.map(t => t._2).reduce( (lhs, rhs) => {
      lhs.zip(rhs).map(x => max(x._1, x._2))
    })
    println(f"calcWidth : maximum         : ${maxWidth}")
  }
}
