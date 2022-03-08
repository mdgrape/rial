//% @file CrossMultPackedFPGeneric.scala

package rial.simd

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._
import rial.arith._
import rial.util._

// bonded-engine/emulate/fpu.hpp
//
// ```cpp
// inline void select_source12_vregfx(int vfregindex1,int vfregindex2)
// {
//   vec4f_t v1 = vfreg.read(vfregindex1);
//   vec4f_t v2 = vfreg.read(vfregindex2);
//   sfpu[0].set_source1(v1[1]); sfpu[0].set_source2(v2[2]);
//   sfpu[1].set_source1(v1[2]); sfpu[1].set_source2(v2[0]);
//   sfpu[2].set_source1(v1[0]); sfpu[2].set_source2(v2[1]);
//   sfpu[3].set_source1(v1[3]); sfpu[3].set_source2(v2[3]);
// }
//
// void vfxmulv(int vfd, int vfs1, int vfs2)
// {
//   select_source12_vregfx(vfs1,vfs2);
//   for(int i=0;i<4;i++)sfpu[i].fmul();
//   set_results_to_vregf(vfd);
//   debug_print(std::string("cross*").c_str(),vfd,vfs1,vfs2,true,true,true);
// }
// ```

class CrossMultPackedFPGeneric( // length is fixed: 4.
  xSpec : RealSpec, ySpec : RealSpec, zSpec : RealSpec, // Input / Output floating spec
  roundSpec : RoundSpec, // Rounding spec
  stage : PipelineStageConfig
) extends Module {

  val nStage = stage.total

  def getParam = { (xSpec, ySpec, zSpec, roundSpec, nStage) }

  def getStage = nStage

  val io = IO(new Bundle{
    val x   = Input (Vec(4, UInt(xSpec.W.W)))
    val y   = Input (Vec(4, UInt(ySpec.W.W)))
    val z   = Output(Vec(4, UInt(zSpec.W.W)))
  })

  val mults = (0 until 4).map {i => Module(new MultFPGeneric(
    xSpec, ySpec, zSpec, roundSpec, stage
  ))}

  mults(0).io.x := io.x(1)
  mults(1).io.x := io.x(2)
  mults(2).io.x := io.x(0)
  mults(3).io.x := io.x(3)

  mults(0).io.y := io.y(2)
  mults(1).io.y := io.y(0)
  mults(2).io.y := io.y(1)
  mults(3).io.y := io.y(3)

  io.z(0) := mults(0).io.z
  io.z(1) := mults(1).io.z
  io.z(2) := mults(2).io.z
  io.z(3) := mults(3).io.z
  //printf("x=%x y=%x z=%x\n", io.x, io.y, io.z)
}
