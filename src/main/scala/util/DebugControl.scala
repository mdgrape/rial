//% @file DynamicDebug.scala
//
// Dynamic Debug Control
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.util

import scala.language.reflectiveCalls
import chisel3._
import chisel3.util._
import chisel3.util._
import chisel3.util.experimental.BoringUtils

trait DebugControlMaster extends Module {

  val debugEnableIO = IO(Input(Bool()))

  BoringUtils.addSource( debugEnableIO, "GlobalDebugEnable" )
  
}

trait DebugControlSlave {
  val enableDebug : Boolean

  val globalDynamicDebugEnable = Wire(Bool())
  globalDynamicDebugEnable := false.B
  if (enableDebug)
    BoringUtils.addSink( globalDynamicDebugEnable, "GlobalDebugEnable" )

  def dbgPrintf( fmt : String, data : Bits* ) = {
    if (enableDebug) {
      when (globalDynamicDebugEnable) {
        printf(fmt, data : _*)
      }
    }
  }

}
