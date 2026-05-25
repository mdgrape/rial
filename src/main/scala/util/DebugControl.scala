//% @file DebugControl.scala
//
// Dynamic Debug Control.
//
// Mix `DebugControl` into a Module to get an optional run-time debug-enable
// input plus a guarded `dbgPrintf`. The enable port (and the printf hardware)
// exist only when `enableDebug` is true, so a non-debug build carries no extra
// logic. Whoever owns the module -- a test, or a parent that threads its own
// enable down -- drives `debugEnableIO` directly. No cross-module boring, so no
// dependency on the (deprecated) BoringUtils name-based API.
//
// Copyright (C) Makoto Taiji RIKEN BDR 2020
//
package rial.util

import chisel3._
import chisel3.util._

trait DebugControl extends Module {
  val enableDebug : Boolean

  val debugEnableIO : Option[Bool] =
    if (enableDebug) Some(IO(Input(Bool())).suggestName("debugEnableIO")) else None

  def dbgPrintf( fmt : String, data : Bits* ) : Unit = {
    if (enableDebug) {
      when (debugEnableIO.get) {
        printf(fmt, data : _*)
      }
    }
  }
}
