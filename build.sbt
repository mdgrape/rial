// See README.md for license details.

ThisBuild / scalaVersion     := "2.13.18"
val chiselVersion = "7.11.0"

// Settings shared between the library and the synthesis-support subprojects.
lazy val commonSettings = Seq(
  scalacOptions ++= Seq(
    "-language:reflectiveCalls",
    "-deprecation",
    "-feature",
    "-Xcheckinit",
    ),
  addCompilerPlugin("org.chipsalliance" % "chisel-plugin" % chiselVersion cross CrossVersion.full),
)

lazy val rial = (project in file("."))
  .settings(commonSettings)
  .settings(
    version      := "0.1.0",
    libraryDependencies ++= Seq(
      "org.chipsalliance" %% "chisel" % chiselVersion,
      "org.scalatest" %% "scalatest" % "3.2.19" % "test",
      "org.jfree" % "jfreechart" % "1.5.0",

      "org.apache.commons" % "commons-math3" % "3.2",
      "org.typelevel" %% "spire" % "0.17.0",
      "org.scala-lang" % "scala-compiler" % scalaVersion.value,
      "org.hipparchus" % "hipparchus-optim" % "4.0.3",
      "org.hipparchus" % "hipparchus-core"  % "4.0.3"
      ),

    // Run test suites one at a time: verilator builds break intermittently
    // ("attempted to destroy locked Thread Pool") when several simulations
    // are built concurrently.
    Test / parallelExecution := false,

    // Fat jar for the JNI binding (jni/, rialcpp): the C++ side loads a single
    // jar via -Djava.class.path, so all runtime deps must be bundled.
    assembly / assemblyJarName := "rial.jar",
    assembly / assemblyMergeStrategy := {
      case PathList("META-INF", _ @ _*)          => MergeStrategy.discard
      case x if x.endsWith("module-info.class")  => MergeStrategy.discard
      case _                                     => MergeStrategy.first
    },
  )

// SystemVerilog emitter used by the FPGA resource-estimation flow (syn/fpga).
// Kept out of `rial` so that `rial/assembly` and `rial/test` are unaffected.
lazy val synEmit = (project in file("syn/emit"))
  .dependsOn(rial)
  .settings(commonSettings)
  .settings(
    name    := "rial-syn-emit",
    version := "0.1.0",
    libraryDependencies ++= Seq(
      "org.chipsalliance" %% "chisel" % chiselVersion,
      ),
    // Forked runs must see the build root so that relative --out paths
    // (e.g. syn/fpga/rtl) resolve the same way as from the shell.
    Compile / run / fork := true,
    Compile / run / baseDirectory := (LocalRootProject / baseDirectory).value,
    assembly / assemblyJarName := "rial-syn-emit.jar",
  )
