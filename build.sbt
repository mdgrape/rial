// See README.md for license details.

lazy val rial = (project in file("."))
  .settings(
    version      := "0.1.0",
    scalaVersion := "2.12.10",
    libraryDependencies ++= Seq(
      "edu.berkeley.cs" %% "chisel3" % "3.5.6",
      "edu.berkeley.cs" %% "chiseltest" % "0.5.6" % "test",
      "net.java.dev.jna" % "jna" % "5.10.0",
      "net.java.dev.jna" % "jna-platform" % "5.10.0",
      "org.jfree" % "jfreechart" % "1.5.0",
      "org.apache.commons" % "commons-math3" % "3.2",
      "org.typelevel" %% "spire" % "0.17.0",
      "org.scala-lang" % "scala-compiler" % scalaVersion.value,
      ),
    scalacOptions ++= Seq(
      "-language:reflectiveCalls",
      "-deprecation",
      "-feature",
      "-Xcheckinit",
      "-P:chiselplugin:genBundleElements"
      ),
    addCompilerPlugin("edu.berkeley.cs" % "chisel3-plugin" % "3.5.6" cross CrossVersion.full)
  )

Test / parallelExecution := false
