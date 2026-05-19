// See README.md for license details.

lazy val rial = (project in file("."))
  .settings(
    version      := "0.1.0",
    scalaVersion := "2.12.18",
    libraryDependencies ++= Seq(
      "edu.berkeley.cs" %% "chisel3" % "3.6.1",
      "edu.berkeley.cs" %% "chiseltest" % "0.6.1" % "test",
      "net.java.dev.jna" % "jna" % "5.10.0",
      "net.java.dev.jna" % "jna-platform" % "5.10.0",
      "org.jfree" % "jfreechart" % "1.5.0",
      "org.apache.commons" % "commons-math3" % "3.2",
      "org.typelevel" %% "spire" % "0.17.0",
      "org.scala-lang" % "scala-compiler" % scalaVersion.value,
      "org.hipparchus" % "hipparchus-optim" % "4.0.3",
      "org.hipparchus" % "hipparchus-core"  % "4.0.3"
      ),
    scalacOptions ++= Seq(
      "-language:reflectiveCalls",
      "-deprecation",
      "-feature",
      "-Xcheckinit",
      ),
    addCompilerPlugin("edu.berkeley.cs" % "chisel3-plugin" % "3.6.1" cross CrossVersion.full),
    // fix jar file name that is used by rial-jni
    assembly / assemblyJarName := "rial.jar"
  )

Test / parallelExecution := false
