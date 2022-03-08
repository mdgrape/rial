// See README.md for license details.

name := "rial"

version := "0.1.0"

scalaVersion := "2.13.7"

resolvers ++= Seq(
  Resolver.sonatypeRepo("snapshots"),
  Resolver.sonatypeRepo("releases")
)

// from chisel 3.5, compiler plugin is required
// https://github.com/chipsalliance/chisel3#build-your-own-chisel-projects
addCompilerPlugin("edu.berkeley.cs" % "chisel3-plugin" % "3.5.1" cross CrossVersion.full)

scalacOptions ++= Seq("-language:reflectiveCalls",
                      "-deprecation",
                      "-feature",
                      "-Xcheckinit",
                      "-P:chiselplugin:genBundleElements")

libraryDependencies += "edu.berkeley.cs" %% "chisel3" % "3.5.1"
libraryDependencies += "edu.berkeley.cs" %% "chiseltest" % "0.5.1" % "test"

libraryDependencies ++= Seq(
  "net.java.dev.jna" % "jna" % "5.10.0",
  "net.java.dev.jna" % "jna-platform" % "5.10.0")

// Scala-chart
//libraryDependencies += "com.github.wookietreiber" %% "scala-chart" % "latest.integration"
//libraryDependencies += "com.itextpdf" % "itextpdf" % "5.5.6"
//libraryDependencies += "org.jfree" % "jfreesvg" % "3.0"
// JFreeChart
libraryDependencies += "org.jfree" % "jfreechart" % "1.5.0"

// Scala-compiler
libraryDependencies += "org.scala-lang" % "scala-compiler" % "2.13.7"

// ScalaTest
libraryDependencies += "org.scalactic" % "scalactic_2.13" % "3.1.4"
libraryDependencies += "org.scalatest" % "scalatest_2.13" % "3.1.4" % "test"

// https://mvnrepository.com/artifact/org.apache.commons/commons-math3
libraryDependencies += "org.apache.commons" % "commons-math3" % "3.2"

// Spire
libraryDependencies += "org.typelevel" %% "spire" % "0.17.0"

//excludeFilter in unmanagedSources := "testThreefry.scala"
//excludeFilter in unmanagedSources ~= { _ || "ThreefryUnitTest.scala"}
//excludeFilter in unmanagedSources ~= { _ || "ThreefryMain.scala"}
