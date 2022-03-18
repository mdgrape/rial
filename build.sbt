// See README.md for license details.

name := "rial"

version := "0.1.0"

val scala212 = "2.12.10" // chipyard 1.6
val scala213 = "2.13.7"  // chisel 3.5
crossScalaVersions := Seq(scala212, scala213)

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

libraryDependencies += "org.scala-lang" % "scala-compiler" % scalaVersion.value

libraryDependencies ++= {
  CrossVersion.partialVersion(scalaVersion.value) match {
    case Some((major, minor)) => Seq(
      "org.scalactic"  % f"scalactic_${major}.${minor}" % "3.1.+",
      "org.scalatest"  % f"scalatest_${major}.${minor}" % "3.1.+" % "test"
    )
    case _ => Nil
  }
}

// Scala-chart
//libraryDependencies += "com.github.wookietreiber" %% "scala-chart" % "latest.integration"
//libraryDependencies += "com.itextpdf" % "itextpdf" % "5.5.6"
//libraryDependencies += "org.jfree" % "jfreesvg" % "3.0"
// JFreeChart
libraryDependencies += "org.jfree" % "jfreechart" % "1.5.0"

// https://mvnrepository.com/artifact/org.apache.commons/commons-math3
libraryDependencies += "org.apache.commons" % "commons-math3" % "3.2"

// Spire
libraryDependencies += "org.typelevel" %% "spire" % "0.17.0"

//excludeFilter in unmanagedSources := "testThreefry.scala"
//excludeFilter in unmanagedSources ~= { _ || "ThreefryUnitTest.scala"}
//excludeFilter in unmanagedSources ~= { _ || "ThreefryMain.scala"}
