// See README.md for license details.

def scalacOptionsVersion(scalaVersion: String): Seq[String] = {
  Seq() ++ {
    // If we're building with Scala > 2.11, enable the compile option
    //  switch to support our anonymous Bundle definitions:
    //  https://github.com/scala/bug/issues/10047
    CrossVersion.partialVersion(scalaVersion) match {
      case Some((2, scalaMajor: Long)) if scalaMajor < 12 => Seq()
      case _ => Seq("-Xsource:2.11")
    }
  }
}

def javacOptionsVersion(scalaVersion: String): Seq[String] = {
  Seq() ++ {
    // Scala 2.12 requires Java 8. We continue to generate
    //  Java 7 compatible code for Scala 2.11
    //  for compatibility with old clients.
    CrossVersion.partialVersion(scalaVersion) match {
      case Some((2, scalaMajor: Long)) if scalaMajor < 12 =>
        Seq("-source", "1.7", "-target", "1.7")
      case _ =>
        Seq("-source", "1.8", "-target", "1.8")
    }
  }
}

name := "rial"

version := "0.1.0"

scalaVersion := "2.12.10"

crossScalaVersions := Seq("2.12.10", "2.11.12")

resolvers ++= Seq(
  Resolver.sonatypeRepo("snapshots"),
  Resolver.sonatypeRepo("releases")
)

// Provide a managed dependency on X if -DXVersion="" is supplied on the command line.
val defaultVersions = Map(
  "chisel3" -> "3.4-SNAPSHOT",
  "chiseltest" -> "0.3-SNAPSHOT"
  )

libraryDependencies ++= Seq("chisel3","chiseltest").map {
  dep: String => "edu.berkeley.cs" %% dep % sys.props.getOrElse(dep + "Version", defaultVersions(dep)) }

libraryDependencies ++= Seq(
  "net.java.dev.jna" % "jna" % "5.5.0",
  "net.java.dev.jna" % "jna-platform" % "5.5.0")

// Scala-chart
//libraryDependencies += "com.github.wookietreiber" %% "scala-chart" % "latest.integration"
//libraryDependencies += "com.itextpdf" % "itextpdf" % "5.5.6"
//libraryDependencies += "org.jfree" % "jfreesvg" % "3.0"
// JFreeChart
libraryDependencies += "org.jfree" % "jfreechart" % "1.5.0"

// Scala-compiler
libraryDependencies += "org.scala-lang" % "scala-compiler" % "2.12.10"

// ScalaTest
libraryDependencies += "org.scalactic" % "scalactic_2.12" % "3.2.0"
libraryDependencies += "org.scalatest" % "scalatest_2.12" % "3.2.0" % "test"

// https://mvnrepository.com/artifact/org.apache.commons/commons-math3
libraryDependencies += "org.apache.commons" % "commons-math3" % "3.2"

// Spire
libraryDependencies += "org.typelevel" %% "spire" % "0.14.1"

scalacOptions ++= scalacOptionsVersion(scalaVersion.value)
scalacOptions ++= Seq("-feature","-deprecation")

javacOptions ++= javacOptionsVersion(scalaVersion.value)

//excludeFilter in unmanagedSources := "testThreefry.scala"
//excludeFilter in unmanagedSources ~= { _ || "ThreefryUnitTest.scala"}
//excludeFilter in unmanagedSources ~= { _ || "ThreefryMain.scala"}
