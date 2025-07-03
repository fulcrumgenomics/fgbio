import com.typesafe.sbt.SbtGit.GitCommand
import sbt._
import sbtassembly.AssemblyKeys.assembly
import sbtassembly.MergeStrategy
import sbtrelease.ReleasePlugin.autoImport.ReleaseTransformations._
import scoverage.ScoverageKeys._

////////////////////////////////////////////////////////////////////////////////////////////////
// We have the following "settings" in this build.sbt:
// - versioning with sbt-release
// - custom JAR name for the root project
// - settings to publish to Sonatype
// - scaladoc settings
// - custom merge strategy for assembly
////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////
// Use sbt-release to bupm the version numbers.
//
// see: http://blog.byjean.eu/2015/07/10/painless-release-with-sbt.html
////////////////////////////////////////////////////////////////////////////////////////////////

// Release settings
releaseVersionBump := sbtrelease.Version.Bump.Next
releasePublishArtifactsAction := PgpKeys.publishSigned.value
releaseProcess := Seq[ReleaseStep](
  checkSnapshotDependencies,
  inquireVersions,
  runClean,
  runTest,
  setReleaseVersion,
  commitReleaseVersion,
  tagRelease,
  releaseStepCommand("publishSigned"),
  setNextVersion,
  commitNextVersion,
  releaseStepCommand("sonatypeReleaseAll"),
  pushChanges
)

////////////////////////////////////////////////////////////////////////////////////////////////
// For the aggregate (root) jar, override the name.  For the sub-projects,
// see the build.sbt in each project folder.
////////////////////////////////////////////////////////////////////////////////////////////////
assembly / assemblyJarName := "fgbio-" + version.value + ".jar"

////////////////////////////////////////////////////////////////////////////////////////////////
// Sonatype settings
////////////////////////////////////////////////////////////////////////////////////////////////
publishMavenStyle := true
publishTo := {
  val centralSnapshots = "https://central.sonatype.com/repository/maven-snapshots/"
  if (isSnapshot.value) Some("central-snapshots" at centralSnapshots)
  else localStaging.value
}
Test / publishArtifact := false
pomIncludeRepository := { _ => false }
credentials ++= (for {
  username <- Option(System.getenv().get("SONATYPE_USER"))
  password <- Option(System.getenv().get("SONATYPE_PASS"))
} yield Credentials("Sonatype Nexus Repository Manager", "central.sonatype.com", username, password)).toSeq

////////////////////////////////////////////////////////////////////////////////////////////////
// Coverage settings: don't include personal packages in coverage counts
////////////////////////////////////////////////////////////////////////////////////////////////
coverageExcludedPackages := "com.fulcrumgenomics.personal.*;com.fulcrumgenomics.internal.*"
val htmlReportsDirectory: String = "target/test-reports"

////////////////////////////////////////////////////////////////////////////////////////////////
// scaladoc options
////////////////////////////////////////////////////////////////////////////////////////////////
val docScalacOptions = Seq("-groups", "-implicits")

////////////////////////////////////////////////////////////////////////////////////////////////
// Common settings 
////////////////////////////////////////////////////////////////////////////////////////////////

val primaryScalaVersion = "2.13.14"

lazy val commonSettings = Seq(
  organization         := "com.fulcrumgenomics",
  organizationName     := "Fulcrum Genomics LLC",
  organizationHomepage := Some(url("https://fulcrumgenomics.com/")),
  homepage             := Some(url("http://github.com/fulcrumgenomics/fgbio")),
  scmInfo              := Some(ScmInfo(url("http://github.com/fulcrumgenomics/fgbio"), "scm:git@github.com:fulcrumgenomics/fgbio.git")),
  startYear            := Some(2015),
  scalaVersion         := primaryScalaVersion,
  crossScalaVersions   :=  Seq(primaryScalaVersion),
  scalacOptions        ++= Seq(
    "-release:8", // Target the Java 1.8 release
    "-deprecation", // Emit warning and location for usages of deprecated APIs.
    "-explaintypes", // Explain type errors in more detail.
    "-feature", // Emit warning and location for usages of features that should be imported explicitly.
    "-opt-inline-from:com.fulcrumgenomics.**", // Tells the inliner that it is allowed to inline things from these classes.
    "-opt-warnings:at-inline-failed", // Tells you if methods marked with `@inline` cannot be inlined, so you can remove the tag.
    "-opt:inline:com.fulcrumgenomics.**", // Turn on the inliner.
    "-unchecked", // Enable additional warnings where generated code depends on assumptions.
    "-Xcheckinit", // Wrap field accessors to throw an exception on uninitialized access.
    "-Xfatal-warnings", // Fail the compilation if there are any warnings.
    "-Xlint:adapted-args", // Warn if an argument list is modified to match the receiver.
    "-Xlint:constant", // Evaluation of a constant arithmetic expression results in an error.
    "-Xlint:delayedinit-select", // Selecting member of DelayedInit.
    "-Xlint:doc-detached", // A Scaladoc comment appears to be detached from its element.
    "-Xlint:inaccessible", // Warn about inaccessible types in method signatures.
    "-Xlint:infer-any", // Warn when a type argument is inferred to be `Any`.
    "-Xlint:missing-interpolator", // A string literal appears to be missing an interpolator id.
    "-Xlint:nullary-unit", // Warn when nullary methods return Unit.
    "-Xlint:option-implicit", // Option.apply used implicit view.
    "-Xlint:package-object-classes", // Class or object defined in package object.
    "-Xlint:poly-implicit-overload", // Parameterized overloaded implicit methods are not visible as view bounds.
    "-Xlint:private-shadow", // A private field (or class parameter) shadows a superclass field.
    "-Xlint:stars-align", // Pattern sequence wildcard must align with sequence component.
    "-Xlint:type-parameter-shadow", // A local type parameter shadows a type already in scope.
    "-Xlint:unused", // Warn when anything is unused.
    "-Ycache-macro-class-loader:last-modified", // and macro definitions. This can lead to performance improvements.
    "-Ycache-plugin-class-loader:last-modified", // Enables caching of classloaders for compiler plugins
    "-Yopt-inline-heuristics:at-inline-annotated", // Tells the inliner to use your `@inliner` tags.
    "-Yopt-log-inline", "_", // Optional, logs the inliner activity so you know it is doing something.
    "-Ywarn-dead-code", // Warn when dead code is identified.
    "-Ywarn-extra-implicit", // Warn when more than one implicit parameter section is defined.
    "-Ywarn-numeric-widen", // Warn on obsolete octal syntax.
    "-Ywarn-numeric-widen", // Warn when numerics are widened.
    "-Ywarn-unused:implicits", // Warn if an implicit parameter is unused.
    "-Ywarn-unused:imports", // Warn if an import selector is not referenced.
    "-Ywarn-unused:locals", // Warn if a local definition is unused.
    "-Ywarn-unused:params", // Warn if a value parameter is unused.
    "-Ywarn-unused:patvars", // Warn if a variable bound in a pattern is unused.
    "-Ywarn-unused:privates", // Warn if a private member is unused.
    "-Ywarn-value-discard", // Warn when non-Unit expression results are unused.
  ),
  Compile / doc / scalacOptions ++= docScalacOptions,
  Test / doc / scalacOptions    ++= docScalacOptions,
  useCoursier          :=  false,
  Test / testOptions   += Tests.Argument(TestFrameworks.ScalaTest, "-h", Option(System.getenv("TEST_HTML_REPORTS")).getOrElse(htmlReportsDirectory)),
  // uncomment for full stack traces
  // Test / testOptions   += Tests.Argument("-oDF"),
  Test / fork          := true,
  resolvers            += Resolver.sonatypeCentralSnapshots,
  resolvers            += Resolver.sonatypeCentralRepo("releases"),
  resolvers            += Resolver.mavenLocal,
  resolvers            += "broad-snapshots" at "https://broadinstitute.jfrog.io/artifactory/libs-snapshot/",
  shellPrompt          := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), version.value) },
  updateOptions        := updateOptions.value.withCachedResolution(true),
  pomExtra             := <url>https://github.com/fulcrumgenomics/fgbio</url>
    <licenses>
      <license>
        <name>MIT License</name>
        <url>https://www.opensource.org/licenses/mit-license.html</url>
      </license>
    </licenses>
    <developers>
      <developer>
        <id>nh13</id>
        <name>Nils Homer</name>
        <url>https://github.com/nh13</url>
        <email>nils@fulcrumgenomics.com</email>
      </developer>
      <developer>
        <id>tfenne</id>
        <name>Tim Fennell</name>
        <url>https://github.com/tfenne</url>
        <email>tim@fulcrumgenomics.com</email>
      </developer>
      <developer>
        <id>clintval</id>
        <name>Clint Valentine</name>
        <url>https://github.com/clintval</url>
        <email>clint@fulcrumgenomics.com</email>
      </developer>
    </developers>
) ++ Defaults.coreDefaultSettings

////////////////////////////////////////////////////////////////////////////////////////////////
// root project
////////////////////////////////////////////////////////////////////////////////////////////////
lazy val htsjdkExcludes = Seq(
  ExclusionRule(organization="org.apache.ant"),
  ExclusionRule(organization="gov.nih.nlm.ncbi"),
  ExclusionRule(organization="org.testng"),
  ExclusionRule(organization="com.google.cloud.genomics")
)

lazy val assemblySettings = Seq(
  assembly / test     := {},
  assembly / logLevel := Level.Info
)
lazy val root = Project(id="fgbio", base=file("."))
  .settings(commonSettings: _*)
  .settings(assemblySettings: _*)
  .settings(description := "fgbio")
  .settings(mainClass := Some("com.fulcrumgenomics.cmdline.FgBioMain"))  
  .settings(
    libraryDependencies ++= Seq(
      "org.scala-lang"            %  "scala-reflect"  % scalaVersion.value,
      "org.scala-lang"            %  "scala-compiler" % scalaVersion.value,
      "org.scala-lang.modules"    %% "scala-xml"      % "2.1.0",
      "com.fulcrumgenomics"       %% "commons"        % "1.6.0",
      "com.fulcrumgenomics"       %% "sopt"           % "1.1.0",
      "com.github.samtools"       %  "htsjdk"         % "3.0.5" excludeAll(htsjdkExcludes: _*),
      "org.apache.commons"        %  "commons-math3"  % "3.6.1",
      "com.beachape"              %% "enumeratum"     % "1.7.0",
      "com.intel.gkl"             %  "gkl"            % "0.8.10",

      //---------- Test libraries -------------------//
      "org.scalatest"             %% "scalatest"     % "3.1.3"  % "test->*" excludeAll ExclusionRule(organization="org.junit", name="junit")
  ))
  .settings(dependencyOverrides ++= Seq(
      "org.apache.logging.log4j" % "log4j-api"   % "[2.17.0,)",
      "org.apache.logging.log4j" % "log4j-core"  % "[2.17.0,)",
      "org.xerial.snappy"        % "snappy-java" % "[1.1.8.4,)"
  ))


////////////////////////////////////////////////////////////////////////////////////////////////
// Merge strategy for assembly
////////////////////////////////////////////////////////////////////////////////////////////////
val customMergeStrategy: String => MergeStrategy = {
  case x if Assembly.isConfigFile(x) =>
    MergeStrategy.concat
  case PathList(ps@_*) if Assembly.isReadme(ps.last) || Assembly.isLicenseFile(ps.last) =>
    MergeStrategy.rename
  case PathList("META-INF", xs@_*) =>
    xs map {
      _.toLowerCase
    } match {
      case ("manifest.mf" :: Nil) | ("index.list" :: Nil) | ("dependencies" :: Nil) =>
        MergeStrategy.discard
      case ps@(x :: xt) if ps.last.endsWith(".sf") || ps.last.endsWith(".dsa") =>
        MergeStrategy.discard
      case "plexus" :: xt =>
        MergeStrategy.discard
      case "spring.tooling" :: xt =>
        MergeStrategy.discard
      case "com.google.guava" :: xt =>
        MergeStrategy.discard
      case "services" :: xt =>
        MergeStrategy.filterDistinctLines
      case ("spring.schemas" :: Nil) | ("spring.handlers" :: Nil) =>
        MergeStrategy.filterDistinctLines
      case _ => MergeStrategy.deduplicate
    }
  case "asm-license.txt" | "overview.html" =>
    MergeStrategy.discard
  case "logback.xml" =>
    MergeStrategy.first
  case _ => MergeStrategy.deduplicate
}
assembly / assemblyMergeStrategy := customMergeStrategy
