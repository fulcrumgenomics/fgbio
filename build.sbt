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
assemblyJarName in assembly := "fgbio-" + version.value + ".jar"

////////////////////////////////////////////////////////////////////////////////////////////////
// Sonatype settings
////////////////////////////////////////////////////////////////////////////////////////////////
publishMavenStyle := true
publishTo := {
  val nexus = "https://oss.sonatype.org/"
  if (isSnapshot.value)
    Some("snapshots" at nexus + "content/repositories/snapshots")
  else
    Some("releases"  at nexus + "service/local/staging/deploy/maven2")
}
publishArtifact in Test := false
pomIncludeRepository := { _ => false }
// For Travis CI - see http://www.cakesolutions.net/teamblogs/publishing-artefacts-to-oss-sonatype-nexus-using-sbt-and-travis-ci
credentials ++= (for {
  username <- Option(System.getenv().get("SONATYPE_USER"))
  password <- Option(System.getenv().get("SONATYPE_PASS"))
} yield Credentials("Sonatype Nexus Repository Manager", "oss.sonatype.org", username, password)).toSeq

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

lazy val commonSettings = Seq(
  organization         := "com.fulcrumgenomics",
  organizationName     := "Fulcrum Genomics LLC",
  homepage             := Some(url("http://github.com/fulcrumgenomics/fgbio")),
  startYear            := Some(2015),
  scalaVersion         := "2.13.0",
  crossScalaVersions   :=  Seq("2.12.8", "2.13.0"),
  scalacOptions        ++= Seq("-target:jvm-1.8", "-deprecation", "-unchecked"),
  scalacOptions in (Compile, doc) ++= docScalacOptions,
  scalacOptions in (Test, doc) ++= docScalacOptions,
  autoAPIMappings := true,
  testOptions in Test  += Tests.Argument(TestFrameworks.ScalaTest, "-h", Option(System.getenv("TEST_HTML_REPORTS")).getOrElse(htmlReportsDirectory)),
  // uncomment for full stack traces
  //testOptions in Test  += Tests.Argument("-oDF"),
  fork in Test         := true,
  resolvers            += Resolver.sonatypeRepo("public"),
  resolvers            += Resolver.mavenLocal,
  resolvers            += "broad-snapshots" at "https://artifactory.broadinstitute.org/artifactory/libs-snapshot/",
  shellPrompt          := { state => "%s| %s> ".format(GitCommand.prompt.apply(state), version.value) },
  updateOptions        := updateOptions.value.withCachedResolution(true)
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
  test in assembly     := {},
  logLevel in assembly := Level.Info
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
      "org.scala-lang.modules"    %% "scala-xml"      % "1.2.0",
      "org.scala-lang.modules"    %% "scala-collection-compat" % "2.1.1",
      "com.fulcrumgenomics"       %% "commons"        % "1.1.0-43c062d-SNAPSHOT",
      "com.fulcrumgenomics"       %% "sopt"           % "1.0.0",
      "com.github.samtools"       %  "htsjdk"         % "2.20.3" excludeAll(htsjdkExcludes: _*),
      "org.apache.commons"        %  "commons-math3"  % "3.6.1",
      "com.beachape"              %% "enumeratum"     % "1.5.13",
      "com.intel.gkl"             %  "gkl"            % "0.8.6",

      //---------- Test libraries -------------------//
      "org.scalatest"             %% "scalatest"     % "3.0.8"  % "test->*" excludeAll ExclusionRule(organization="org.junit", name="junit")
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
assemblyMergeStrategy in assembly := customMergeStrategy
