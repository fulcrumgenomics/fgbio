val gitHeadCommitSha = settingKey[String]("current git commit SHA")
gitHeadCommitSha in ThisBuild := scala.sys.process.Process("git rev-parse --short HEAD").lineStream.head

// *** IMPORTANT ***
val _version = "2.5.14"
// One of the two "version" lines below needs to be uncommented.
//ThisBuild / version := _version // the release version
ThisBuild / version := s"${_version}-${gitHeadCommitSha.value}-SNAPSHOT" // the snapshot version
