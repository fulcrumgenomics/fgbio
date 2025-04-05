val gitHeadCommitSha = settingKey[String]("current git commit SHA")
gitHeadCommitSha in ThisBuild := scala.sys.process.Process("git rev-parse --short HEAD").lineStream.head

// *** IMPORTANT ***
// One of the two "version" lines below needs to be uncommented.
val _version = "2.5.3"
//ThisBuild / version := _version // the release version
ThisBuild / version := s"${_version}-${gitHeadCommitSha.value}-SNAPSHOT" // the snapshot version
