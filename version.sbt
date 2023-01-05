val gitHeadCommitSha = settingKey[String]("current git commit SHA")
gitHeadCommitSha in ThisBuild := scala.sys.process.Process("git rev-parse --short HEAD").lineStream.head

// *** IMPORTANT ***
// One of the two "version" lines below needs to be uncommented.
ThisBuild / version := "2.1.0" // the release version
// ThisBuild / version := s"2.1.1-${gitHeadCommitSha.value}-SNAPSHOT" // the snapshot version
