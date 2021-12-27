val gitHeadCommitSha = settingKey[String]("current git commit SHA")
gitHeadCommitSha in ThisBuild := scala.sys.process.Process("git rev-parse --short HEAD").lineStream.head

// *** IMPORTANT ***
// One of the two "version" lines below needs to be uncommented.
//version in ThisBuild := "1.5.0" // the release version
version in ThisBuild := s"1.5.0-${gitHeadCommitSha.value}-SNAPSHOT" // the snapshot version
