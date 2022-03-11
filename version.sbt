val gitHeadCommitSha = settingKey[String]("current git commit SHA")
gitHeadCommitSha in ThisBuild := scala.sys.process.Process("git rev-parse --short HEAD").lineStream.head

// *** IMPORTANT ***
// One of the two "version" lines below needs to be uncommented.
version in ThisBuild := "2.0.0-beta1" // the release version
// version in ThisBuild := s"2.0.0-${gitHeadCommitSha.value}-SNAPSHOT" // the snapshot version
