# CONTRIBUTING

# Publishing to Sonatype

## Requirements

### Set up GPG and your Sonatype Credentials

Set up an access token at https://central.sonatype.org/publish/generate-token/ to use with the username and password below.

The simplest way to setup your Sonatype credentials before doing a release is just:

```console
$ export SONATYPE_USER=...
$ export SONATYPE_PASS=xxx
```

You may need to run the following before using GPG:

```console
export GPG_TTY=$(tty)
```

## Perform the release

1. Use the _release_ version in `version.sbt` and perform a commit.

2. Make a commit and tag, and push the commit and tag triggering the release workflow:


```console
$ git commit -a -m "chore: release X.Y.Z
$ git tag X.Y.Z"
$ git push -f origin main && git push origin X.Y.Z`
```

3. Bump the _patch_ version number and switch to the `-SNAPSHOT` version in `version.sbt`:

```scala
ThisBuild / version := s"X.Y.{Z+1}-${gitHeadCommitSha.value}-SNAPSHOT"
```

4. Make and push a commit:

```console
$ git commit -a -m "chore: bump to X.Y.{Z+1}-SNAPSHOT"
$ git push -f origin main
```

This will kick off the release and docs workflows, and automatically create a release on GitHub.

5. Review release and add context to [the release notes](https://github.com/fulcrumgenomics/fgbio/releases).


6. Update Bioconda

Bioconda **should** auto-detect version changes.
