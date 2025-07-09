



# CONTRIBUTING

# Publishing to Sonatype

## Perform the release

1. Use the _release_ version in `version.sbt` and perform a commit.

2. Make a commit and tag, and push the commit and tag triggering the release workflow:


```console
$ git commit -a -m "chore: release X.Y.Z
$ git tag X.Y.Z"
$ git push -f origin main && git push origin X.Y.Z`
```

3. Bump the version number (in the `_version` value) and switch to the `-SNAPSHOT` version in `version.sbt`:

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
