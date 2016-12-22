[![Build Status](https://travis-ci.org/fulcrumgenomics/fgbio.svg?branch=master)](https://travis-ci.org/fulcrumgenomics/fgbio)
[![codecov](https://codecov.io/gh/fulcrumgenomics/fgbio/branch/master/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/fgbio)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/fc4f5fe8dbe34bf784114435b202fab4)](https://www.codacy.com/app/contact_32/fgbio?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=fulcrumgenomics/fgbio&amp;utm_campaign=Badge_Grade)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/com.fulcrumgenomics/fgbio_2.11/badge.svg)](https://maven-badges.herokuapp.com/maven-central/com.fulcrumgenomics/fgbio_2.11)
[![Dependency Status](https://www.versioneye.com/user/projects/57a1584c3d8eb6002dc1e812/badge.svg)](https://www.versioneye.com/user/projects/57a1584c3d8eb6002dc1e812#dialog_dependency_badge) 
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/fgbio/blob/master/LICENSE)
[![Language](http://img.shields.io/badge/language-scala-brightgreen.svg)](http://www.scala-lang.org/)

fgbio
====

A set of tools to analyze genomic data.

<!---toc start-->
  * [Goals](#goals)
  * [Building](#building)
  * [Command line](#command-line)
  * [Include fgbio in your project](#include-fgbio-in-your-project)
  * [Overview](#overview)
  * [List of tools](#list-of-tools)
  * [Contributing](#contributing)
  * [Authors](#authors)
  * [License](#license)

<!---toc end-->


# Goals

There are many toolkits available for analyzing genomic data; fgbio does not aim to be all things to all people but is specifically focused on providing:

* Robust, well-tested tools.
* An easy to use command-line.
* Documentation for each tool.
* Tools not found anywhere else.
* Open source development for the benefit of the community and our clients.

## Building 
### Cloning the Repository

[Git LFS](https://git-lfs.github.com/) is used to store large files used in testing fgbio.  In order to compile and run tests it is necessary to [install git lfs](https://git-lfs.github.com/).  To retrieve the large files either:

1. Clone the repository _after_ installing git lfs, or
2. In a previously cloned repository run `git lfs pull` once

After initial setup regular git commands (e.g. `pull`, `fetch`, `push`) will also operate on large files and no special handling is needed.

To clone the repository: `git clone https://github.com/fulcrumgenomics/fgbio.git`

### Running the build
fgbio is built using [sbt](http://www.scala-sbt.org/).

Use ```sbt assembly``` to build an executable jar in ```target/scala-2.11/```.  
Tests may be run with ```sbt test```.
Java SE 8 is required.


## Command line

`java -jar target/scala-2.11/fgbio-0.1.3-SNAPSHOT.jar` to see the commands supported.  Use `java -jar target/scala-2.11/fgbio-0.1.3-SNAPSHOT.jar <command>` to see the help message for a particular command.

## Include fgbio in your project

You can include `fgbio` in your project:

```
"com.fulcrumgenomics" %% "fgbio" % "0.1.3-SNAPSHOT"
```

## Overview

Fgbio is a command line tool to perform bioinformatic genomic data analysis. 
The collection of tools within `fgbio` are used by our customers for both ad-hoc data analysis and within their production pipelines.
These tools typically operate on read-level data (ex. FASTQ, SAM, or BAM) or variant-level data (ex. VCF or BCF).
They range from simple tools to filter reads in a BAM file, to tools to compute consensus reads from reads with the same molecular index/tag.
See the [list of tools](#list-of-tools) for more detail on the tools

## List of tools

Below we highlight a few tools that you may find useful.
Please see the help message for a full list of tools available.
In no particular order ...

* Tools to work with unique molecular tags/indexes (Umis). 
  * Annotating/Extract Umis from read-level data: `AnnotateBamWithUmis` and `ExtractUmisFromBam`.
  * Tools to manipulate read-level data containing Umis: `CallMolecularConsensusReads` and `GroupReadsByUmi`
* Tools to manipulate read-level data:
	* Filter read-level data: `FilterBam`.
	* Randomize the order of read-level data: `RandomizeBam`.
	* Update read-level metadata: `SetMateInformation` and `UpdateReadGroups`.
* Miscellaneous tools:
	* Pick molecular indices (ex. sample barcodes, or molecular indexes): `PickIlluminaIndices`.
	* Convert the output of HAPCUT (a tool for phasing variants): `HapCutToVcf`.
	* Find technical or synthetic sequences in read-level data: `FindTechnicalReads`.
    * Assess phased variant calls: `AssessPhasing`.


## Contributing

Contributions are welcome and encouraged.
We will do our best to provide an initial response to any pull request or issue within one-week.
For urgent matters, please contact us directly.

## Authors

* [Tim Fennell](https://github.com/tfenne) (maintainer)
* [Nils Homer](https://github.com/nh13) (maintainer)

## License

`fgbio` is open source software released under the [MIT License](https://github.com/fulcrumgenomics/fgbio/blob/master/LICENSE).

