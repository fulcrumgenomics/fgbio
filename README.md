[![Build Status](https://github.com/fulcrumgenomics/fgbio/workflows/unit%20tests/badge.svg)](https://github.com/fulcrumgenomics/fgbio/actions?query=workflow%3A%22unit+tests%22)
[![codecov](https://codecov.io/gh/fulcrumgenomics/fgbio/branch/master/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/fgbio)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/fc4f5fe8dbe34bf784114435b202fab4)](https://www.codacy.com/app/contact_32/fgbio?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=fulcrumgenomics/fgbio&amp;utm_campaign=Badge_Grade)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/com.fulcrumgenomics/fgbio_2.13/badge.svg)](https://maven-badges.herokuapp.com/maven-central/com.fulcrumgenomics/fgbio_2.13)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/fgbio.svg?label=Bioconda)](http://bioconda.github.io/recipes/fgbio/README.html)
[![Javadocs](http://javadoc.io/badge/com.fulcrumgenomics/fgbio_2.13.svg)](http://javadoc.io/doc/com.fulcrumgenomics/fgbio_2.13)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/fgbio/blob/master/LICENSE)
[![Language](http://img.shields.io/badge/language-scala-brightgreen.svg)](http://www.scala-lang.org/)

fgbio
====

A set of tools to analyze genomic data with a focus on Next Generation Sequencing.  This readme document is mostly for developers/contributors and those attempting to build the project from source.
Detailed user documentation is available on the [project website](http://fulcrumgenomics.github.io/fgbio/) including [tool usage](http://fulcrumgenomics.github.io/fgbio/tools/latest) and [documentation of metrics produced](http://fulcrumgenomics.github.io/fgbio/metrics/latest).  Detailed developer documentation can be found [here](http://javadoc.io/doc/com.fulcrumgenomics/fgbio_2.13).

<!---toc start-->
  * [Goals](#goals)
  * [Overview](#overview)
  * [List of tools](#list-of-tools)
  * [Building](#building)
  * [Command line](#command-line)
  * [Include fgbio in your project](#include-fgbio-in-your-project)
  * [Contributing](#contributing)
  * [Authors](#authors)
  * [License](#license)
  * [Sponsorship](#sponsorship)

<!---toc end-->

# Goals

There are many toolkits available for analyzing genomic data; fgbio does not aim to be all things to all people but is specifically focused on providing:

* Robust, well-tested tools.
* An easy to use command-line.
* Clear and thorough documentation for each tool.
* Open source development for the benefit of the community and our clients.

## Overview

Fgbio is a set of command line tools to perform bioinformatic/genomic data analysis. 
The collection of tools within `fgbio` are used by our customers and others both for ad-hoc data analysis and within production pipelines.
These tools typically operate on read-level data (ex. FASTQ, SAM, or BAM) or variant-level data (ex. VCF or BCF).
They range from simple tools to filter reads in a BAM file, to tools to compute consensus reads from reads with the same molecular index/tag.
See the [list of tools](#list-of-tools) for more detail on the tools

## List of tools

For a full list of available tools please see the [tools section](http://fulcrumgenomics.github.io/fgbio/tools/latest) of the project website.

Below we highlight a few tools that you may find useful.

* Tools for working with Unique Molecular Indexes (UMIs, aka Molecular IDs or MIDs). 
  * Annotating/Extract Umis from read-level data: `AnnotateBamWithUmis` and `ExtractUmisFromBam`.
  * Tools to manipulate read-level data containing Umis: `CorrectUmis`, `GroupReadsByUmi`, `CallMolecularConsensusReads` and `CallDuplexConsensusReads`
* Tools to manipulate read-level data:
  * FastqManipulation: `DemuxFastqs` and `FastqToBam`
  * Filter read-level data: `FilterBam`.
  * Clipping of reads: `ClipBam`.
  * Randomize the order of read-level data: `RandomizeBam`.
  * Update read-level metadata: `SetMateInformation` and `UpdateReadGroups`.
* Quality assessment tools:
  * Detailed substitution error rate evaluation: `ErrorRateByReadPosition`
  * Sample pooling QC: `EstimatePoolingFractions`
  * Splice-aware insert size QC for RNA-seq libraries: `EstimateRnaSeqInsertSize`
  * Assessment of duplex sequencing experiments: `CollectDuplexSeqMetrics`
* Miscellaneous tools:
  * Pick molecular indices (ex. sample barcodes, or molecular indexes): `PickIlluminaIndices` and `PickLongIndices`.
  * Convert the output of HAPCUT (a tool for phasing variants): `HapCutToVcf`.
  * Find technical or synthetic sequences in read-level data: `FindTechnicalReads`.
  * Assess phased variant calls: `AssessPhasing`.

## Building 
### Cloning the Repository

[Git LFS](https://git-lfs.github.com/) is used to store large files used in testing fgbio.  In order to compile and run tests it is necessary to [install git lfs](https://git-lfs.github.com/).  To retrieve the large files either:

1. Clone the repository _after_ installing git lfs, or
2. In a previously cloned repository run `git lfs pull` once

After initial setup regular git commands (e.g. `pull`, `fetch`, `push`) will also operate on large files and no special handling is needed.

To clone the repository: `git clone https://github.com/fulcrumgenomics/fgbio.git`

### Running the build
fgbio is built using [sbt](http://www.scala-sbt.org/).

Use ```sbt assembly``` to build an executable jar in ```target/scala-2.13/```.

Tests may be run with ```sbt test```.

Java SE 8 is required.


## Command line

`java -jar target/scala-2.13/fgbio-<version>.jar` to see the commands supported.  Use `java -jar target/scala-2.13/fgbio-<version>.jar <command>` to see the help message for a particular command.

## Include fgbio in your project

You can include `fgbio` in your project using:

```
"com.fulcrumgenomics" %% "fgbio" % "1.0.0"
```

for the latest released version or (buyer beware):

```
"com.fulcrumgenomics" %% "fgbio" % "0.9.0-<commit-hash>-SNAPSHOT"
```

for the latest development snapshot.

## Contributing

Contributions are welcome and encouraged.
We will do our best to provide an initial response to any pull request or issue within one-week.
For urgent matters, please contact us directly.

## Authors

* [Tim Fennell](https://github.com/tfenne) (maintainer)
* [Nils Homer](https://github.com/nh13) (maintainer)

## License

`fgbio` is open source software released under the [MIT License](https://github.com/fulcrumgenomics/fgbio/blob/master/LICENSE).

## Sponsorship

### Become a sponsor

As a free and open source project, `fgbio` relies on the support of the community of users for its development. If you work for an organization that uses and benefits from `fgbio`, please consider supporting `fgbio`. There are different ways, such as employing people to work on `fgbio`, funding the project, or becoming a [sponsor](https://github.com/sponsors/fulcrumgenomics) to support the broader ecosystem. Please [contact@fulcrumgenomics.com](https://www.fulcrumgenomics.com/contact/) to discuss.

### Sponsors

Sponsors provide support for `fgbio` through direct funding or employing contributors.
Public sponsors include:

<p>
<a href float="left"="https://fulcrumgenomics.com"><img src=".github/logos/fulcrumgenomics.svg" alt="Fulcrum Genomics" height="35"/></a>
&nbsp;
<a href float="left"="https://twinstrandbio.com/"><img src=".github/logos/twinstrandbio.svg" alt="TwinStrand Biosciences" height="45"/></a>
&nbsp;
<a href float="left"="https://www.jumpcodegenomics.com//"><img src=".github/logos/jumpcodegenomics.png" alt="Jumpcode Genomics" height="30"/></a>
&nbsp;
<a href float="left"="https://www.igenomx.com//"><img src=".github/logos/igenomx.png" alt="iGenomX" height="30"/></a>
&nbsp;
<a href float="left"="https://myriad.com"><img src=".github/logos/myriad.png" alt="Myriad Genetics" height="35"/></a>
&nbsp;
<a href float="left"="https://missionbio.com"><img src=".github/logos/missionbio.svg" alt="Mission Bio" height="30"/></a>
&nbsp;
<a href float="left"="https://singulargenomics.com"><img src=".github/logos/singulargenomics.svg" alt="Singular Genomics" height="30"/></a>
&nbsp;
<a href float="left"="https://verogen.com"><img src=".github/logos/verogen.jpg" alt="Verogen" height="30"/></a>
&nbsp;
<a href float="left"="https://idtdna.com"><img src=".github/logos/idtdna.png" alt="Integrated DNA Technologies" height="30"/></a>
&nbsp;
<a href float="left"="https://strataoncology.com"><img src=".github/logos/strataoncology.png" alt="Strata Oncology" height="30"/></a>
</p>

The full list of sponsors supporting `fgbio` is available in the [sponsor](https://github.com/sponsors/fulcrumgenomics) page.

