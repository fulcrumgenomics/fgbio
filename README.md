[![Build Status](https://github.com/fulcrumgenomics/fgbio/workflows/unit%20tests/badge.svg)](https://github.com/fulcrumgenomics/fgbio/actions?query=workflow%3A%22unit+tests%22)
[![codecov](https://codecov.io/gh/fulcrumgenomics/fgbio/branch/main/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/fgbio)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/com.fulcrumgenomics/fgbio_2.13/badge.svg)](https://maven-badges.herokuapp.com/maven-central/com.fulcrumgenomics/fgbio_2.13)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/fgbio.svg?label=Bioconda)](http://bioconda.github.io/recipes/fgbio/README.html)
[![Javadocs](http://javadoc.io/badge/com.fulcrumgenomics/fgbio_2.13.svg)](http://javadoc.io/doc/com.fulcrumgenomics/fgbio_2.13)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/fgbio/blob/main/LICENSE)
[![Language](http://img.shields.io/badge/language-scala-brightgreen.svg)](http://www.scala-lang.org/)
[![DOI](https://zenodo.org/badge/53011104.svg)](https://zenodo.org/doi/10.5281/zenodo.10456900)

fgbio
====

A set of tools to analyze genomic data with a focus on Next Generation Sequencing.

<p>
<a href float="left"="https://fulcrumgenomics.com"><img src=".github/logos/fulcrumgenomics.svg" alt="Fulcrum Genomics" height="100"/></a>
</p>


[Visit us at Fulcrum Genomics](www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with fgbio and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>





This readme document is mostly for developers/contributors and those attempting to build the project from source.
Detailed user documentation is available on the [project website](http://fulcrumgenomics.github.io/fgbio/) including [tool usage](http://fulcrumgenomics.github.io/fgbio/tools/latest) and [documentation of metrics produced](http://fulcrumgenomics.github.io/fgbio/metrics/latest).  Detailed developer documentation can be found [here](http://javadoc.io/doc/com.fulcrumgenomics/fgbio_2.13).

<!---toc start-->
  * [Quick Installation](#quick-installation)
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

## Quick Installation

The [conda](https://conda.io/) package manager (configured with [bioconda channels](https://bioconda.github.io/)) can be used to quickly install fgbio:

```
conda install fgbio
```

To install fgbio without extra dependencies (e.g. [R](https://www.r-project.org/)), use the command:

```
conda install fgbio-minimal
```

## Goals

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

-   Tools for working with Unique Molecular Indexes (UMIs, aka Molecular IDs or Molecular Barcodes):
    -   Annotate/Extract Umis from read-level data: [`FastqToBam`][fgbio-fastqtobam-link], [`AnnotateBamWithUmis`][fgbio-annotatebamwithumis-link], [`ExtractUmisFromBam`][fgbio-extractumisfrombam-link], and [`CopyUmiFromReadName`][fgbio-copyumifromreadname-link].
    -   Manipulate read-level data containing Umis: [`CorrectUmis`][fgbio-correctumis-link], [`GroupReadsByUmi`][fgbio-groupreadsbyumi-link], [`CallMolecularConsensusReads`][fgbio-callmolecularconsensusreads-link], [`CallDuplexConsensusReads`][fgbio-callduplexconsensusreads-link], and [`FilterConsensusReads`][fgbio-filterconsensusreads-link].
    -   Collect metrics and review consensus reads: [`CollectDuplexSeqMetrics`][fgbio-collectduplexseqmetrics-link] and [`ReviewConsensusVariants`][fgbio-reviewconsensusvariants-link].
-   Tools to manipulate read-level data:
    -   Fastq Manipulation: [`FastqToBam`][fgbio-fastqtobam-link], [`ZipperBams`][fgbio-zipperbams-link], and [`DemuxFastqs`][fgbio-demuxfastqs-link] (see `[fqtk`][fqtk-link], our rust re-implementation for sample demultiplexing).
    -   Filter, clip, randomize, sort, and update metadata for read-level data: [`FilterBam`][fgbio-filterbam-link], [`ClipBam`][fgbio-clipbam-link], [`RandomizeBam`][fgbio-randomizebam-link], [`SortBam`][fgbio-sortbam-link], [`SetMateInformation`][fgbio-setmateinformation-link] and [`UpdateReadGroups`][fgbio-updatereadgroups-link].
-   Tools for quality control assessment:
    -   Detailed substitution error rate evaluation: [`ErrorRateByReadPosition`][fgbio-errorratebyreadposition-link].
    -   Sample pooling QC: [`EstimatePoolingFractions`]: [fgbio-estimatepoolingfractions-link].
    -   Splice-aware insert size QC for RNA-seq libraries: [`EstimateRnaSeqInsertSize`][fgbio-estimaternaseqinsertsize-link].
-   Tools for adding or manipulating alternate contig names:
    -   Extract contig names from an NCBI Assembly Report: [`CollectAlternateContigNames`][fgbio-collectalternatecontignames-link].
    -   Update contig names in common file formats: [`UpdateFastaContigNames`][fgbio-updatefastacontignames-link], [`UpdateVcfContigNames`][fgbio-updatevcfcontignames-link], [`UpdateGffContigNames`][fgbio-updategffcontignames-link], [`UpdateIntervalListContigNames`][fgbio-updateintervallistcontignames-link], [`UpdateDelimitedFileContigNames`][fgbio-updatedelimitedfilecontignames-link].
-   Miscellaneous tools:
    -   Pick molecular indices (ex. sample barcodes, or molecular indexes): [`PickIlluminaIndices`][fgbio-pickilluminaindices-link] and [`PickLongIndices`][fgbio-picklongindices-link].
    -   Find technical/synthetic, or switch-back sequences in read-level data: [`FindTechnicalReads`][fgbio-findtechnicalreads-link] and [`FindSwitchbackReads`][fgbio-findswitchbackreads-link].
    -   Make synthetic mixture VCFs: [`MakeMixtureVcf`][fgbio-makemixturevcf-link] and [`MakeTwoSampleMixtureVcf`][fgbio-maketwosamplemixturevcf-link].

[fgbio-fastqtobam-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/FastqToBam.html
[fgbio-annotatebamwithumis-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/AnnotateBamWithUmis.html
[fgbio-extractumisfrombam-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/ExtractUmisFromBam.html
[fgbio-copyumifromreadname-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/CopyUmiFromReadName.html
[fgbio-correctumis-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/CorrectUmis.html
[fgbio-groupreadsbyumi-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/GroupReadsByUmi.html
[fgbio-callmolecularconsensusreads-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/CallMolecularConsensusReads.html
[fgbio-callduplexconsensusreads-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/CallDuplexConsensusReads.html
[fgbio-filterconsensusreads-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/FilterConsensusReads.html
[fgbio-collectduplexseqmetrics-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/CollectDuplexSeqMetrics.html
[fgbio-reviewconsensusvariants-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/ReviewConsensusVariants.html
[fgbio-fastqtobam-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/FastqToBam.html
[fgbio-zipperbams-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/ZipperBams.html
[fgbio-demuxfastqs-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/DemuxFastqs.html
[fgbio-filterbam-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/FilterBam.html
[fgbio-clipbam-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/ClipBam.html
[fgbio-randomizebam-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/RandomizeBam.html
[fgbio-setmateinformation-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/SetMateInformation.html
[fgbio-updatereadgroups-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/UpdateReadGroups.html
[fgbio-collectalternatecontignames-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/CollectAlternateContigNames.html
[fgbio-updatefastacontignames-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/UpdateFastaContigNames.html
[fgbio-updatevcfcontignames-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/UpdateVcfContigNames.html
[fgbio-updategffcontignames-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/UpdateGffContigNames.html
[fgbio-updateintervallistcontignames-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/UpdateIntervalListContigNames.html
[fgbio-updatedelimitedfilecontignames-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/UpdateDelimitedFileContigNames.html
[fgbio-errorratebyreadposition-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/ErrorRateByReadPosition.html
[fgbio-estimatepoolingfractions-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/EstimatePoolingFractions.html
[fgbio-estimaternaseqinsertsize-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/EstimateRnaSeqInsertSize.html
[fgbio-pickilluminaindices-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/PickIlluminaIndices.html
[fgbio-picklongindices-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/PickLongIndices.html
[fgbio-findtechnicalreads-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/FastqToBam.html
[fgbio-sortbam-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/SortBam.html
[fgbio-makemixturevcf-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/MakeMixtureVcf.html
[fgbio-maketwosamplemixturevcf-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/MakeTwoSampleMixtureVcf.html
[fgbio-findswitchbackreads-link]: https://fulcrumgenomics.github.io/fgbio/tools/latest/FindSwitchbackReads.html

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

`fgbio` is open source software released under the [MIT License](https://github.com/fulcrumgenomics/fgbio/blob/main/LICENSE).

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

