---
title: FilterSomaticVcf
---

# FilterSomaticVcf

## Overview
**Group:** VCF/BCF

Applies one or more filters to a VCF of somatic variants. The VCF must contain genotype information
for the tumor sample.  If the VCF also contains genotypes for one or more other samples, the
`--sample` option must be provided to specify the sample whose genotypes to examine and whose reads
are present in the BAM file.

Various options are available for filtering the reads coming from the BAM file, including
`--min-mapping-quality`, `--min-base-quality` and `--paired-reads-only`.  The latter filters to only
paired end reads where both reads are mapped.  Reads marked as duplicates, secondary alignments
and supplemental alignments are all filtered out.

Each available filter may generate annotations in the `INFO` field of the output VCF and optionally,
if a threshold is specified, may apply one or more `FILTER`s to applicable variants.

## Available Filters

### End Repair Artifact Filter

The end repair artifact filter attempts to measure the probability that a variant is caused by
errors in the template generated during the end-repair and A-base addition steps that are common
to many Illumina library preparation protocols.  The artifacts occur if/when the end repair creates
a recessed 3' end which is subsequently and incorrectly filled in with As during A-base addition.
This presents specifically as errors to T at the beginning of reads (and in very short templates,
as matching errors to A at the ends of reads).

The filter creates the `ERAP` info attribute on SNVs with an A or T alternate allele, to record
the p-value for rejecting the possibility that the variant is due to an end repair artifact.

Two options are available:

* `--end-repair-distance`  allows control over how close to the ends of reads/templates errors can be
                           considered to be candidates for the artifact. Higher values decrease the
                           power of the test, so this should be set as low as possible given observed
                           errors.
* `--end-repair-p-value`   the p-value below which a filter should be applied. If no value is supplied
                           only the annotation is produced and no filtering is performed.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToVcf|Input VCF of somatic variant calls.|Required|1||
|output|o|PathToVcf|Output VCF of filtered somatic variants.|Required|1||
|bam|b|PathToBam|BAM file for the tumor sample.|Required|1||
|sample|s|String|Sample name in VCF if `> 1` sample present.|Optional|1||
|min-mapping-quality|m|Int|Minimum mapping quality for reads.|Optional|1|30|
|min-base-quality|q|Int|Minimum base quality.|Optional|1|20|
|paired-reads-only|p|Boolean|Use only paired reads mapped in pairs.|Optional|1|false|
|end-repair-distance||Int|Distance from end of read to implicate end repair.|Optional|1|2|
|end-repair-p-value||Double|Minimum acceptable p-value for end repair test.|Optional|1||

