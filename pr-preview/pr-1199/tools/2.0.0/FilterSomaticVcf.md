---
title: FilterSomaticVcf
---

# FilterSomaticVcf

## Overview
**Group:** VCF/BCF

Applies one or more filters to a VCF of somatic variants. The VCF must contain genotype information for the
tumor sample. If the VCF also contains genotypes for one or more other samples, the `--sample` option must be
provided to specify the sample whose genotypes to examine and whose reads are present in the BAM file.

Various options are available for filtering the reads coming from the BAM file, including
`--min-mapping-quality`, `--min-base-quality` and `--paired-reads-only`. The latter filters to only paired end
reads where both reads are mapped. Reads marked as duplicates, secondary alignments and supplemental alignments
are all filtered out.

Each available filter may generate annotations in the `INFO` field of the output VCF and optionally, if a
threshold is specified, may apply one or more `FILTER`s to applicable variants.

In previous versions of this tool, the only available filter was specific to A-base addition artifacts and was
referred to as the 'End Repair Artifact Filter.' This filter has been renamed to 'A-tailing Artifact Filter', but
its functionality is unchanged. The filter's associated command-line parameters, `INFO` field key, and `FILTER`
tag have also been renamed accordingly, as described below.

## Available Filters

### A-tailing Artifact Filter (previously 'End Repair Artifact Filter')

The A-tailing artifact filter attempts to measure the probability that a single-nucleotide mismatch is the
product of errors in the template generated during the A-base addition steps that are common to many Illumina
library preparation protocols. The artifacts occur if/when a recessed 3' end is incorrectly filled in with one\
or more adenines during A-base addition. Incorrect adenine incorporation presents specifically as errors to T at
the beginning of reads (and in very short templates, as matching errors to A at the ends of reads).

The filter adds the `INFO` field `ATAP` (previously `ERAP`) to SNVs with an A or T alternate allele. This field
records the p-value representing the probability of the null hypothesis that the variant is a true mutation, so
lower p-values indicate that the variant is more likely an A-tailing artifact. If a threshold p-value is
specified, the `FILTER` tag `ATailingArtifact` (previously `EndRepairArtifact`) will be applied to variants with
p-values less than or equal to the threshold.

Two options are available:

* `--a-tailing-distance`  (previously `--end-repair-distance`) allows control over how close to the ends of
                          reads/templates errors can be considered to be candidates for the A-tailing artifact.
                          Higher values decrease the power of the test, so this should be set as low as possible
                          given observed errors.
* `--a-tailing-p-value`   (previously `--end-repair-p-value`) the p-value at or below which a filter should be
                          applied. If no value is supplied only the `INFO` annotation is produced and no `FILTER`
                          is applied.

### End Repair Fill-in Artifact Filter

The end repair fill-in artifact filter attempts to measure the probability that a single-nucleotide mismatch is
the product of an error in the template generated during the end repair fill-in step that is common to many
Illumina library preparation protocols, in which single-stranded 3' overhangs are filled in to create a blunt
end. These artifacts originate from single-stranded templates containing damaged bases, often as a consequence
of oxidative damage. These DNA lesions, for example 8-oxoguanine, undergo mismatched pairing, which after PCR
appear as mutations at the ends of reads.

The filter adds the `INFO` field `ERFAP` to records SNVs. This field records the p-value representing the
probability of the null hypothesis (e.g. that the variant is a true mutation), so lower p-values indicate that
the variant is more likely an end repair fill-in artifact. If a threshold p-value is specified, then the `FILTER`
tag `EndRepairFillInArtifact` will be applied to variants with p-values less than or equal to the threshold.

Two options are available:

* `--end-repair-fill-in-distance`  allows control over how close to the ends of reads/templates errors can be
                                   considered to be candidates for the artifact. Higher values decrease the
                                   power of the test, so this should be set as low as possible given observed
                                   errors.
* `--end-repair-fill-in-p-value`   the p-value below which a filter should be applied. If no value is supplied
                                   only the annotation is produced and no filtering is performed.

## Performance Expectations

By default `--access-pattern` will be set to `RandomAccess` and the input BAM will be queried using index-based
random access. Random access is mandatory if the input VCF is not coordinate sorted. If random access is not
requested and the input VCF is not coordinate sorted, then an exception will be raised on the first
non-coordinate increasing VCF record found. The BAM must be coordinate sorted in all cases and additionally be
indexed if random access is requested.

Often, a VCF file will contain a sparse set of records that are scattered across a given territory within a
genome (or the records will be sparsely scattered genome-wide). If the territory of the VCF records is markedly
smaller than the territory of all aligned SAM records in the BAM file, then random access may be the most
efficient BAM access pattern. However, there are cases where random access will be less efficient such as when
the VCF is coordinate sorted and the variant call records are very densely packed across a similar territory as
compared to all aligned SAM records. Such a case is common in deeply sequenced hybrid selection NGS experiments
and setting `--access-pattern` to `Streaming` will often be the most efficient BAM access pattern.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToVcf|Input VCF of somatic variant calls.|Required|1||
|output|o|PathToVcf|Output VCF of filtered somatic variants.|Required|1||
|bam|b|PathToBam|BAM file for the tumor sample.|Required|1||
|sample|s|String|Sample name in VCF if `> 1` sample present.|Optional|1||
|min-mapping-quality|m|Int|Minimum mapping quality for reads.|Optional|1|20|
|min-base-quality|q|Int|Minimum base quality.|Optional|1|20|
|paired-reads-only|p|Boolean|Use only paired reads mapped in pairs.|Optional|1|false|
|access-pattern|A|BamAccessPattern|The type of BAM access pattern to use.|Optional|1|RandomAccess|
|a-tailing-distance||Int|Distance from 5-prime end of read to implicate A-base addition artifacts. Set to :none: to deactivate the filter.|Optional|1|2|
|a-tailing-p-value||Double|Minimum acceptable p-value for the A-base addition artifact test.|Optional|1||
|end-repair-fill-in-distance||Int|Distance from 5-prime end of read to implicate end repair fill-in artifacts. Set to :none: to deactivate the filter.|Optional|1|15|
|end-repair-fill-in-p-value||Double|Minimum acceptable p-value for the end repair fill-in artifact test.|Optional|1||

