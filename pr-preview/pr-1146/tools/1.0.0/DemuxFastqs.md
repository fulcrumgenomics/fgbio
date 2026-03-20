---
title: DemuxFastqs
---

# DemuxFastqs

## Overview
**Group:** FASTQ

Performs sample demultiplexing on FASTQs.

The sample barcode for each sample in the sample sheet will be compared against the sample barcode bases extracted from
the FASTQs, to assign each read to a sample.  Reads that do not match any sample within the given error tolerance
will be placed in the 'unmatched' file.

The type of output is specified with the `--output-type` option, and can be BAM (`--output-type BamOnly`),
gzipped FASTQ (`--output-type FastqOnly`), or both (`--output-type Both`).

For BAM output, the output directory will contain one BAM file per sample in the sample sheet or metadata CSV file,
plus a BAM for reads that could not be assigned to a sample given the criteria.  The output file names will be the
concatenation of sample id, sample name, and sample barcode bases (expected not observed), delimited by `-`.  A
metrics file will also be output providing analogous information to the metric described
[SampleBarcodeMetric](http://fulcrumgenomics.github.io/fgbio/metrics/latest/#samplebarcodemetric).

For gzipped FASTQ output, one or more gzipped FASTQs per sample in the sample sheet or metadata CSV file will be
written to the output directory. For paired end data, the output will have the suffix `_R1.fastq.gz` and
`_R2.fastq.gz` for read one and read two respectively.  The sample barcode and molecular barcodes (concatenated)
will be appended to the read name and delimited by a colon.  If the `--illumina-standards` option is given, then
the output read names and file names will follow the
[Illumina standards described here](https://help.basespace.illumina.com/articles/tutorials/upload-data-using-web-uploader/).

The output base qualities will be standardized to Sanger/SAM format.

FASTQs and associated read structures for each read should be given:

- a single fragment read should have one FASTQ and one read structure
- paired end reads should have two FASTQs and two read structures
- a dual-index sample with paired end reads should have four FASTQs and four read structures given: two for the
  two index reads, and two for the template reads.

Read structures are made up of `<number><operator>` pairs much like the `CIGAR` string in BAM files. Four kinds of
operators are recognized:

1. `T` identifies a template read
2. `B` identifies a sample barcode read
3. `M` identifies a unique molecular index read
4. `S` identifies a set of bases that should be skipped or ignored

The last `<number><operator>` pair may be specified using a `+` sign instead of number to denote "all remaining
bases". This is useful if, e.g., fastqs have been trimmed and contain reads of varying length. Both reads must
have template bases.  Any molecular identifiers will be concatenated using
the `-` delimiter and placed in the given SAM record tag (`RX` by default).  Similarly, the sample barcode bases
from the given read will be placed in the `BC` tag.

Metadata about the samples should be given in either an Illumina Experiment Manager sample sheet or a metadata CSV
file.  Formats are described in detail below.

The read structures will be used to extract the observed sample barcode, template bases, and molecular identifiers
from each read.  The observed sample barcode will be matched to the sample barcodes extracted from the bases in
the sample metadata and associated read structures.

## Sample Sheet
The read group's sample id, sample name, and library id all correspond to the similarly named values in the
sample sheet.  Library id will be the sample id if not found, and the platform unit will be the sample name
concatenated with the sample barcode bases delimited by a `.`.

The sample section of the sample sheet should contain information related to each sample with the following columns:

  * Sample_ID:      The sample identifier unique to the sample in the sample sheet.
  * Sample_Name:    The sample name.
  * Library_ID:     The library Identifier.  The combination sample name and library identifier should be unique
                    across the samples in the sample sheet.
  * Description:    The description of the sample, which will be placed in the description field in the output BAM's
                    read group.  This column may be omitted.
  * Sample_Barcode: The sample barcode bases unique to each sample. The name of the column containing the sample barcode 
                    can be changed using the `--column-for-sample-barcode` option.  If the sample barcode is present 
                    across multiple reads (ex. dual-index, or inline in both reads of a pair), then the expected
                    barcode bases from each read should be concatenated in the same order as the order of the reads'
                    FASTQs and read structures given to this tool.

## Metadata CSV

In lieu of a sample sheet, a simple CSV file may be provided with the necessary metadata.  This file should
contain the same columns as described above for the sample sheet (`Sample_ID`, `Sample_Name`, `Library_ID`, and
`Description`).

## Example Command Line

As an example, if the sequencing run was 2x100bp (paired end) with two 8bp index reads both reading a sample
barcode, as well as an in-line 8bp sample barcode in read one, the command line would be

```
--inputs r1.fq i1.fq i2.fq r2.fq --read-structures 8B92T 8B 8B 100T \
    --sample-sheet SampleSheet.csv --metrics metrics.txt --output output_folder
```

## Illumina FASTQ File Naming

The output FASTQs will be written according to Illumina naming standards when the `--illumina-standards` option
is specified.  This requires that the input fastq files contain Illumina format read names with comments containing
`<ReadNum>:<FilterFlag>:<0>:<SampleNumber>` (no angle brackets).  The output create output files with Illumina
names:

1. The file extension will be `_R1_001.fastq.gz` for read one, and `_R2_001.fastq.gz` for read two (if paired end).
2. The per-sample output prefix will be `<SampleName>_S<SampleOrdinal>_L<LaneNumber>` (without angle brackets).

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|inputs|i|PathToFastq|One or more input fastq files each corresponding to a sub-read (ex. index-read, read-one, read-two, fragment).|Required|Unlimited||
|output|o|DirPath|The output directory in which to place sample BAMs.|Required|1||
|metadata|x|FilePath|A file containing the metadata about the samples.|Required|1||
|read-structures|r|ReadStructure|The read structure for each of the FASTQs.|Required|Unlimited||
|metrics|m|FilePath|The file to which per-barcode metrics are written.  If none given, a file named `demux_barcode_metrics.txt` will be written to the output directory.|Optional|1||
|column-for-sample-barcode|c|String|The column name in the sample sheet or metadata CSV for the sample barcode.|Optional|1|Sample_Barcode|
|unmatched|u|String|Output BAM file name for the unmatched records.|Optional|1|unmatched.bam|
|quality-format|q|QualityEncoding|A value describing how the quality values are encoded in the FASTQ. Either Solexa for pre-pipeline 1.3 style scores (solexa scaling + 66), Illumina for pipeline 1.3 and above (phred scaling + 64) or Standard for phred scaled scores with a character shift of 33.  If this value is not specified, the quality format will be detected automatically.|Optional|1||
|threads|t|Int|The number of threads to use while de-multiplexing. The performance does not increase linearly with the # of threads and seems not to improve beyond 2-4 threads.|Optional|1|1|
|max-mismatches||Int|Maximum mismatches for a barcode to be considered a match.|Optional|1|1|
|min-mismatch-delta||Int|Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match.|Optional|1|2|
|max-no-calls||Int|Maximum allowable number of no-calls in a barcode read before it is considered unmatchable.|Optional|1|2|
|sort-order||SortOrder|The sort order for the output sam/bam file (typically unsorted or queryname).|Optional|1|queryname|
|umi-tag||String|The SAM tag for any molecular barcode.  If multiple molecular barcodes are specified, they will be concatenated and stored here.|Optional|1|RX|
|platform-unit||String|The platform unit (typically `<flowcell-barcode>-<sample-barcode>.<lane>`)|Optional|1||
|sequencing-center||String|The sequencing center from which the data originated|Optional|1||
|predicted-insert-size||Integer|Predicted median insert size, to insert into the read group header|Optional|1||
|platform-model||String|Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)|Optional|1||
|comments||String|Comment(s) to include in the merged output file's header.|Optional|Unlimited||
|run-date||Iso8601Date|Date the run was produced, to insert into the read group header|Optional|1||
|output-type||OutputType|The type of outputs to produce.|Optional|1||
|illumina-standards||Boolean|Output FASTQs according to Illumina naming standards, for example, for upload to the BaseSpace Sequence Hub|Optional|1|false|
|include-all-bases-in-fastqs||Boolean|Output all bases (i.e. all sample barcode, molecular barcode, skipped,         and template bases) for every read with template bases (ex. read one         and read two) as defined by the corresponding read structure(s).|Optional|1|false|

