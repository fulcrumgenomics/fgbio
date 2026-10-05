---
title: FastqToBam
---

# FastqToBam

## Overview
**Group:** FASTQ

Generates an unmapped BAM (or SAM or CRAM) file from fastq files.  Takes in one or more fastq files (optionally
gzipped), each representing a different sequencing read (e.g. R1, R2, I1 or I2) and can use a set of read
structures to allocate bases in those reads to template reads, sample indices, unique molecular indices, or to
designate bases to be skipped over.

Read structures are made up of `<number><operator>` pairs much like the CIGAR string in BAM files. Four kinds of
operators are recognized:

1. `T` identifies a template read
2. `B` identifies a sample barcode read
3. `M` identifies a unique molecular index read
4. `S` identifies a set of bases that should be skipped or ignored

The last `<number><operator>` pair may be specified using a `+` sign instead of number to denote "all remaining
bases". This is useful if, e.g., fastqs have been trimmed and contain reads of varying length.  For example
to convert a paired-end run with an index read and where the first 5 bases of R1 are a UMI and the second
five bases are monotemplate you might specify:

```
--input r1.fq r2.fq i1.fq --read-structures 5M5S+T +T +B
```

Alternative if you know your reads are of fixed length you could specify:

```
--input r1.fq r2.fq i1.fq --read-structures 5M5S65T 75T 8B
```

For more information on read structures see the
[Read Structure Wiki Page](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)

The same number of input files and read structures must be provided, with one exception: if supplying exactly
1 or 2 fastq files, both of which are solely template reads, no read structures need be provided.

The output file can be sorted by queryname using the `--sort-order` option; the default is to produce a BAM
with reads in the same order as they appear in the fastq file.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToFastq|Fastq files corresponding to each sequencing read (e.g. R1, I1, etc.).|Required|Unlimited||
|output|o|PathToBam|The output SAM or BAM file to be written.|Required|1||
|read-structures|r|ReadStructure|Read structures, one for each of the FASTQs.|Optional|Unlimited||
|sort|s|Boolean|If true, queryname sort the BAM file, otherwise preserve input order.|Optional|1|false|
|umi-tag|u|String|Tag in which to store molecular barcodes/UMIs.|Optional|1|RX|
|umi-qual-tag|q|String|Tag in which to store molecular barcode/UMI qualities.|Optional|1||
|read-group-id||String|Read group ID to use in the file header.|Optional|1|A|
|sample||String|The name of the sequenced sample.|Required|1||
|library||String|The name/ID of the sequenced library.|Required|1||
|platform||String|Sequencing Platform.|Optional|1|illumina|
|platform-unit||String|Platform unit (e.g. '<flowcell-barcode>.<lane>.<sample-barcode>')|Optional|1||
|platform-model||String|Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)|Optional|1||
|sequencing-center||String|The sequencing center from which the data originated|Optional|1||
|predicted-insert-size||Integer|Predicted median insert size, to insert into the read group header|Optional|1||
|description||String|Description of the read group.|Optional|1||
|comment||String|Comment(s) to include in the output file's header.|Optional|Unlimited||
|run-date||Iso8601Date|Date the run was produced, to insert into the read group header|Optional|1||

