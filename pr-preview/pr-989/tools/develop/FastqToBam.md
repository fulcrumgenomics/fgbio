---
title: FastqToBam
---

# FastqToBam

## Overview
**Group:** FASTQ

Generates an unmapped BAM (or SAM or CRAM) file from fastq files.  Takes in one or more fastq files (optionally
gzipped), each representing a different sequencing read (e.g. R1, R2, I1 or I2) and can use a set of read
structures to allocate bases in those reads to template reads, sample indices, unique molecular indices, cell
barcodes, or to designate bases to be skipped over.

Read structures are made up of `<number><operator>` pairs much like the CIGAR string in BAM files. Five kinds of
operators are recognized:

1. `T` identifies a template read
2. `B` identifies a sample barcode read
3. `M` identifies a unique molecular index read
4. `C` identifies a cell barcode read
5. `S` identifies a set of bases that should be skipped or ignored

At most one `<number><operator>` pair in a read structure may use a `+` sign in place of the number to denote
"all remaining bases", and it may appear at any position in the read structure.  This is useful if, e.g., FASTQs
have been trimmed and contain reads of varying length.  For example to convert a paired-end run with an index
read and where the first 5 bases of R1 are a UMI and the next five bases should be skipped (e.g. because they
are monotemplate) you might specify:

```
--input r1.fq r2.fq i1.fq --read-structures 5M5S+T +T +B
```

Alternative if you know your reads are of fixed length you could specify:

```
--input r1.fq r2.fq i1.fq --read-structures 5M5S65T 75T 8B
```

For more information on read structures see the
[Read Structure Wiki Page](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures)

UMIs may be extracted from the read sequences, from the FASTQ header (either the read name or the read comment),
or both.  If UMI segments are present in the read structures those will be extracted.  UMIs in the FASTQ header
are extracted when one of the following (mutually exclusive) options is specified:

1. `--extract-umis-from-read-names`: the read name (the text before the first whitespace) is expected to be
   `:`-separated with 7 or 8 fields, with any UMI in the 8th field.  A read name with 7 fields, or an empty 8th
   field, has no UMI.  Any other number of fields is an error.
2. `--extract-umis-from-read-comment`: the read comment (the text after the first whitespace) is expected to be
   `:`-separated with exactly 4 fields, with any UMI in the 4th field.  Illumina documents the comment as
   `<read>:<is filtered>:<control number>:<index>`, where the last field is normally the sample index, so only
   use this option when the UMIs were written there instead.  A read with an empty 4th field has no UMI.  A read
   with no comment, or with any other number of fields, is an error.

For example, both of the following FASTQ headers yield the UMI `ACGTACGT-TTGCGGCT`, the first with
`--extract-umis-from-read-names` and the second with `--extract-umis-from-read-comment`:

```
@A00123:8:HFLV2DSXX:1:1101:1000:2000:ACGTACGT+TTGCGGCT 1:N:0:GATTACA
@A00123:8:HFLV2DSXX:1:1101:1000:2000 1:N:0:ACGTACGT+TTGCGGCT
```

In both cases the UMI is taken from the header of the first FASTQ given to `--input`.  Multiple UMIs in the
header must be separated by `+` and are joined with `-`, UMIs prefixed with `r` are reverse complemented, and
the UMI is upper-cased.  It is an error if the resulting UMI contains any character other than `ACGTN-`.
Neither option may be used with `--umi-qual-tag` as qualities are not available for UMIs in the header.  If UMIs
are present in both the header and the reads, the final UMIs are constructed by first taking the UMIs from the
header, then adding a hyphen, then the UMIs extracted from the reads.

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
|cell-tag|c|String|Tag in which to store the cellular barcodes.|Optional|1|CB|
|cell-qual-tag|C|String|Tag in which to store the cellular barcodes qualities.|Optional|1||
|store-sample-barcode-qualities|Q|Boolean|Store the sample barcode qualities in the QT Tag.|Optional|1|false|
|extract-umis-from-read-names|n|Boolean|Extract UMI(s) from read names and prepend to UMIs from reads.|Optional|1|false|
|extract-umis-from-read-comment||Boolean|Extract UMI(s) from the read comment and prepend to UMIs from reads.|Optional|1|false|
|read-group-id||String|Read group ID to use in the file header.|Optional|1|A|
|sample||String|The name of the sequenced sample.|Required|1||
|library||String|The name/ID of the sequenced library.|Required|1||
|barcode|b|String|Library or Sample barcode sequence.|Optional|1||
|platform||String|Sequencing Platform.|Optional|1|illumina|
|platform-unit||String|Platform unit (e.g. '<flowcell-barcode>.<lane>.<sample-barcode>')|Optional|1||
|platform-model||String|Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)|Optional|1||
|sequencing-center||String|The sequencing center from which the data originated|Optional|1||
|predicted-insert-size||Integer|Predicted median insert size, to insert into the read group header|Optional|1||
|description||String|Description of the read group.|Optional|1||
|comment||String|Comment(s) to include in the output file's header.|Optional|Unlimited||
|run-date||Iso8601Date|Date the run was produced, to insert into the read group header|Optional|1||

