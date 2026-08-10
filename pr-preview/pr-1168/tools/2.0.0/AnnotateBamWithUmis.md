---
title: AnnotateBamWithUmis
---

# AnnotateBamWithUmis

## Overview
**Group:** SAM/BAM

Annotates existing BAM files with UMIs (Unique Molecular Indices, aka Molecular IDs,
Molecular barcodes) from separate FASTQ files. Takes an existing BAM file and either
one FASTQ file with UMI reads or multiple FASTQs if there are multiple UMIs per template,
matches the reads between the files based on read names, and produces an output BAM file
where each record is annotated with an optional tag (specified by `attribute`) that
contains the read sequence of the UMI.  Trailing read numbers (`/1` or `/2`) are
removed from FASTQ read names, as is any text after whitespace, before matching.
If multiple UMI segments are specified (see `--read-structure`) across one or more FASTQs,
they are delimited in the same order as FASTQs are specified on the command line.
The delimiter is controlled by the `--delimiter` option.

The `--read-structure` option may be used to specify which bases in the FASTQ contain UMI
bases.  Otherwise it is assumed the FASTQ contains only UMI bases.

The `--sorted` option may be used to indicate that the FASTQ has the same reads and is
sorted in the same order as the BAM file.

At the end of execution, reports how many records were processed and how many were
missing UMIs. If any read from the BAM file did not have a matching UMI read in the
FASTQ file, the program will exit with a non-zero exit status.  The `--fail-fast` option
may be specified to cause the program to terminate the first time it finds a records
without a matching UMI.

In order to avoid sorting the input files, the entire UMI fastq file(s) is read into
memory. As a result the program needs to be run with memory proportional the size of
the (uncompressed) fastq(s).  Use the `--sorted` option to traverse the UMI fastq and BAM
files assuming they are in the same order.  More precisely, the UMI fastq file will be
traversed first, reading in the next set of BAM reads with same read name as the
UMI's read name.  Those BAM reads will be annotated.  If no BAM reads exist for the UMI,
no logging or error will be reported.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|The input SAM or BAM file.|Required|1||
|fastq|f|PathToFastq|Input FASTQ(s) with UMI reads.|Required|Unlimited||
|output|o|PathToBam|Output BAM file to write.|Required|1||
|attribute|t|String|The BAM attribute to store UMI bases in.|Optional|1|RX|
|qual-attribute|q|String|The BAM attribute to store UMI qualities in.|Optional|1||
|read-structure|r|ReadStructure|The read structure for the FASTQ, otherwise all bases will be used.|Required|Unlimited|+M|
|sorted|s|Boolean|Whether the FASTQ file is sorted in the same order as the BAM.|Optional|1|false|
|fail-fast||Boolean|If set, fail on the first missing UMI.|Optional|1|false|

