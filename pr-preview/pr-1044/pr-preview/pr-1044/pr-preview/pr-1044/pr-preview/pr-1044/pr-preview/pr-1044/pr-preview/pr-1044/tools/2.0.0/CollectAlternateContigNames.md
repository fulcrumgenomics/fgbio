---
title: CollectAlternateContigNames
---

# CollectAlternateContigNames

## Overview
**Group:** FASTA

Collates the alternate contig names from an NCBI assembly report.

The input is to be the `*.assembly_report.txt` obtained from NCBI.

The output will be a "sequence dictionary", which is a valid SAM file, containing the version header line and one
line per contig.  The primary contig name (i.e. `@SQ.SN`) is specified with `--primary` option, while alternate
names (i.e. aliases) are specified with the `--alternates` option.

The `Assigned-Molecule` column, if specified as an `--alternate`, will only be used for sequences with
`Sequence-Role` `assembled-molecule`.

When updating an existing sequence dictionary with `--existing` the primary contig names must match.  I.e. the
contig name from the assembly report column specified by `--primary` must match the contig name in the existing
sequence dictionary (`@SQ.SN`).  All contigs in the existing sequence dictionary must be present in the assembly
report.  Furthermore, contigs in the assembly report not found in the sequence dictionary will be ignored.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|FilePath|Input NCBI assembly report file.|Required|1||
|output|o|PathToSequenceDictionary|Output sequence dictionary file.|Required|1||
|primary|p|AssemblyReportColumn|The assembly report column for the primary contig name.|Optional|1|RefSeqAccession|
|alternates|a|AssemblyReportColumn|The assembly report column(s) for the alternate contig name(s)|Required|Unlimited||
|sequence-roles|s|SequenceRole|Only output sequences with the given sequence roles.  If none given, all sequences will be output.|Optional|Unlimited||
|existing|d|PathToSequenceDictionary|Update an existing sequence dictionary file.  The primary names must match.|Optional|1||
|allow-mismatching-lengths|x|Boolean|Allow mismatching sequence lengths when using an existing sequence dictionary file.|Optional|1|false|
|skip-missing-alternates||Boolean|Skip contigs that have no alternates|Optional|1|true|

