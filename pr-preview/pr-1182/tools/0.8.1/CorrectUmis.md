---
title: CorrectUmis
---

# CorrectUmis

## Overview
**Group:** Unique Molecular Identifiers (UMIs)

Corrects UMIs stored in BAM files when a set of fixed UMIs is in use.  If the set of UMIs used in
an experiment is known and is a _subset_ of the possible randomers of the same length, it is possible
to error-correct UMIs prior to grouping reads by UMI.  This tool takes an input BAM with UMIs in a
tag (`RX` by default) and set of known UMIs (either on the command line or in a file) and produces:

  1. A new BAM with corrected UMIs in the same tag the UMIs were found in
  2. Optionally a set of metrics about the representation of each UMI in the set
  3. Optionally a second BAM file of reads whose UMIs could not be corrected within the specific parameters

All of the fixed UMIs must be of he same length, and all UMIs in the BAM file must also have the same
length.  Multiple UMIs that are concatenated with hyphens (e.g. `AACCAGT-AGGTAGA`) are split apart,
corrected individually and then re-assembled.  A read is accepted only if all the UMIs can be corrected.

Correction is controlled by two parameters that are applied per-UMI:

  1. _--max-mismatches_ controls how many mismatches (no-calls are counted as mismatches) are tolerated
         between a UMI as read and a fixed UMI.
  2. _--min-distance_ controls how many more mismatches the next best hit must have

For example, with two fixed UMIs `AAAAA` and `CCCCC` and `--max-mismatches=3` and `--min-distance=2` the
following would happen:

  - AAAAA would match to AAAAA
  - AAGTG would match to AAAAA with three mismatches because CCCCCC has six mismatches and 6 >= 3 + 2
  - AACCA would be rejected because it is 2 mismatches to AAAAA and 3 to CCCCCC and 3 <= 2 + 2

The set of fixed UMIs may be specified on the command line using `--umis umi1 umi2 ...` or via one or
more files of UMIs with a single sequence per line using `--umi-files umis.txt more_umis.txt`.  If there
are multiple UMIs per template, leading to hyphenated UMI tags, the values for the fixed UMIs should
be single, non-hyphenated UMIs (e.g. if a record has `RX:Z:ACGT-GGCA`, you would use `--umis ACGT GGCA`).

Records which have their UMIs corrected (i.e. the UMI is not identical to one of the expected UMIs but is
close enough to be corrected) will by default have their original UMI stored in the `OX` tag. This can be
disabled with the `--dont-store-original-umis` option.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input SAM or BAM file.|Required|1||
|output|o|PathToBam|Output SAM or BAM file.|Required|1||
|rejects|r|PathToBam|Reject BAM file to save unassigned reads.|Optional|1||
|metrics|M|FilePath|Metrics file to write.|Optional|1||
|max-mismatches|m|Int|Maximum number of mismatches between a UMI and an expected UMI.|Required|1||
|min-distance|d|Int|Minimum distance (in mismatches) to next best UMI.|Required|1||
|umis|u|String|Expected UMI sequences.|Optional|Unlimited||
|umi-files|U|FilePath|File of UMI sequences, one per line.|Optional|Unlimited||
|umi-tag|t|String|Tag in which UMIs are stored.|Optional|1|RX|
|dont-store-original-umis|x|Boolean|Don't store original UMIs upon correction.|Optional|1|false|

