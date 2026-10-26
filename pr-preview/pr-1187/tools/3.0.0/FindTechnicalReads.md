---
title: FindTechnicalReads
---

# FindTechnicalReads

## Overview
**Group:** SAM/BAM

Find reads that are from technical or synthetic sequences in a BAM file. Takes in
a BAM file, extracts the read pairs and fragment reads that are unmapped, and tests
them to see if they are likely generated from a technical sequence (e.g. adapter
dimer).

The identification of reads is done by testing the first N bases (controlled by the
match-length parameter) of each read against all sub-sequences of length N from the
technical sequences.  Sub-sequences are generated from both the sequences and the
reverse complement of the sequences, ignoring any sub-sequences that include `N`s.

By default the output BAM file will contain all reads that matched to a sub-sequence of the
technical sequences and, if the read is paired, the read's mate pair.  An option is
available to apply a tag to matched reads (--tag/-t), and if specified each matching
read will be tagged with the 0-based index of the sequence to which it matched. In
combination with tagging it is possible to output all reads (-a/--all-reads) which will
re-create the input BAM with the addition of tags on matching reads.

The default set of sequences include a range of different Illumina adapter sequences
with the sample index/barcode region masked to Ns.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input SAM or BAM file|Required|1||
|output|o|PathToBam|Output SAM or BAM file|Required|1||
|match-length|m|Int|The number of bases at the start of the read to match against.|Optional|1|15|
|max-errors|e|Int|The maximum number of errors in the matched region.|Optional|1|1|
|sequences|s|String|The set of technical sequences to look for.|Required|Unlimited|AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG, CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGCAGACCGNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG, CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG|
|all-reads|a|Boolean|Output all reads.|Optional|1|false|
|tag|t|String|Tag to set to indicate a read is a technical sequence.|Optional|1||

