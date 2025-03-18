---
title: PickIlluminaIndices
---

# PickIlluminaIndices

## Overview
**Group:** Utilities

Picks a set of molecular indices that should work well together.

The indexes to evaluate are generated randomly, unless the `--candidates` option is given.  The latter
specifies a path to a file with one index per line from which indices are picked.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|length|l|Int|The length of each barcode sequence.|Optional|1|8|
|indices|n|Int|The number of indices desired.|Required|1||
|edit-distance|e|Int|The minimum edit distance between two indices in the set.|Optional|1|3|
|output|o|FilePath|File to write indices to.|Required|1||
|allow-reverses||Boolean|Allow indices that are lexical reverses of one another|Optional|1|false|
|allow-reverse-complements||Boolean|Allow indices that are reverse complements of one another|Optional|1|false|
|allow-palindromes||Boolean|Allow indices that are palindromic (`bases == rev(bases)`).|Optional|1|false|
|max-homopolymer||Int|Reject indices with a homopolymer of greater than this length.|Optional|1|2|
|min-gc||Double|The minimum GC fraction for a barcode to be accepted.|Optional|1|0.0|
|max-gc||Double|The maximum GC fraction for a barcode to be accepted.|Optional|1|0.7|
|threads|t|Int|Number of threads to use.|Optional|1|4|
|vienna-rna-dir||DirPath|The installation directory for `ViennaRNA`.|Optional|1||
|min-delta-g||Double|The lowest acceptable secondary structure `deltaG`.|Optional|1|-10.0|
|adapters||String|The indexed adapter sequence into which the indices will be integrated.|Required|Unlimited|AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG|
|avoid-sequence||String|Sequences that should be avoided.  Any kmer of `length` that appears in these sequences and their reverse complements will be thrown out.|Required|Unlimited|AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG, CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGCAGACCGNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG, CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG|
|candidates||FilePath|The candidate indices from which to choose with one index per line (nothing else), otherwise generate all possible indices|Optional|1||

