---
title: PickLongIndices
---

# PickLongIndices

## Overview
**Group:** Utilities

Picks a set of molecular indices that have at least a given number of mismatches between
them. Whereas `PickIlluminaIndices` attempts to pick a near-optimal set of indices,
`PickLongIndices` implements a significantly more efficient method based on generation of
random indices that can generate a large set of satisfactory indices in a small amount of
time and memory even for index lengths `>> 10bp`.

Many options exist for controlling aspects of the indices picked, including length, edit
distance (mismatches only), gc range, homopolymer content, secondary structure etc.

Secondary structure is predicted using ViennaRNA's `RNAfold` in DNA mode. To enable structure
checking both the `--vienna-rna-dir` and `--adapters` must be specified.  Adapters must be
strings of A, C, G, and T with a single block of Ns (e.g. `ACGTNNNN` or `ACNNNNGT`).  At runtime
the Ns are replaced with indices, and `deltaG` of the index-containing sequence is calculated.

The number of indices requested may not be possible to produce given other constraints.
When this is the case the tool will output as many indices as possible, though less than
the requested number.  In such cases it may be useful to try different values for `--attempt`.
This parameter controls how many attempts are made to find the next valid index before
quitting and outputting the accumulated indices.  Higher values will yield incrementally more
indices but require significantly longer runtimes.

A file of existing indices may be provided. Existing indices must be of the same length as
requested indices and composed of A, C, G and Ts, but are subject to no other constraints.
Index picking will then built a set comprised of the existing indices, and new indices which
satisfy all constraints.  Existing indices are included in the generated output file.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|length|l|Int|The length of each index sequence.|Optional|1|8|
|number-of-indices|n|Int|The number of indices desired.|Required|1||
|edit-distance|e|Int|The minimum edit distance between two indices in the set.|Optional|1|3|
|output|o|FilePath|File to write indices to.|Required|1||
|allow-reverses||Boolean|Allow indices that are lexical reverses of one another|Optional|1|false|
|allow-reverse-complements||Boolean|Allow indices that are reverse complements of one another|Optional|1|false|
|allow-palindromes||Boolean|Allow indices that are palindromic (`index == revcomp(index)`).|Optional|1|false|
|max-homopolymer||Int|Reject indices with a homopolymer of greater than this length.|Optional|1|2|
|min-gc|g|Double|The minimum GC fraction for an index to be accepted.|Optional|1|0.2|
|max-gc|G|Double|The maximum GC fraction for an index to be accepted.|Optional|1|0.8|
|existing||FilePath|File of existing index sequences to integrate, one per line.|Optional|1||
|seed|s|Int|Random seed value.|Optional|1|1|
|attempts|a|Int|Attempts to pick the next index before quitting.|Optional|1|100000|
|vienna-rna-dir||DirPath|The installation directory for `ViennaRNA`.|Optional|1||
|temperature|t|Double|The temperature at which to predict secondary structure.|Optional|1|25.0|
|min-delta-g||Double|The lowest acceptable secondary structure `deltaG`.|Optional|1|-10.0|
|adapters||String|Adapter sequence(s) into which indices will be inserted.|Optional|Unlimited||
|avoid-sequence||String|Any index sequence that appears in an avoid sequence or its reverse complement will be discarded.|Required|Unlimited|AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACGCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG, CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGCAGACCGNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG, CTGTCTCTTATACACATCTCCGAGCCCACGAGACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG, AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT, AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG|

