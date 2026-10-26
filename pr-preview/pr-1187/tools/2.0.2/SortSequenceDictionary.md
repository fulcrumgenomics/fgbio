---
title: SortSequenceDictionary
---

# SortSequenceDictionary

## Overview
**Group:** FASTA

Sorts a sequence dictionary file in the order of another sequence dictionary.

The inputs are to two `*.dict` files.  One to be sorted, and the other to provide the order for the sorting.

If there is a contig in the input dictionary that is not in the sorting dictionary, that contig will be appended
to the end of the sequence dictionary in the same relative order to other appended contigs as in the input dictionary.
Missing contigs can be omitted by setting `--skip-missing-contigs` to true.

If there is a contig in the sorting dictionary that is not in the input dictionary, that contig will be ignored.

The output will be a sequence dictionary, containing the version header line and one
line per contig.  The fields of the entries in this dictionary will be the same as in input, but in the order of
`--sort-dictionary`.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToSequenceDictionary|Input sequence dictionary file to be sorted.|Required|1||
|sort-dictionary|d|PathToSequenceDictionary|Input sequence dictionary file containing contigs in the desired sort order.|Required|1||
|output|o|PathToSequenceDictionary|Output sequence dictionary file.|Required|1||
|skip-missing-contigs||Boolean|Skip input contigs that have no matching contig in the sort dictionary rather than appending to the end of the output dictionary.|Optional|1|false|

