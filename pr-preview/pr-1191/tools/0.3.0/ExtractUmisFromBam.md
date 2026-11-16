---
title: ExtractUmisFromBam
---

# ExtractUmisFromBam

## Overview
**Group:** SAM/BAM

Extracts unique molecular indexes from reads in a BAM file into tags.

Currently only unmapped reads are supported.

Only template bases will be retained as read bases (stored in the `SEQ` field) as specified by the read structure.

A read structure should be provided for each read of a template.  For example, paired end reads should have two
read structures specified.  The tags to store the molecular indices will be associated with the molecular index
segment(s) in the read structure based on the order specified.  If only one molecular index tag is given, then the
molecular indices will be concatenated and stored in that tag. Otherwise the number of molecular indices in the
read structure should match the number of tags given. In the resulting BAM file each end of a pair will contain
the same molecular index tags and values. Additionally, when multiple molecular indices are present the
`--single-tag` option may be used to write all indices, concatenated, to a single tag in addition to the tags
specified in `--molecular-index-tags`.

Optionally, the read names can be annotated with the molecular indices directly.  In this case, the read name
will be formatted `<NAME>+<UMIs1><UMIs2>` where `<UMIs1>` is the concatenation of read one's molecular indices.
Similarly for `<UMIs2>`.

Mapping information will not be adjusted, as such, this tool should not be used on reads that have been mapped since
it will lead to an BAM with inconsistent records.

The read structure describes the structure of a given read as one or more read segments. A read segment describes
a contiguous stretch of bases of the same type (ex. template bases) of some length and some offset from the start
of the read.  Read structures are made up of `<number><operator>` pairs much like the CIGAR string in BAM files.
Four kinds ofoperators are recognized:

1. `T` identifies a template read
2. `B` identifies a sample barcode read
3. `M` identifies a unique molecular index read
4. `S` identifies a set of bases that should be skipped or ignored

The last `<number><operator>` pair may be specified using a '+' sign instead of number to denote "all remaining
bases". This is useful if, e.g., fastqs have been trimmed and contain reads of varying length.

An example would be `10B3M7S100T` which describes 120 bases, with the first ten bases being a sample barcode,
bases 11-13 being a molecular index, bases 14-20 ignored, and bases 21-120 being template bases. See
[Read Structures](https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures) for more information.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathToBam|Output BAM file.|Required|1||
|read-structure|r|ReadStructure|The read structure, one per read in a template.|Required|Unlimited||
|molecular-barcode-tags|b|String|**[DEPRECATED]** SAM tags in which to store the molecular barcodes (one-per segment).|Optional|Unlimited||
|molecular-index-tags|t|String|SAM tag(s) in which to store the molecular indices.|Optional|Unlimited||
|single-tag|s|String|Single tag into which to concatenate all molecular indices.|Optional|1||
|annotate-read-names|a|Boolean|Annotate the read names with the molecular indices. See usage for more details.|Optional|1|false|
|clipping-attribute|c|String|The SAM tag with the position in read to clip adapters (e.g. `XT` as produced by Picard's `MarkIlluminaAdapters`).|Optional|1||

