---
title: EstimateRnaSeqInsertSize
---

# EstimateRnaSeqInsertSize

## Overview
**Group:** RNA-Seq

Computes the insert size for RNA-Seq experiments.

Computes the insert size by counting the # of bases sequenced in transcript space.  The insert size is defined
as the distance between the first bases sequenced of each pair respectively (5' sequencing ends).

This tool skips reads that overlap multiple genes, reads that aren't fully enclosed in a gene, and reads where the
insert size would disagree across transcripts from the same gene.  Also skips reads that are unpaired, failed QC,
secondary, supplementary, pairs without both ends mapped, duplicates, and pairs whose reads map to different
chromosomes. Finally, skips transcripts where too few mapped read bases overlap exonic sequence.

This tool requires each mapped pair to have the mate cigar (`MC`) tag.  Use `SetMateInformation` to add the mate cigar.

The output metric file will have the extension `.rna_seq_insert_size.txt` and the output histogram file will have
the extension `.rna_seq_insert_size_histogram.txt`.  The histogram file gives for each orientation (`FR`, `RF`, `tandem`),
the number of read pairs that had the given insert size.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|ref-flat|r|FilePath|Input gene annotations in [RefFlat](http://genome.ucsc.edu/goldenPath/gbdDescriptionsOld.html#RefFlat) form|Required|1||
|prefix|p|PathPrefix|Output prefix file.  The file will have the extension `.rna_seq_insert_size.txt` if not given|Optional|1||
|include-duplicates|d|Boolean|Include duplicates|Optional|1|false|
|deviations|D|Double|Generate mean and standard deviation by filtering to `median + deviations*median_absolute_deviation`. This is done because insert size data typically includes enough anomalous values from chimeras and other artifacts to make the mean and sd grossly misleading regarding the real distribution.   "|Optional|1|10.0|
|minimum-mapping-quality|q|Int|Ignore reads with mapping quality less than this value.|Optional|1|30|
|minimum-overlap|m|Double|The minimum fraction of read bases that must overlap exonic sequence in a transcript|Optional|1|0.95|

