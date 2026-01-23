---
title: CallOverlappingConsensusBases
---

# CallOverlappingConsensusBases

## Overview
**Group:** SAM/BAM

Consensus calls overlapping bases in read pairs.

## Inputs and Outputs

In order to correctly correct reads by template, the input BAM must be either `queryname` sorted or `query` grouped.
The sort can be done in streaming fashion with:

```
samtools sort -n -u in.bam | fgbio CallOverlappingConsensusBases -i /dev/stdin ...
```

The output sort order may be specified with `--sort-order`.  If not given, then the output will be in the same
order as input.

The reference FASTA must be given so that any existing `NM`, `UQ` and `MD` tags can be repaired.

## Correction

Only mapped read pairs with overlapping bases will be eligible for correction.

Each read base from the read and its mate that map to same position in the reference will be used to create
a consensus base as follows:

1. If the base agree, then the chosen agreement strategy (`--agreement-strategy`) will be used.
2. If the base disagree, then the chosen disagreement strategy (`--disagreement-strategy`) will be used.

The agreement strategies are as follows:

* Consensus:   Call the consensus base and return a new base quality that is the sum of the two base qualities.
* MaxQual:     Call the consensus base and return a new base quality that is the maximum of the two base qualities.
* PassThrough: Leave the bases and base qualities unchanged.

In the context of disagreement strategies, masking a base will make the base an "N" with base quality phred-value "2".
The disagreement strategies are as follows:

* MaskBoth:      Mask both bases.
* MaskLowerQual: Mask the base with the lowest base quality, with the other base unchanged.  If the base qualities
                 are the same, mask both bases.
* Consensus:     Consensus call the base.  If the base qualities are the same, mask both bases.  Otherwise, call the
                 base with the highest base quality and return a new base quality that is the difference between the
                 highest and lowest base quality.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input SAM or BAM file of aligned reads.|Required|1||
|output|o|PathToBam|Output SAM or BAM file.|Required|1||
|metrics|m|FilePath|Output metrics file.|Required|1||
|ref|r|PathToFasta|Reference sequence fasta file.|Required|1||
|threads||Int|The number of threads to use while consensus calling.|Optional|1|1|
|sort-order|S|SamOrder|The sort order of the output. If not given, output will be in the same order as input if the input.|Optional|1||
|agreement-strategy||AgreementStrategy|The strategy to consensus call when both bases agree.  See the usage for more details|Optional|1|Consensus|
|disagreement-strategy||DisagreementStrategy|The strategy to consensus call when both bases disagree.  See the usage for more details|Optional|1|Consensus|

