---
title: CallMolecularConsensusReads
---

# CallMolecularConsensusReads

## Overview
**Group:** Unique Molecular Identifiers (UMIs)

Calls consensus sequences from reads with the same unique molecular tag.

Reads with the same unique molecular tag are examined base-by-base to assess the likelihood of each base in the
source molecule.  The likelihood model is as follows:

1. First, the base qualities are adjusted. The base qualities are assumed to represent the probability of a
   sequencing error (i.e. the sequencer observed the wrong base present on the cluster/flowcell/well). The base
   quality scores are converted to probabilities incorporating a probability representing the chance of an error
   from the time the unique molecular tags were integrated to just prior to sequencing.  The resulting probability
   is the error rate of all processes from right after integrating the molecular tag through to the end of
   sequencing.
2. Next, a consensus sequence is called for all reads with the same unique molecular tag base-by-base.  For a
   given base position in the reads, the likelihoods that an A, C, G, or T is the base for the underlying
   source molecule respectively are computed by multiplying the likelihood of each read observing the base
   position being considered.  The probability of error (from 1.) is used when the observed base does not match
   the hypothesized base for the underlying source molecule, while one minus that probability is used otherwise.
   The computed likelihoods are normalized by dividing them by the sum of all four likelihoods to produce a
   posterior probability, namely the probability that the source molecule was an A, C, G, or T from just after
   integrating molecular tag through to sequencing, given the observations.  The base with the maximum posterior
   probability as the consensus call, and the posterior probability is used as its raw base quality.
3. Finally, the consensus raw base quality is modified by incorporating the probability of an error prior to
   integrating the unique molecular tags.  Therefore, the probability used for the final consensus base
   quality is the posterior probability of the source molecule having the consensus base given the observed
   reads with the same molecular tag, all the way from sample extraction and through sample and library
   preparation, through preparing the library for sequencing (e.g. amplification, target selection), and finally,
   through sequencing.

This tool assumes that reads with the same tag are grouped together (consecutive in the file). Also, this tool
calls each end of a pair independently, and does not jointly call bases that overlap within a pair.  Insertion or
deletion errors in the reads are not considered in the consensus model.

Particular attention should be paid to setting the `--min-reads` parameter as this can have a dramatic effect on
both results and runtime.  For libraries with low duplication rates (e.g. 100-300X exomes libraries) in which it
is desirable to retain singleton reads while making consensus reads from sets of duplicates, `--min-reads=1` is
appropriate.  For libraries with high duplication rates where it is desirable to only produce consensus reads
supported by 2+ reads to allow error correction, `--min-reads=2` or higher is appropriate.  After generation,
consensus reads can be further filtered using the _FilterConsensusReads_ tool.  As such it is always safe to run
with `--min-reads=1` and filter later, but filtering at this step can improve performance significantly.

Consensus reads have a number of additional optional tags set in the resulting BAM file.  The tags break down into
those that are single-valued per read:

```
consensus depth      [cD] (int)  : the maximum depth of raw reads at any point in the consensus read
consensus min depth  [cM] (int)  : the minimum depth of raw reads at any point in the consensus read
consensus error rate [cE] (float): the fraction of bases in raw reads disagreeing with the final consensus calls
```

And those that have a value per base:

```
consensus depth  [cd] (short[]): the count of bases contributing to the consensus read at each position
consensus errors [ce] (short[]): the number of bases from raw reads disagreeing with the final consensus base
```

The per base depths and errors are both capped at 32,767. In all cases no-calls (`N`s) and bases below the
`--min-input-base-quality` are not counted in tag value calculations.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|The input SAM or BAM file.|Required|1||
|output|o|PathToBam|Output SAM or BAM file to write consensus reads.|Required|1||
|rejects|r|PathToBam|Optional output SAM or BAM file to write reads not used.|Optional|1||
|tag|t|String|The SAM attribute with the unique molecule tag.|Optional|1|MI|
|read-name-prefix|p|String|The Prefix all consensus read names|Optional|1||
|read-group-id|R|String|The new read group ID for all the consensus reads.|Optional|1|A|
|error-rate-pre-umi|1|PhredScore|The Phred-scaled error rate for an error prior to the UMIs being integrated.|Optional|1|45|
|error-rate-post-umi|2|PhredScore|The Phred-scaled error rate for an error post the UMIs have been integrated.|Optional|1|40|
|min-input-base-quality|m|PhredScore|Ignore bases in raw reads that have Q below this value.|Optional|1|10|
|min-consensus-base-quality|N|PhredScore|Deprecated: will be removed in future versions; use FilterConsensusReads to filter consensus bases on quality instead. Mask (make 'N') consensus bases with quality less than this threshold.|Optional|1|2|
|min-reads|M|Int|The minimum number of reads to produce a consensus base.|Required|1||
|max-reads||Int|The maximum number of reads to use when building a consensus. If more than this many reads are present in a tag family, the family is randomly downsampled to exactly max-reads reads.|Optional|1||
|output-per-base-tags|B|Boolean|If true produce tags on consensus reads that contain per-base information.|Optional|1|true|
|sort-order|S|SamOrder|The sort order of the output, the same as the input if not given.|Optional|1||
|debug|D|Boolean|Turn on debug logging.|Optional|1|false|
|threads||Int|The number of threads to use while consensus calling.|Optional|1|1|

