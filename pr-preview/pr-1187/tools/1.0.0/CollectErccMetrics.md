---
title: CollectErccMetrics
---

# CollectErccMetrics

## Overview
**Group:** RNA-Seq

Collects metrics for ERCC spike-ins for RNA-Seq experiments.

Currently calculates per-transcript ERCC metrics and summarizes dose response, but does not calculate fold-change
response.

The input BAM should contain reads mapped to a reference containing the ERCC transcripts.  The reference may have
additional contigs, for example, when concatenating the sample's reference genome and the ERCC transcripts.  The
BAM should have sequence lines in the header matching the ERCC ids (ex. ERCC-00130 or ERCC-00004).

The standard ERCC transcripts, including their unique IDs and concentrations, are taken from
[here](https://assets.thermofisher.com/TFS-Assets/LSG/manuals/cms_095046.txt).  The second column lists the ERCC
transcript identifier which should be present as a sequence line in the BAM header, while columns four and five
give the concentration of the transcript for mixtures #1 and #2.

The choice of mixture to use can be specified with the `--mixture-name` option (either `Mix1` or `Mix2`), or with
a file containing a custom list of transcripts and concentrations using the `--custom-mixture` option as follows.
The custom mixture file should be tab-delimited file containing the following columns:
  1. ERCC ID - each ERCC ID should match a contig/reference-sequence name in the input SAM/BAM header.
  2. Concentration - the concentration (in `attomoles/ul`).
The custom mixture file should contain a header line with names for each column, though the actual values will be ignored.

Three outputs will be produced:
  1. <output>.ercc_summary_metrics.txt - summary statistics for total # of reads mapping to the ERCC transcripts and dose
                                response metrics.
  2. <output>.ercc_detailed_metrics.txt - gives a per-ERCC-transcript expected concentration and observed fragment count.
  3. <output>.ercc_metrics.pdf - plots the expected concentration versus the observed fragment count.

Secondary andsupplementary reads will be ignored.  A read pair mapping to an ERCC transcript is counted only if both
ends of the pair map to the same ERCC transcript.  A minimum mapping quality can be required for reads to be
counted as mapped to an ERCC transcript.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input SAM or BAM file of aligned reads in coordinate order.|Required|1||
|output|o|PathPrefix|Output prefix.|Required|1||
|mixture-name|m|ErccMixture|The name of the standard ERCC mixture.|Optional|1||
|custom-mixture||FilePath|Tab-delimited file containing ERCC IDs and expected concentrations.|Optional|1||
|min-transcript-count|c|Int|Minimum # of counts required to include an ERCC transcript.|Optional|1|3|
|minimum-mapping-quality|M|Int|The minimum mapping quality|Optional|1|10|

