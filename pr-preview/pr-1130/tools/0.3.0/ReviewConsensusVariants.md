---
title: ReviewConsensusVariants
---

# ReviewConsensusVariants

## Overview
**Group:** Unique Molecular Identifiers (UMIs)

Extracts data to make reviewing of variant calls from consensus reads easier. Creates
a list of variant sites from the input VCF (SNPs only) or IntervalList then extracts all
the consensus reads that do not contain a reference allele at the the variant sites, and
all raw reads that contributed to those consensus reads.  This will include consensus
reads that carry the alternate allele, a third allele, a no-call or a spanning
deletion at the variant site.

Reads are correlated between consensus and grouped BAMs using a molecule ID stored
in an optional attribute, `MI` by default.  In order to support paired molecule IDs
where two or more molecule IDs are related (e.g. see the Paired assignment strategy
in _GroupReadsByUmi_) the molecule ID is truncated at the last `/` if present
(e.g. `1/A => 1` and `2 => 2`).

Both input BAMs must be coordinate sorted and indexed.

A pair of output BAMs named `<output>.consensus.bam` and `<output>.grouped.bam` are created
with the relevant reads from each input BAM, and a review file `<output>.txt` is
created.  The review file contains details on each variant position along with detailed
information on each consensus read that supports the variant.  If the sample-name argument
is supplied and the input is VCF, genotype information for that sample will be retrieved.
If the sample-name isn't supplied and the VCF contains only a single sample then those
genotypes will be used.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|FilePath|Input VCF or IntervalList of variant locations.|Required|1||
|sample|s|String|Name of the sample being reviewed.|Optional|1||
|consensus-bam|c|PathToBam|BAM file of consensus reads used to call variants.|Required|1||
|grouped-bam|g|PathToBam|BAM file of grouped raw reads used to build consensuses.|Required|1||
|ref|r|PathToFasta|Reference fasta file.|Required|1||
|output|o|PathPrefix|Basename of output files to create.|Required|1||
|ignore-ns-in-consensus-reads|N|Boolean|Ignore N bases in the consensus reads.|Optional|1|false|
|maf|m|Double|Only output detailed information for variants at maf and below.|Optional|1|0.05|

