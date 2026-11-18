
# fgbio Metrics Descriptions

This page contains descriptions of all metrics produced by all fgbio tools.  Within the descriptions
the type of each field/column is given, including two commonly used types:

* `Count` is an integer representing the count of some item
* `Proportion` is a real number with a value between 0 and 1 representing a proportion or fraction

## Table of Contents

|Metric Type|Description|
|-----------|-----------|
|[PhaseBlockLengthMetric](#phaseblocklengthmetric)|Metrics produced by `AssessPhasing` describing the number of phased blocks of a given length|
|[AssessPhasingMetric](#assessphasingmetric)|Metrics produced by `AssessPhasing` describing various statistics assessing the performance of phasing variants relative to a known set of phased variant calls|
|[SampleBarcodeMetric](#samplebarcodemetric)|Metrics for matching templates to sample barcodes primarily used in com|
|[TagFamilySizeMetric](#tagfamilysizemetric)|Metrics produced by `GroupReadsByUmi` to describe the distribution of tag family sizes observed during grouping|
|[ConsensusVariantReviewInfo](#consensusvariantreviewinfo)|Detailed information produced by `ReviewConsensusVariants` on variants called in consensus reads|
|[UmiCorrectionMetrics](#umicorrectionmetrics)|Metrics produced by `CorrectUmis` regarding the correction of UMI sequences to a fixed set of known UMIs|
|[DuplexUmiMetric](#duplexumimetric)|Metrics produced by `CollectDuplexSeqMetrics` describing the set of observed duplex UMI sequences and the frequency of their observations|
|[UmiMetric](#umimetric)|Metrics produced by `CollectDuplexSeqMetrics` describing the set of observed UMI sequences and the frequency of their observations|
|[DuplexYieldMetric](#duplexyieldmetric)|Metrics produced by `CollectDuplexSeqMetrics` that are sampled at various levels of coverage, via random downsampling, during the construction of duplex metrics|
|[DuplexFamilySizeMetric](#duplexfamilysizemetric)|Metrics produced by `CollectDuplexSeqMetrics` to describe the distribution of double-stranded (duplex) tag families in terms of the number of reads observed on each strand|
|[FamilySizeMetric](#familysizemetric)|Metrics produced by `CollectDuplexSeqMetrics` to quantify the distribution of different kinds of read family sizes|
|[InsertSizeMetric](#insertsizemetric)|Metrics produced by `EstimateRnaSeqInsertSize` to describe the distribution of insert sizes within an RNA-seq experiment|
|[ErccSummaryMetrics](#erccsummarymetrics)|Metrics produced by `CollectErccMetrics` describing various summary metrics related to the spike-in of ERCC (External RNA Controls Consortium) into an RNA-Seq experiment|
|[ErccDetailedMetric](#erccdetailedmetric)|Metrics produced by `CollectErccMetrics` describing various per-transcript metrics related to the spike-in of ERCC (External RNA Controls Consortium) into an RNA-Seq experiment|
|[RunInfo](#runinfo)|Stores the result of parsing the run info (RunInfo|
|[PoolingFractionMetric](#poolingfractionmetric)|Metrics produced by `EstimatePoolingFractions` to quantify the estimated proportion of a sample mixture that is attributable to a specific sample with a known set of genotypes|
|[ErrorRateByReadPositionMetric](#errorratebyreadpositionmetric)|Metrics produced by `ErrorRateByReadPosition` describing the number of base observations and substitution errors at each position within each sequencing read|
|[ClippingMetrics](#clippingmetrics)|Metrics produced by ClipBam that detail how many reads and bases are clipped respectively|

## Metric File Descriptions


### PhaseBlockLengthMetric

Metrics produced by `AssessPhasing` describing the number of phased blocks of a given length.  The output will have
multiple rows, one for each observed phased block length.


|Column|Type|Description|
|------|----|-----------|
|dataset|String|The name of the dataset being assessed (i.e. "truth" or "called").|
|length|Long|The length of the phased block.|
|count|Long|The number of phased blocks of the given length.|


### AssessPhasingMetric

Metrics produced by `AssessPhasing` describing various statistics assessing the performance of phasing variants
relative to a known set of phased variant calls.  Included are methods for assessing sensitivity and accuracy from
a number of previous papers (ex. http://dx.doi.org/10.1038%2Fng.3119).The N50, N90, and L50 statistics are defined as follows:
- The N50 is the longest block length such that the bases covered by all blocks this length and longer are at least
50% of the # of bases covered by all blocks.
- The N90 is the longest block length such that the bases covered by all blocks this length and longer are at least
90% of the # of bases covered by all blocks.
- The L50 is the smallest number of blocks such that the sum of the lengths of the blocks is `>=` 50% of the sum of
the lengths of all blocks.


|Column|Type|Description|
|------|----|-----------|
|num_called|Long|The number of variants called.|
|num_phased|Long|The number of variants called with phase.|
|num_truth|Long|The number of variants with known truth genotypes.|
|num_truth_phased|Long|The number of variants with known truth genotypes with phase.|
|num_called_with_truth_phased|Long|The number of variants called that had a known phased genotype.|
|num_phased_with_truth_phased|Long|The number of variants called with phase that had a known phased genotype.|
|num_truth_phased_in_called_block|Long|The number of known phased variants that were in a called phased block.|
|num_both_phased_in_called_block|Long|The number of called phase variants that had a known phased genotype in a called phased block.|
|num_short_switch_errors|Long|The number of short switch errors (isolated switch errors).|
|num_long_switch_errors|Long|The number of long switch errors (# of runs of consecutive switch errors).|
|num_switch_sites|Long|The number of sites that could be (short or long) switch errors (i.e. the # of sites with both known and called phased variants).|
|num_illumina_point_switch_errors|Long|The number of point switch errors (defined in http://dx.doi.org/10.1038%2Fng.3119).|
|num_illumina_long_switch_errors|Long|The number of long switch errors (defined in http://dx.doi.org/10.1038%2Fng.3119).|
|num_illumina_switch_sites|Long|The number of sites that could be (point or long) switch errors (defined in http://dx.doi.org/10.1038%2Fng.3119).|
|frac_phased|Double|The fraction of called variants with phase.|
|frac_phased_with_truth_phased|Double|The fraction of known phased variants called with phase.|
|frac_truth_phased_in_called_block|Double|The fraction of phased known genotypes in a called phased block.|
|frac_phased_with_truth_phased_in_called_block|Double|The fraction of called phased variants that had a known phased genotype in a called phased block.|
|short_accuracy|Double|The fraction of switch sites without short switch errors (`1 - (num_short_switch_errors / num_switch_sites)`).|
|long_accuracy|Double|The fraction of switch sites without long switch errors (`1 - (num_long_switch_errors / num_switch_sites)`).|
|illumina_point_accuracy|Double|The fraction of switch sites without point switch errors according to the Illumina                                method defining switch sites and errors (`1 - (num_illumina_point_switch_errors / num_illumina_switch_sites )`).|
|illumina_long_accuracy|Double|The fraction of switch sites wihtout long switch errors  according to the Illumina                               method defining switch sites and errors (`1 - (num_illumina_long_switch_errors / num_illumina_switch_sites )`).|
|mean_called_block_length|Double|The mean phased block length in the callset.|
|median_called_block_length|Double|The median phased block length in the callset.|
|stddev_called_block_length|Double|The standard deviation of the phased block length in the callset.|
|n50_called_block_length|Double|The N50 of the phased block length in the callset.|
|n90_called_block_length|Double|The N90 of the phased block length in the callset.|
|l50_called|Double|The L50  of the phased block length in the callset.|
|mean_truth_block_length|Double|The mean phased block length in the truth.|
|median_truth_block_length|Double|The median phased block length in the truth.|
|stddev_truth_block_length|Double|The standard deviation of the phased block length in the truth.|
|n50_truth_block_length|Double|The N50 of the phased block length in the truth.|
|n90_truth_block_length|Double|The N90 of the phased block length in the callset.|
|l50_truth|Double|The L50 of the phased block length in the callset.|


### SampleBarcodeMetric

Metrics for matching templates to sample barcodes primarily used in com.fulcrumgenomics.fastq.DemuxFastqs.The number of templates will match the number of reads for an Illumina single-end sequencing run, while the number
of templates will be half the number of reads for an Illumina paired-end sequencing run (i.e. R1 & R2 observe the
same template).


|Column|Type|Description|
|------|----|-----------|
|barcode_name|String|The name for the sample barcode, typically the sample name from the SampleSheet.|
|library_name|String|The name of the library, typically the library identifier from the SampleSheet.|
|barcode|String|The sample barcode bases.  Dual index barcodes will have two sample barcode sequences delimited by a                dash.|
|templates|Count|The total number of templates matching the given barcode.|
|pf_templates|Count|The total number of pass-filter templates matching the given barcode.|
|perfect_matches|Count|The number of templates that match perfectly the given barcode.|
|pf_perfect_matches|Count|The number of pass-filter templates that match perfectly the given barcode.|
|one_mismatch_matches|Count|The number of pass-filter templates that match the given barcode with exactly one                             mismatch.|
|pf_one_mismatch_matches|Count|The number of pass-filter templates that match the given barcode with exactly                                one mismatch.|
|fraction_matches|Proportion|The fraction of all templates that match the given barcode.|
|ratio_this_barcode_to_best_barcode|Proportion|The rate of all templates matching this barcode to all template                                               reads matching the most prevalent barcode. For the most prevalent                                               barcode this will be 1, for all others it will be less than 1 (except                                               for the possible exception of when there are more unmatched templates                                               than for any other barcode, in which case the value may be arbitrarily                                               large).  One over the lowest number in this column gives you the                                               fold-difference in representation between barcodes.|
|pf_fraction_matches|Proportion|The fraction of all pass-filter templates that match the given barcode.|
|pf_ratio_this_barcode_to_best_barcode|Proportion|The rate of all pass-filter templates matching this barcode to                                                  all templates matching the most prevalent barcode. For the                                                  most prevalent barcode this will be 1, for all others it will be                                                  less than 1 (except for the possible exception of when there are                                                  more unmatched templates than for any other barcode, in which                                                  case the value may be arbitrarily large).  One over the lowest                                                  number in this column gives you the fold-difference in                                                  representation between barcodes.|
|pf_normalized_matches|Proportion|The "normalized" matches to each barcode. This is calculated as the number of                              pass-filter templates matching this barcode over the mean of all pass-filter                              templates matching any barcode (excluding unmatched). If all barcodes are                              represented equally this will be                              1.|


### TagFamilySizeMetric

Metrics produced by `GroupReadsByUmi` to describe the distribution of tag family sizes
observed during grouping.


|Column|Type|Description|
|------|----|-----------|
|family_size|Int|The family size, or number of templates/read-pairs belonging to the family.|
|count|Count|The number of families (or source molecules) observed with `family_size` observations.|
|fraction|Proportion|The fraction of all families of all sizes that have this specific `family_size`.|
|fraction_gt_or_eq_family_size|Proportion|The fraction of all families that have `>= family_size`.|


### ConsensusVariantReviewInfo

Detailed information produced by `ReviewConsensusVariants` on variants called in consensus reads. Each
row contains information about a consensus _read_ that carried a variant or non-reference allele at a
particular variant site.The first 10 columns (up to `N`) contain information about the variant site and are repeated for each
consensus read reported at that site.  The remaining fields are specific to the consensus read.


|Column|Type|Description|
|------|----|-----------|
|chrom|String|The chromosome on which the variant exists.|
|pos|Int|The position of the variant.|
|ref|String|The reference allele at the position.|
|genotype|String|The genotype of the sample in question.|
|filters|String|The set of filters applied to the variant in the VCF.|
|A|Int|The count of A observations at the variant locus across all consensus reads.|
|C|Int|The count of C observations at the variant locus across all consensus reads.|
|G|Int|The count of G observations at the variant locus across all consensus reads.|
|T|Int|The count of T observations at the variant locus across all consensus reads.|
|N|Int|The count of N observations at the variant locus across all consensus reads.|
|consensus_read|String|The consensus read name for which the following fields contain values.|
|consensus_insert|String|A description of the insert that generated the consensus read.|
|consensus_call|Char|The base call from the consensus read.|
|consensus_qual|Int|The quality score from the consensus read.|
|a|Int|The number of As in raw-reads contributing to the consensus base call at the variant site.|
|c|Int|The number of Cs in raw-reads contributing to the consensus base call at the variant site.|
|g|Int|The number of Gs in raw-reads contributing to the consensus base call at the variant site.|
|t|Int|The number of Ts in raw-reads contributing to the consensus base call at the variant site.|
|n|Int|The number of Ns in raw-reads contributing to the consensus base call at the variant site.|


### UmiCorrectionMetrics

Metrics produced by `CorrectUmis` regarding the correction of UMI sequences to a fixed set of known UMIs.


|Column|Type|Description|
|------|----|-----------|
|umi|String|The corrected UMI sequence (or all `N`s for unmatched).|
|total_matches|Count|The number of UMI sequences that matched/were corrected to this UMI.|
|perfect_matches|Count|The number of UMI sequences that were perfect matches to this UMI.|
|one_mismatch_matches|Count|The number of UMI sequences that matched with a single mismatch.|
|two_mismatch_matches|Count|The number of UMI sequences that matched with two mismatches.|
|other_matches|Count|The number of UMI sequences that matched with three or more mismatches.|
|fraction_of_matches|Proportion|The fraction of all UMIs that matched or were corrected to this UMI.|
|representation|Double|The `total_matches` for this UMI divided by the _mean_ `total_matches` for all UMIs.|


### DuplexUmiMetric

Metrics produced by `CollectDuplexSeqMetrics` describing the set of observed duplex UMI sequences and the
frequency of their observations.  The UMI sequences reported may have been corrected using information
within a double-stranded tag family.  For example if a tag family is comprised of three read pairs with
UMIs `ACGT-TGGT`, `ACGT-TGGT`, and `ACGT-TGGG` then a consensus UMI of `ACGT-TGGT` will be generated.UMI pairs are normalized within a tag family so that observations are always reported as if they came
from a read pair with read 1 on the positive strand (F1R2). Another way to view this is that for FR or RF
read pairs, the duplex UMI reported is the UMI from the positive strand read followed by the UMI from the
negative strand read.  E.g. a read pair with UMI `AAAA-GGGG` and with R1 on the negative strand and R2 on
the positive strand, will be reported as `GGGG-AAAA`.


|Column|Type|Description|
|------|----|-----------|
|umi|String|The duplex UMI sequence, possibly-corrected.|
|raw_observations|Count|The number of read pairs in the input BAM that observe the duplex UMI (after correction).|
|raw_observations_with_errors|Count|The subset of raw observations that underwent any correction.|
|unique_observations|Count|The number of double-stranded tag families (i.e unique double-stranded molecules)                            that observed the duplex UMI.|
|fraction_raw_observations|Proportion|The fraction of all raw observations that the duplex UMI accounts for.|
|fraction_unique_observations|Proportion|The fraction of all unique observations that the duplex UMI accounts for.|
|fraction_unique_observations_expected|Proportion|The fraction of all unique observations that are expected to be                                              attributed to the duplex UMI based on the `fraction_unique_observations`                                              of the two individual UMIs.|


### UmiMetric

Metrics produced by `CollectDuplexSeqMetrics` describing the set of observed UMI sequences and the
frequency of their observations.  The UMI sequences reported may have been corrected using information
within a double-stranded tag family.  For example if a tag family is comprised of three read pairs with
UMIs `ACGT-TGGT`, `ACGT-TGGT`, and `ACGT-TGGG` then a consensus UMI of `ACGT-TGGT` will be generated,
and three raw observations counted for each of `ACGT` and `TGGT`, and no observations counted for `TGGG`.


|Column|Type|Description|
|------|----|-----------|
|umi|String|The UMI sequence, possibly-corrected.|
|raw_observations|Count|The number of read pairs in the input BAM that observe the UMI (after correction).|
|raw_observations_with_errors|Count|The subset of raw-observations that underwent any correction.|
|unique_observations|Count|The number of double-stranded tag families (i.e unique double-stranded molecules)                            that observed the UMI.|
|fraction_raw_observations|Proportion|The fraction of all raw observations that the UMI accounts for.|
|fraction_unique_observations|Proportion|The fraction of all unique observations that the UMI accounts for.|


### DuplexYieldMetric

Metrics produced by `CollectDuplexSeqMetrics` that are sampled at various levels of coverage, via random
downsampling, during the construction of duplex metrics.  The downsampling is done in such a way that the
`fraction`s are approximate, and not exact, therefore the `fraction` field should only be interpreted as a guide
and the `read_pairs` field used to quantify how much data was used.See `FamilySizeMetric` for detailed definitions of `CS`, `SS` and `DS` as used below.


|Column|Type|Description|
|------|----|-----------|
|fraction|Proportion|The approximate fraction of the full dataset that was used to generate the remaining values.|
|read_pairs|Count|The number of read pairs upon which the remaining metrics are based.|
|cs_families|Count|The number of _CS_ (Coordinate & Strand) families present in the data.|
|ss_families|Count|The number of _SS_ (Single-Strand by UMI) families present in the data.|
|ds_families|Count|The number of _DS_ (Double-Strand by UMI) families present in the data.|
|ds_duplexes|Count|The number of _DS_ families that had the minimum number of observations on both strands to be                    called duplexes (default = 1 read on each strand).|
|ds_fraction_duplexes|Proportion|The fraction of _DS_ families that are duplexes (`ds_duplexes / ds_families`).|
|ds_fraction_duplexes_ideal|Proportion|The fraction of _DS_ families that should be duplexes under an idealized model                                   where each strand, `A` and `B`, have equal probability of being sampled, given                                   the observed distribution of _DS_ family sizes.|


### DuplexFamilySizeMetric

Metrics produced by `CollectDuplexSeqMetrics` to describe the distribution of double-stranded (duplex)
tag families in terms of the number of reads observed on each strand.We refer to the two strands as `ab` and `ba` because we identify the two strands by observing the same pair of
UMIs (A and B) in opposite order (A->B vs B->A). Which strand is `ab` and which is `ba` is largely arbitrary, so
to make interpretation of the metrics simpler we use a definition here that for a given tag family
`ab` is the sub-family with more reads and `ba` is the tag family with fewer reads.


|Column|Type|Description|
|------|----|-----------|
|ab_size|Int|The number of reads in the `ab` sub-family (the larger sub-family) for this double-strand tag family.|
|ba_size|Int|The number of reads in the `ba` sub-family (the smaller sub-family) for this double-strand tag family.|
|count|Count|The number of families with the `ab` and `ba` single-strand families of size `ab_size` and `ba_size`.|
|fraction|Proportion|The fraction of all double-stranded tag families that have `ab_size` and `ba_size`.|
|fraction_gt_or_eq_size|Proportion|The fraction of all double-stranded tag families that have                               `ab reads >= ab_size` and `ba reads >= ba_size`.|


### FamilySizeMetric

Metrics produced by `CollectDuplexSeqMetrics` to quantify the distribution of different kinds of read family
sizes.  Three kinds of families are described:1. _CS_ or _Coordinate & Strand_: families of reads that are grouped together by their unclipped 5'
   genomic positions and strands just as they are in traditional PCR duplicate marking
2. _SS_ or _Single Strand_: single-strand families that are each subsets of a CS family create by
   also using the UMIs to partition the larger family, but not linking up families that are
   created from opposing strands of the same source molecule.
3. _DS_ or _Double Strand_: families that are created by combining single-strand families that are from
   opposite strands of the same source molecule. This does **not** imply that all DS families are composed
   of reads from both strands; where only one strand of a source molecule is observed a DS family is
   still counted.


|Column|Type|Description|
|------|----|-----------|
|family_size|Int|The family size, i.e. the number of read pairs grouped together into a family.|
|cs_count|Count|The count of families with `size == family_size` when grouping just by coordinates and strand information.|
|cs_fraction|Proportion|The fraction of all _CS_ families where `size == family_size`.|
|cs_fraction_gt_or_eq_size|Proportion|The fraction of all _CS_ families where `size >= family_size`.|
|ss_count|Count|The count of families with `size == family_size` when also grouping by UMI to create single-strand families.|
|ss_fraction|Proportion|The fraction of all _SS_ families where `size == family_size`.|
|ss_fraction_gt_or_eq_size|Proportion|The fraction of all _SS_ families where `size >= family_size`.|
|ds_count|Count|The count of families with `size == family_size`when also grouping by UMI and merging single-strand                 families from opposite strands of the same source molecule.|
|ds_fraction|Proportion|The fraction of all _DS_ families where `size == family_size`.|
|ds_fraction_gt_or_eq_size|Proportion|The fraction of all _DS_ families where `size >= family_size`.|


### InsertSizeMetric

Metrics produced by `EstimateRnaSeqInsertSize` to describe the distribution of insert sizes within an
RNA-seq experiment.  The insert sizes are computed in "transcript space", accounting for spliced
alignments, in order to get a true estimate of the size of the DNA fragment, not just it's span on
the genome.


|Column|Type|Description|
|------|----|-----------|
|pair_orientation|PairOrientation|The orientation of the reads within a read-pair relative to each other.                         Possible values are FR, RF and TANDEM.|
|read_pairs|Long|The number of read pairs observed with the `pair_orientation`.|
|mean|Double|The mean insert size of the read pairs.|
|standard_deviation|Double|The standard deviation of the insert size of the read pairs.|
|median|Double|The median insert size of the read pairs.|
|min|Long|The minimum observed insert size of the read pairs.|
|max|Long|The maximum observed insert size of the read pairs.|
|median_absolute_deviation|Double|The median absolution deviation of the read pairs.|


### ErccSummaryMetrics

Metrics produced by `CollectErccMetrics` describing various summary metrics related to the spike-in of ERCC
(External RNA Controls Consortium) into an RNA-Seq experiment.The correlation coefficients and linear regression are calculated based on the log2 observed read pair count normalized
by ERCC transcript length versus the log2 expected concentration.


|Column|Type|Description|
|------|----|-----------|
|total_reads|Long|The total number of reads considered.|
|ercc_reads|Long|The total number of reads mapping to an ERCC transcript.|
|fraction_ercc_reads|Double|The fraction of total reads that map to an ERCC transcript.|
|ercc_templates|Long|The total number of read pairs (or single end reads) mapping to an ERCC transcript.|
|total_transcripts|Int|The total number of ERCC transcripts with at least one read observed.|
|passing_filter_transcripts|Int|The total number of ERCC transcripts with at least the user-set minimum # of reads observed.|
|pearsons_correlation|Option[Double]|Pearson's correlation coefficient for correlation of concentration and normalized counts.|
|spearmans_correlation|Option[Double]|Spearman's correlation coefficient for correlation of concentration and normalized counts.|
|intercept|Option[Double]|The intercept of the linear regression.|
|slope|Option[Double]|The slope of the linear regression.|
|r_squared|Option[Double]|The r-squared of the linear regression.|


### ErccDetailedMetric

Metrics produced by `CollectErccMetrics` describing various per-transcript metrics related to the spike-in of ERCC
(External RNA Controls Consortium) into an RNA-Seq experiment.  One metric per ERCC transcript will be present.


|Column|Type|Description|
|------|----|-----------|
|name|String|The name (or ID) of the ERCC transcript.|
|concentration|Double|The expected concentration as input to `CollectErccMetrics`.|
|count|Long|The observed count of the number of read pairs (or single end reads) .|
|normalized_count|Double|The observed count of the number of read pairs (or single end reads) normalized by the ERCC transcript length.|


### RunInfo

Stores the result of parsing the run info (RunInfo.xml) file from an Illumina run folder.


|Column|Type|Description|
|------|----|-----------|
|run_barcode|String|The unique identifier for the sequencing run and flowcell, stored as                   "<instrument-name>_<flowcell-barcode>".|
|flowcell_barcode|String|The flowcell barcode.|
|instrument_name|String|The instrument name.|
|run_date|Iso8601Date|The date of the sequencing run.|
|read_structure|ReadStructure|The description of the logical structure of cycles within the sequencing run.  This will only                       contain template and sample barcode segments, as the RunInfo.xml does not contain information                       about other segments (i.e. molecular barcodes and skips).|
|num_lanes|Int|The number of lanes in the flowcell.|


### PoolingFractionMetric

Metrics produced by `EstimatePoolingFractions` to quantify the estimated proportion of a sample
mixture that is attributable to a specific sample with a known set of genotypes.


|Column|Type|Description|
|------|----|-----------|
|sample|String|The name of the sample within the pool being reported on.|
|variant_sites|Count|How many sites were examined at which the reported sample is known to be variant.|
|singletons|Count|How many of the variant sites were sites at which only this sample was variant.|
|estimated_fraction|Proportion|The estimated fraction of the pool that comes from this sample.|
|standard_error|Double|The standard error of the estimated fraction.|
|ci99_low|Proportion|The lower bound of the 99% confidence interval for the estimated fraction.|
|ci99_high|Proportion|The upper bound of the 99% confidence interval for the estimated fraction.|


### ErrorRateByReadPositionMetric

Metrics produced by `ErrorRateByReadPosition` describing the number of base observations and
substitution errors at each position within each sequencing read.  Error rates are given for
the overall substitution error rate and also for each kind of substitution separately. Instead
of reporting 12 substitution rates, 6 are reported where complementary substitutions are grouped
together, e.g. `T>G` substitutions are reported as `A>C`.


|Column|Type|Description|
|------|----|-----------|
|read_number|Int|The read number (0 for fragments, 1 for first of pair, 2 for second of pair).|
|position|Int|The position or cycle within the read (1-based).|
|bases_total|Count|The total number of bases observed at this position.|
|errors|Count|The total number of errors or non-reference basecalls observed at this position.|
|error_rate|Double|The overall error rate at position.|
|a_to_c_error_rate|Double|The rate of `A>C` (and `T>G`) errors at the position.|
|a_to_g_error_rate|Double|The rate of `A>G` (and `T>C`) errors at the position.|
|a_to_t_error_rate|Double|The rate of `A>T` (and `T>A`) errors at the position.|
|c_to_a_error_rate|Double|The rate of `C>A` (and `G>T`) errors at the position.|
|c_to_g_error_rate|Double|The rate of `C>G` (and `G>C`) errors at the position.|
|c_to_t_error_rate|Double|The rate of `C>T` (and `G>A`) errors at the position.|


### ClippingMetrics

Metrics produced by ClipBam that detail how many reads and bases are clipped respectively.


|Column|Type|Description|
|------|----|-----------|
|read_type|ReadType|The type of read (i.e. Fragment, ReadOne, ReadTwo).|
|reads|Long|The number of reads examined.|
|reads_unmapped|Long|The number of reads that became unmapped due to clipping.|
|reads_clipped_pre|Long|The number of reads with any type of clipping prior to clipping with ClipBam.|
|reads_clipped_post|Long|The number of reads with any type of clipping after clipping with ClipBam, including reads that became unmapped.|
|reads_clipped_five_prime|Long|The number of reads with the 5' end clipped.|
|reads_clipped_three_prime|Long|The number of reads with the 3' end clipped.|
|reads_clipped_overlapping|Long|The number of reads clipped due to overlapping reads.|
|bases|Long|The number of aligned bases after clipping.|
|bases_clipped_pre|Long|The number of bases clipped prior to clipping with ClipBam.|
|bases_clipped_post|Long|The number of bases clipped after clipping with ClipBam, including bases from reads that became unmapped.|
|bases_clipped_five_prime|Long|The number of bases clipped on the 5' end of the read.|
|bases_clipped_three_prime|Long|The number of bases clipped on the 3 end of the read.|
|bases_clipped_overlapping|Long|The number of bases clipped due to overlapping reads.|
