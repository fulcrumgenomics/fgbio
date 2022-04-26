################################################################################
# Example snakemake file that implements the R&D Best Practices for fgbio's
# Fastq -> Filtered Consensus Pipeline
#
# Required software:
#   snakemake (!)
#   bwa
#   samtools v1.14 or higher
#   fgbio v2.0.0 of higher
################################################################################

# Adjust these parameters to match your dataset
r1_fastq          = "r1.fq.gz"
r2_fastq          = "r2.fq.gz"
r1_read_structure = "8M+T"
r2_read_structure = "8M+T"
genome            = "hs38DH/hs38DH.fa"

# Since both rules can generate a uBam...
ruleorder: call_consensus_reads > fastq_to_ubam

rule all:
    input:
        "my_sample.mapped.bam",
        "my_sample.grouped.bam",
        "my_sample.cons.filtered.bam",


rule fastq_to_ubam:
    """Generates a uBam from R1 and R2 fastq files."""
    input:
        r1 = r1_fastq,
        r2 = r2_fastq
    params:
        rs1 = r1_read_structure,
        rs2 = r2_read_structure,
    output:
        bam = temp("{sample}.unmapped.bam")
    resources:
        mem_gb = 1
    log:
        "logs/fastq_to_ubam.{sample}.log"
    shell:
        " fgbio -Xmx1g --compression 1 --async-io FastqToBam "
        "   --input {input.r1} {input.r2} "
        "   --read-structures {params.rs1} {params.rs2} "
        "   --sample {wildcards.sample} "
        "   --library {wildcards.sample} "
        "   --platform-unit flowcell.lane "
        "   --output {output.bam} &> {log} "


rule align_bam:
    """Takes an unmapped BAM and generates an aligned BAM using bwa and ZipperBams."""
    input:
        bam = "{prefix}.unmapped.bam",
        ref = genome
    output:
        bam = "{prefix}.mapped.bam"
    threads:
        16
    resources:
        mem_gb = 14
    log:
        "logs/align_bam.{prefix}.log"
    shell:
        " ( "
        " samtools fastq {input.bam} "
        "   | bwa mem -t {threads} -p -K 150000000 -Y {input.ref} - "
        "   | fgbio -Xmx4g --compression 1 --async-io ZipperBams "
        "       --unmapped {input.bam} "
        "       --ref {input.ref} "
        "       --output {output.bam} "
        "       --tags-to-reverse Consensus "
        "       --tags-to-revcomp Consensus "
        " ) &> {log}"
        

rule group_reads:
    """Group the raw reads by UMI and position ready for consensus calling."""
    input:
        bam = "{sample}.mapped.bam",
    output:
        bam = "{sample}.grouped.bam",
        stats = "{sample}.grouped-family-sizes.txt"
    params:
        allowed_edits = 1,
    threads:
        2
    resources:
        mem_gb = 8
    log:
        "logs/group_reads.{sample}.log"
    shell:
        "fgbio -Xmx8g --compression 1 --async-io GroupReadsByUmi "
        "  --input {input.bam} "
        "  --strategy Adjacency "
        "  --edits {params.allowed_edits} "
        "  --output {output.bam} "
        "  --family-size-histogram {output.stats} "
        "  &> {log} "


rule call_consensus_reads:
    """Call consensus reads from the grouped reads."""
    input:
        bam = "{sample}.grouped.bam",
    output:
        bam = temp("{sample}.cons.unmapped.bam"),
    params:
        min_reads = 1,
        min_base_qual = 20
    threads:
        4
    resources:
        mem_gb = 8
    log:
        "logs/call_consensus_reads.{sample}.log"
    shell:
        "fgbio -Xmx4g --compression 1 CallMolecularConsensusReads "
        "  --input {input.bam} "
        "  --output {output.bam} "
        "  --min-reads {params.min_reads} "
        "  --min-input-base-quality {params.min_base_qual} "
        "  --threads {threads} "
        "  &> {log}"


rule filter_consensus_reads:
    """Filters the consensus reads and then sorts into coordinate order."""
    input:
        bam = "{sample}.cons.mapped.bam",
        ref = genome,
    output:
        bam = "{sample}.cons.filtered.bam",
    params:
        min_reads = 3,
        min_base_qual = 40,
        max_error_rate = 0.2
    threads:
        8
    resources:
        mem_gb = 8
    log:
        "logs/filter_consensus_reads.{sample}.log"
    shell:
        " ( "
        " fgbio -Xmx8g --compression 0 FilterConsensusReads "
        "   --input {input.bam} "
        "   --output /dev/stdout "
        "   --ref {input.ref} "
        "   --min-reads {params.min_reads} "
        "   --min-base-quality {params.min_base_qual} "
        "   --max-base-error-rate {params.max_error_rate} "
        "   | samtools sort --threads {threads} -o {output.bam}##idx##{output.bam}.bai --write-index "
        " ) &> {log} "
