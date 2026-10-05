---
title: ExtractBasecallingParamsForPicard
---

# ExtractBasecallingParamsForPicard

## Overview
**Group:** Basecalling

Extracts sample and library information from an sample sheet for a given lane.

The sample sheet should be an Illumina Experiment Manager sample sheet. The tool writes two files to the output
directory: a barcode parameter file and a library parameter file.

The barcode parameter file is used by Picard's `ExtractIlluminaBarcodes` and `CollectIlluminaBasecallingMetrics` to
determine how to match sample barcodes to each read.  The parameter file will be written to the output directory
with name `barcode_params.<lane>.txt`.

The library parameter file is used by Picard's `IlluminaBasecallsToSam` to demultiplex samples and name the output
BAM file path for each sample output BAM file.  The parameter file will be written to the output directory with name
`library_params.<lane>.txt`.  The path to each sample's BAM file will be specified in the library parameter
file.  Each BAM file will have path `<output>/<sample-name>.<barcode-sequence>.<lane>.bam`.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|FilePath|The input sample sheet.|Required|1||
|output|o|DirPath|The output folder to where per-lane parameter files should be written.|Required|1||
|bam-output|b|DirPath|Optional output folder to where per-lane BAM files should be written, otherwise the output directory will be used.|Optional|1||
|lanes|l|Int|The lane(s) (1-based) for which to write per-lane parameter files.|Required|Unlimited||

