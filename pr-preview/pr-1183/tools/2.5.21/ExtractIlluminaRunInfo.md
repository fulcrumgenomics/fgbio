---
title: ExtractIlluminaRunInfo
---

# ExtractIlluminaRunInfo

## Overview
**Group:** Basecalling

Extracts information about an Illumina sequencing run from the RunInfo.xml.

The output file will contain a header column and a single column containing the following rows:

1. `run_barcode:` the unique identifier for the sequencing run and flowcell, stored as `<instrument-name>_<flowcell-barcode>`.
2. `flowcell_barcode:` the flowcell barcode.
3. `instrument_name`: the instrument name.
4. `run_date`: the date of the sequencing run.
5. `read_structure`: the description of the logical structure of cycles within the sequencing run, including which cycles
   correspond to sample barcodes, molecular barcodes, template bases, and bases that should be skipped.
6. `number_of_lanes`: the number of lanes in the flowcell.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|FilePath|The input RunInfo.xml typically found in the run folder.|Required|1||
|output|o|FilePath|The output file.|Required|1||

