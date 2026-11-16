---
title: CopyUmiFromReadName
---

# CopyUmiFromReadName

## Overview
**Group:** Unique Molecular Identifiers (UMIs)

Copies the UMI at the end of the BAM's read name to the RX tag.

The read name is split on `:` characters with the last field is assumed to be the UMI sequence.  The UMI
will be copied to the `RX` tag as per the SAM specification.  If any read does not have a UMI composed of
valid bases (ACGTN), the program will report the error and fail.

If a read name contains multiple UMIs they may be delimited by either hyphens (`-`) or pluses (`+`). The
resulting UMI in the `RX` tag will always be hyphen delimited.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|The input BAM file|Required|1||
|output|o|PathToBam|The output BAM file|Required|1||
|remove-umi||Boolean|Remove the UMI from the read name|Optional|1|false|

