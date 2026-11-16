---
title: CopyUmiFromReadName
---

# CopyUmiFromReadName

## Overview
**Group:** Unique Molecular Identifiers (UMIs)

Copies the UMI at the end of the BAM's read name to the RX tag.

The read name is split on `:` characters with the last field assumed to be the UMI sequence.  The UMI
will be copied to the `RX` tag as per the SAM specification.  If any read does not have a UMI composed of
valid bases (ACGTN), the program will report the error and fail.

If a read name contains multiple UMIs they may be delimited (typically by a hyphen (`-`) or plus (`+`)).
The `--umi-delimiter` option specifies the delimiter on which to split.  The resulting UMI in the `RX` tag
will always be hyphen delimited.

Some tools (e.g. BCL Convert) may reverse-complement UMIs on R2 and add an 'r' prefix to indicate that the sequence
has been reverse-complemented.  By default, the 'r' prefix is removed and the sequence is reverse-complemented
back to the forward orientation.   The `--override-reverse-complement-umis` disables the latter behavior, such that
the 'r' prefix is removed but the UMI sequence is left as reverse-complemented.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|The input BAM file.|Required|1||
|output|o|PathToBam|The output BAM file.|Required|1||
|remove-umi||Boolean|Remove the UMI from the read name.|Optional|1|false|
|field-delimiter||Char|Delimiter between the read name and UMI.|Optional|1|:|
|umi-delimiter||Char|Delimiter between UMI sequences.|Optional|1|+|
|override-reverse-complement-umis||Boolean|Do not reverse-complement UMIs prefixed with 'r'.|Optional|1|false|

