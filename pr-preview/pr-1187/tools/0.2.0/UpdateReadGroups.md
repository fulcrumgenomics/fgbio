---
title: UpdateReadGroups
---

# UpdateReadGroups

## Overview
**Group:** SAM/BAM

Updates one or more read groups and their identifiers.

This tool will replace each read group with a new read group, including a new read group identifier.  If the read
group identifier is not to be changed, it is recommended to use `samtools reheader` or Picard's
`ReplaceSamHeader` instead as in this case only the header needs modification. If all read groups are to be
assigned to one read group, it is recommended to use Picard's `AddOrReplaceReadGroups`.  Nonetheless, if the read
group identifier also needs to be changed, use this tool.

Each read group in the input file will be mapped to one and only one new read group identifier, unless
`--ignore-missing-read-groups` is set.  A SAM header file should be given with the new read groups and the ID field
foreach read group containing the new read group identifier.  An additional attribute (`FR`) should be provided
that gives the original read group identifier (`ID`) to which this new read group corresponds.

If `--keep-read-group-attributes` is true, then any read group attribute not replaced will be kept in the new read
group.  Otherwise, only the attributes in the provided SAM header file will be used.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input BAM file.|Required|1||
|output|o|PathToBam|Output BAM file.|Required|1||
|read-groups-file|r|PathToBam|A SAM header file with the replacement read groups (see detailed usage).|Required|1||
|keep-read-group-attributes|k|Boolean|Keep all read group attributes that are not replaced.|Optional|1|false|
|ignore-missing-read-groups|g|Boolean|Keep all read groups not found in the replacement header, otherwise throw an error.|Optional|1|false|

