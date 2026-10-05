---
title: RemoveSamTags
---

# RemoveSamTags

## Overview
**Group:** SAM/BAM

Removes SAM tags from a SAM or BAM file.  If no tags to remove are given, the original file is produced.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|PathToBam|Input SAM or BAM.|Required|1||
|output|o|PathToBam|Output SAM or BAM.|Required|1||
|tags-to-remove|t|String|The tags to remove.|Optional|Unlimited||

