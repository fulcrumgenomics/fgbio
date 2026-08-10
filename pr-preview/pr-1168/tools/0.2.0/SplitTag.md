---
title: SplitTag
---

# SplitTag

## Overview
**Group:** Personal

Splits an optional tag in a SAM or BAM into multiple optional tags.

## Arguments

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input||PathToBam|Input SAM or BAM.|Required|1||
|output||PathToBam|Output SAM or BAM.|Required|1||
|tag-to-split||String|Tag to split.|Required|1||
|tags-to-output||String|Tag(s) to output.  There should be one per produced token.|Required|Unlimited||
|delimiter||String|The delimiter used to split the string.|Optional|1|-|

