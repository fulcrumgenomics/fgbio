---
title: AutoGenerateReadGroupsByName
---

# AutoGenerateReadGroupsByName

## Overview
**Group:** SAM/BAM

Adds read groups to a BAM file for a single sample by parsing the read names.

Will add one or more read groups by parsing the read names.  The read names should be of the form:

```
<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<xpos>:<y-pos>
```

Each unique combination of `<instrument>:<run number>:<flowcell ID>:<lane>` will be its own read group. The ID of the
read group will be an integer and the platform unit will be `<flowcell-id>.<lane>`.

The input is assumed to contain reads for one sample and library.  Therefore, the sample and library must be given
and will be applied to all read groups.  Read groups will be replaced if present.

Two passes will be performed on the input: first to gather all the read groups, and second to write the output BAM
file.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToBam|Input SAM or BAM file|Required|1||
|output|o|PathToBam|Output SAM or BAM file|Required|1||
|sample|s|String|The sample to insert into the read groups|Required|1||
|library|l|String|The library to insert into the read groups|Required|1||
|sequencing-center||String|The sequencing center from which the data originated|Optional|1||
|predicted-insert-size||Integer|Predicted median insert size, to insert into the read groups|Optional|1||
|program-group||String|Program group to insert into the read groups|Optional|1||
|platform-model||String|Platform model to insert into the groups (free-form text providing further details of the platform/technology used)|Optional|1||
|description||String|Description inserted into the read groups|Optional|1||
|run-date||Iso8601Date|Date the run was produced (ISO 8601: `YYYY-MM-DD` ), to insert into the read groups|Optional|1||
|comments||String|Comment(s) to include in the merged output file's header.|Optional|Unlimited||
|sort-order||SamOrder|The sort order for the output sam/bam file.|Optional|1||

