---
title: FixVcfPhaseSet
---

# FixVcfPhaseSet

## Overview
**Group:** VCF/BCF

Adds/fixes the phase set (PS) genotype field.

The VCF specification allows phased genotypes to be annotated with the `PS` (phase set) `FORMAT` field.  The value
should be a non-negative integer, corresponding to the position of the first variant in the phase set.  Some tools
will output a non-integer value, as well as describe this field as having non-Integer type in the VCF header.  This
tool will update the phase set (`PS`) `FORMAT` field to be VCF spec-compliant.

The handling of unphased genotypes with phase-sets is controlled by the `-x` option:
- If `-x` is used the genotype will be converted to a phased genotype (e.g. `0/1` => `0|1`)
- Otherwise the phase-set (`PS`) will be removed from the genotype

The `--keep-original` option may be used to store the original `PS` value in a new `OPS` field.  The type
described in the header will match the original.

This tool cannot fix phased variants without a phase set, or phased variant sets who have different phase set
values.

In some cases, VCFs (e.g. from GIAB/NIST or Platinum Genomes) have illegal header lines, for example, a `PEDIGREE`
header line without a `ID` key-value field.  The `-z` option can be used to remove those lines.  This option is
included in this tool for convenience as those example VCFs in some cases have these illegal header lines, and it
is convenient to fix the phase set in addition to removing those illegal header lines.

## Arguments

|Name|Flag|Type|Description|Required?|Max # of Values|Default Value(s)|
|----|----|----|-----------|---------|---------------|----------------|
|input|i|PathToVcf|Input VCF.|Required|1||
|output|o|PathToVcf|Output VCF.|Required|1||
|keep-original|k|Boolean|Store the original phase set in the `OPS` field.|Optional|1|false|
|phase-genotypes-with-phase-set|x|Boolean|Set unphased genotypes with a PS FORMAT value to be phased.|Optional|1|false|
|remove-no-id-header-lines|z|Boolean|Remove header lines that do not contain an ID key-value, which is required in VCF.|Optional|1|false|

