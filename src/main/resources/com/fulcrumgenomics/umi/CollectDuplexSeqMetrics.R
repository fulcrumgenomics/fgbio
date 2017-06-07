#################################################################################
# The MIT License                                                               #
# Copyright (c) 2017 Fulcrum Genomics LLC                                       #
# Permission is hereby granted, free of charge, to any person obtaining a copy  #
# of this software and associated documentation files (the "Software"), to deal #
# in the Software without restriction, including without limitation the rights  #
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     #
# copies of the Software, and to permit persons to whom the Software is         #
# furnished to do so, subject to the following conditions:                      #
# The above copyright notice and this permission notice shall be included in    #
# all copies or substantial portions of the Software.                           #
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, #
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN     #
# THE SOFTWARE.                                                                 #
#################################################################################

#################################################################################
# R script to generate QC plots for duplex sequencing data QC's with
# CollectDuplexSeqMetrics tools.  All arguments are positional and required:
#   1.  The family size file with data on families by start/stop, ss and ds
#   2.  The duplex family size file
#   3.  The duplex yield metrics file
#   4.  The UMI metrics file
#   5.  An output PDF to write
#   6+. One or more strings to use in plot titles to ID the sample/dataset
#################################################################################

options(warn = -1) # Don't emit warnings, only errors
library(ggplot2)
library(scales)  # Require by ggplot

# Standard colors
fgblue  = "#2E9DD7"
fggreen = "#155936"
fgred   = "firebrick3"
fgcolors = c(fgblue, fggreen, fgred)

args             = commandArgs(trailingOnly=T)
familyData       = read.table(args[1], sep="\t", header=T)
duplexFamilyData = read.table(args[2], sep="\t", header=T)
yieldData        = read.table(args[3], sep="\t", header=T)
umiData          = read.table(args[4], sep="\t", header=T)
outputFile       = args[5]
sampleInfo       = paste(args[6:length(args)])

pdf(outputFile, width=11, height=8.5)

# Plot #1 - Family Sizes for different kinds of familes
ggplot(familyData) + aes(x=family_size) +
  geom_line(aes(y=ds_count, color="DS Families")) +
  geom_line(aes(y=ss_count, color="SS Families")) + 
  geom_line(aes(y=cs_count, color="By Coord+Strand")) +
  scale_x_continuous(trans="log2", minor_breaks=seq(0,max(familyData$family_size), by=2)) +
  scale_color_manual(values=fgcolors) +
  labs(x="Family Size (log2 scaled)", y="Count of Families", title=paste("Family Size Distributions for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #2- Cumulative Family Sizes for different kinds of familes
ggplot(familyData) + aes(x=family_size) +
  geom_line(aes(y=cs_fraction_gt_or_eq_size, color="By Coord+Strand"), alpha=0.5) +
  geom_line(aes(y=ss_fraction_gt_or_eq_size, color="SS Families"), alpha=0.5) + 
  geom_line(aes(y=ds_fraction_gt_or_eq_size, color="DS Families"), alpha=0.5) +
  scale_x_continuous(trans="log2", minor_breaks=seq(0,max(familyData$family_size), by=2)) +
  scale_color_manual(values=fgcolors) +
  labs(x="Family Size (log2 scaled)", y="Fraction of Families at >= Family Size", title=paste("Cumulative Family Size Distributions for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #3a - Duplex Family Size Distribution
ggplot(duplexFamilyData) +
  aes(x=ab_size, y=ba_size, z=count) +
  stat_summary_2d(fun=sum, binwidth=1) +
  scale_fill_gradient2(low=fgblue, mid="white", high=fggreen, space="Lab",
                       midpoint=log10(max(duplexFamilyData$count)), guide = "colourbar", trans="log2") +
  labs(x="AB Reads", y="BA Reads", title=paste("Duplex Tag Family Size Distribution for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5))

# Plot #3b - Duplex Family Size Distribution for families with ab > 1 and ba > 1
ggplot(subset(duplexFamilyData, duplexFamilyData$ba_size > 0)) +
  aes(x=ab_size, y=ba_size, z=count) +
  stat_summary_2d(fun=sum, binwidth=1) +
  scale_fill_gradient2(low=fgblue, mid="white", high=fggreen, space="Lab",
                       midpoint=log10(max(duplexFamilyData$count)), guide = "colourbar", trans="log2") +
  labs(x="AB Reads", y="BA Reads", title=paste("Duplex Tag Family Size Distribution for ", sampleInfo, " ; Only Families with AB>0 & BA>0", sep="")) +
  theme(plot.title = element_text(hjust = 0.5))

# Plot #4 - Duplex Yield
ggplot(yieldData) +
  aes(x=read_pairs) +
  geom_area(aes(y=ds_families, fill="DS Families")) +
  geom_area(aes(y=ds_families * ds_fraction_duplexes_ideal, fill="Duplexes - Ideal")) +
  geom_area(aes(y=ds_duplexes, fill="Duplexes - Actual")) +
  scale_fill_manual(values=fgcolors) +
  labs(x="Read Pairs", y="Count of Double-Strand Families and Duplexes", title=paste("Duplex Yield by Input Read Pairs for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #5 UMI distribution (for UMIs that don't contain any no-calls)
ggplot(subset(umiData, !grepl("N", umiData$umi, fixed=T))) +
  aes(x=raw_observations, y=unique_observations) +
  geom_point(color=fggreen) +
  labs(x="Observations of UMI in Raw Reads", y="Unique Observations (Tag Families w/UMI)", title=paste("UMI Representation in", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

# Plot #6 - Distribution of reads within tag families
ggplot(familyData) + aes(x=family_size) +
  geom_line(aes(y=cs_count * family_size, color="By Coord+Strand")) +
  geom_line(aes(y=ss_count * family_size, color="SS Families")) + 
  geom_line(aes(y=ds_count * family_size, color="DS Families")) +
  scale_x_continuous(trans="log2", minor_breaks=seq(0,max(familyData$family_size), by=2)) +
  scale_color_manual(values=fgcolors) +
  labs(x="Family Size (log2 scaled)", y="Reads Allocated to Families of Size N", title=paste("Read Distribution Among Families for", sampleInfo)) +
  theme(plot.title = element_text(hjust = 0.5), legend.title=element_blank())

dev.off()
