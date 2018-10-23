#################################################################################
# The MIT License                                                               #
# Copyright (c) 2018 Fulcrum Genomics LLC                                       #
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

# R script to generate QC plots from the FindSwitchbackReads tool

options(warn = -1) # Don't emit warnings, only errors
library(ggplot2)

args        = commandArgs(trailingOnly=T)
name        = args[1]
lengthsFile = args[2]
offsetsFile = args[3]
gapsFile    = args[4]
outputFile  = args[5]

pdf(outputFile, width=11, height=8.5)

# Standard colors
files  = c(lengthsFile, offsetsFile, gapsFile)
labels = c("Length of Switchback", "Switchback Offset", "Tandem Switchback Gap")
colors = c("#2E9DD7", "#155936", "firebrick3")

for (i in c(1,2,3)) {
  file  = files[i]
  color = colors[i]
  xlab  = labels[i]

  data = read.table(file, sep="\t", header=T)
  names(data)[1] = "key"
  plot = ggplot(data) + aes(x=key, y=count) + 
    scale_y_continuous(limits=c(0, NA)) +
    geom_bar(stat="identity", color=color, fill=color) +
    labs(x=xlab, y="Number of Templates", title=paste("Distribution of", xlab, "in", name)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  print(plot)
}

dev.off()
