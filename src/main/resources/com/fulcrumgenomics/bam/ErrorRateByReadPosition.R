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

# R script to generate QC plots from the ErrorRateByReadPosition tool

options(warn = -1) # Don't emit warnings, only errors
library(ggplot2)

args    = commandArgs(trailingOnly=T)
metrics = args[1]
output  = args[2]
name    = args[3]

data = read.table(metrics, header=T, sep="\t")

pdf(output, width=11, height=8.5)

ggplot(data) + aes(x=position) +
  geom_line(aes(y=error_rate)) +
  facet_wrap("read_number") +
  scale_y_sqrt() + 
  labs(x="Position in Read", y="Error Rate", title=paste(name, "-", "Error Rate by Position In Read")) +
  theme(plot.title = element_text(hjust = 0.5))

# The error rate plot is conditional on the presence and value of the collapsed column, which is not present in previous
# versions of fgbio
error_rate_plot = ggplot(data) + aes(x=position) +
  geom_line(aes(y=a_to_c_error_rate, color="A>C")) +
  geom_line(aes(y=a_to_g_error_rate, color="A>G")) +
  geom_line(aes(y=a_to_t_error_rate, color="A>T")) +
  geom_line(aes(y=c_to_a_error_rate, color="C>A")) +
  geom_line(aes(y=c_to_g_error_rate, color="C>G")) +
  geom_line(aes(y=c_to_t_error_rate, color="C>T"))
if ("collapsed" %in% names(data) && data$collapsed[[1]] == "false") {
  error_rate_plot = error_rate_plot +
    geom_line(aes(y=g_to_a_error_rate, color="G>A")) +
    geom_line(aes(y=g_to_c_error_rate, color="G>C")) +
    geom_line(aes(y=g_to_t_error_rate, color="G>T")) +
    geom_line(aes(y=t_to_a_error_rate, color="T>A")) +
    geom_line(aes(y=t_to_c_error_rate, color="T>C")) +
    geom_line(aes(y=t_to_g_error_rate, color="T>G"))
}
error_rate_plot = error_rate_plot +
  facet_wrap("read_number") +
  scale_y_sqrt() + 
  labs(x="Position in Read", y="Error Rate", title=paste(name, "-", "Error Rate by Position in Read and Type")) +
  theme(plot.title = element_text(hjust = 0.5))
print(error_rate_plot)

# The cumulative plot is more awkward because cumsum will do the wrong thing if we use ggplots aes(color=read_number),
# the R2 numbers will include the total at the end of R1!!  So instead we plot each read number separately!
cumulative_plot = ggplot(data) + aes(x=position, color=factor(read_number)) +
  labs(x="Position in Read", y="Cumulative Number of Errors", title=paste(name, "-", "Cumulative Error Count by Position in Read")) +
  guides(color=guide_legend(title="Read")) +
  theme(plot.title = element_text(hjust = 0.5))

for (readnum in unique(data$read_number)) {
  subdata = subset(data, data$read_number == readnum, select=c("read_number", "position", "bases_total", "error_rate"))
  subdata = rbind(c(readnum, 0, 0, 0), subdata) # adds in a point at 0,0
  cumulative_plot = cumulative_plot + geom_step(aes(y=cumsum(bases_total * error_rate)), data=subdata)
}

print(cumulative_plot)

dev.off()
