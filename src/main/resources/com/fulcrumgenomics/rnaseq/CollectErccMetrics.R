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

# R script to generate QC plots from the CollectErccMetrics tool
# Three inputs are required:
# 1. The per-ERCC-transcript name/concentration/count
# 2. The output PDF to which to write.
# 3. The name of the sample/BAM.
# 4. The minimum # of counts to include a transcript

options(warn = -1) # Don't emit warnings, only errors
library(ggplot2)

args                 = commandArgs(trailingOnly=T)
metrics              = args[1]
output               = args[2]
name                 = args[3]
min_transcript_count = as.numeric(args[4])

data = read.table(metrics, header=T, sep="\t")

# Returns log2 of the value.  If the return value is infinite, then returns zero.  The supplied value should be greater
# than or equal to zero.
log2_positive <- function(x) {
  stopifnot(x >= 0)
  log2x = log2(x)
  if (is.infinite(log2x)) {
    return (log2(1e-10))
  }
  else {
    return (log2x)
  }
}

# enforce a minimum # of counts per-transcript
data = subset(data, data$count >= min_transcript_count)

# scale the data into log2 and perform a simple linear regression to get fit parameters.
log2data = data.frame(
  concentration = sapply(data$concentration, log2_positive),
  count = sapply(data$count, log2_positive),
  normalized_count = sapply(data$normalized_count, log2_positive)
)

pdf(output, width=11, height=8.5)

fit = lm(concentration ~ count, log2data)
ggplot(data) + aes(x=concentration, y=count) +
  geom_point(shape=1) +
  geom_smooth(method=lm, col="red") +
  scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='log2') +
  labs(title=paste("ERCC Correlation (raw counts)", name),
    x="Expected concentration (log2)", y="Observed counts (log2)",
    subtitle = paste("R2 = ", signif(summary(fit)$adj.r.squared, 5),
      " Intercept =", signif(fit$coef[[1]],5 ),
      " Slope =", signif(fit$coef[[2]], 5))) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.title=element_blank())

fit = lm(concentration ~ normalized_count, log2data)
ggplot(data) + aes(x=concentration, y=normalized_count) +
  geom_point(shape=1) +
  geom_smooth(method=lm, col="red") +
  scale_x_continuous(trans='log2') +
  scale_y_continuous(trans='log2') +
  labs(title=paste("ERCC Correlation (normalized counts)", name),
  x="Expected concentration (log2)", y="Observed counts normalized by transcript length (log2)",
  subtitle = paste("R2 = ", signif(summary(fit)$adj.r.squared, 5),
  " Intercept =", signif(fit$coef[[1]],5 ),
  " Slope =", signif(fit$coef[[2]], 5))) +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), legend.title=element_blank())

dev.off()
