#!/usr/bin/env Rscript

# BORREL Alexandre
# 04-2012

args <- commandArgs(TRUE)

file = args[1]
type = args[2]
brk = as.integer (args[3])


# histograms

# open both case with header and without header
d = read.table (file, header = FALSE)
rownames (d) = d[,1]
d = d[,-1]

png (paste (file, ".png", sep = ""), 400, 400)

hist (d, xlim = c(min(d), max(d)), breaks = brk, main = type, col = "grey")
dev.off()



