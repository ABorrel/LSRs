#!/usr/bin/env Rscript

# BORREL Alexandre
# 04-2012

args <- commandArgs(TRUE)

file = args[1]
type = args[2]
brk = as.integer (args[3])

# histograms

# open both case with header and without header
d = read.table (file, header = FALSE, sep = "\t")
print (dim (d))

rownames (d) = d[,1]
d = d[,-1]

png (paste (file, ".png", sep = ""), 800, 800)
par (mar = c(5,5,5,5))
hist (d, xlim = c(min(d), max(d)), breaks = brk, main = "", col = "grey", xlab = "Distance (Å)", ylab = "Number of occurences", cex.lab = 2, cex.axis = 2)
dev.off()

svg (paste (file, ".svg", sep = ""), 22, 22)
par(mar=c(8,8,8,8))
hist (d, xlim = c(min(d), max(d)), breaks = brk, main = "", col = "grey", xlab = "Distance (Å)", ylab = "Number of occurences", cex.lab = 2.5, cex.axis = 2)
dev.off()
