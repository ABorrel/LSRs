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

png (paste (file, ".png", sep = ""), 600, 600)
par (mar = c(5,5,5,5))
par (oma = c(1,1,1,1))
hist (d, xlim = c(2,6), breaks = brk, main = "", col = "grey", xlab = "Distance (Å)", ylab = "Number of occurences", cex.lab = 1.6, cex.axis = 1.2)
dev.off()

svg (paste (file, ".svg", sep = ""), 6, 6)
par(mar=c(4,4.2,1,1))
par (oma = c(0.5,0.5,0.5,0.5))
hist (d, xlim = c(2, 6), breaks = brk, main = "", col = "grey", xlab = "Distance (Å)", ylab = "Number of occurences", cex.lab = 1.6, cex.axis = 1.2)
dev.off()
