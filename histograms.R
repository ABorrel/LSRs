#!/usr/bin/env Rscript

# BORREL Alexandre
# 04-2012

args <- commandArgs(TRUE)

file = args[1]
brk = as.integer (args[2])


# histograms

# open both case with header and without header
d = read.table (file, header = TRUE, sep = "\t")


# cut function number col
nb_hist = dim (d)[2] - 1

png (paste (file, ".png", sep = ""), 400, 400*nb_hist)
par (mfrow = c(nb_hist, 1))

for (i in seq (1, nb_hist)){
	hist (d[,i+1], xlim = c(min(d[,i+1]), max(d[,i+1])), breaks = brk, main = colnames (d)[i+1], col = "grey")
}
dev.off()



