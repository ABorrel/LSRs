#!/usr/bin/env Rscript

# BORREL Alexandre
# 04-2012

args <- commandArgs(TRUE)

file = args[1]
type = args[2]
brk = as.integer (args[3])


# histograms

# open both case with header and without header
d = read.table (file, header = TRUE)

# cut function number col
nb_hist = dim (d)[2]


png (paste (file, ".png", sep = ""), 400, 400*nb_hist)
par (mfrow = c(nb_hist, 1))

for (i in seq (1, nb_hist)){
	hist (d[,i], xlim = c(min(d[,i]), max(d[,i])), breaks = brk, main = colnames (d)[i], col = "grey")
}
dev.off()



