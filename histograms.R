#!/usr/bin/env Rscript

# BORREL Alexandre
# 04-2012

args <- commandArgs(TRUE)

file = args[1]
type = args[2]
brk = as.integer (args[3])
organised = as.integer (args[4])


# histograms

d = read.table (file, header = FALSE)

# cut function number col
nb_hist = dim (d)[2] - 1


png (paste (file, ".png", sep = ""), 400, 400*nb_hist)
par (mfrow = c(nb_hist, 1))


for (i in seq (2, nb_hist+1)){
	hist (d[,i], xlim = c(min(d[,i]), max(d[,i])), breaks = brk, main = type, col = "grey")
}
dev.off()



