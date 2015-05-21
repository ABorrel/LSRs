#!/usr/bin/env Rscript

# BORREL Alexandre
# 04-2012

args <- commandArgs(TRUE)

file = args[1]
brk = as.integer (args[2])
max_x = as.double(args[3])

# histograms

# open both case with header and without header
d = read.table (file, header = TRUE, sep = "\t")


# cut function number col
nb_hist = dim (d)[2] - 1

pdf (paste (file, ".pdf", sep = ""))
for (i in seq (1, nb_hist)){
	hist (d[,i+1], xlim = c(min(d[,i+1]), max_x), breaks = brk, main = colnames (d)[i+1], col = "grey")
}
dev.off()

for (i in seq (1, nb_hist)){
	
        if (colnames (d)[i+1]=="D_max"){
        	name_xlab = "Distance (Å)"

        	png (paste (file, colnames (d)[i+1], ".png", sep = ""), 600, 600)
		par(mar=c(5,5,5,2))
	        par (oma = c(1,1,1,1))
		hist (d[,i+1], breaks = brk, ylab = "Number of occurencies", xlab = name_xlab, col = "grey", cex.lab = 1.6, cex.axis = 1.2, main = "")
		dev.off()

		svg (paste (file, colnames (d)[i+1], ".svg", sep = ""), 6, 6)
		par(mar=c(4,4.2,1,1))
	        par (oma = c(0.5,0.5,0.5,0.5))	
		hist (d[,i+1], breaks = brk, ylab = "Number of occurencies", xlab = name_xlab, col = "grey", cex.lab = 1.6, cex.axis = 1.2, main = "")
		dev.off()

		png (paste (file, "identic_", colnames (d)[i+1], ".png", sep = ""), 600, 600)
		par(mar=c(5,5,5,2))
		par (oma = c(1,1,1,1))
	        i_identic = which (d[,"identic"]==1)
		hist (d[i_identic,i+1], breaks = brk, ylab = "Number of occurencies", xlab = name_xlab, col = "grey", cex.lab = 1.6, cex.axis = 1.2, main = "")
		dev.off()

		svg (paste (file, "identic_", colnames (d)[i+1], ".svg", sep = ""), 6, 6)
		par(mar=c(4,4.2,1,1))
		par (oma = c(0.5,0.5,0.5,0.5))
       		i_identic = which (d[,"identic"]==1)
		hist (d[i_identic,i+1], breaks = brk, ylab = "Number of occurencies", xlab = name_xlab, col = "grey", cex.lab = 1.6, cex.axis = 1.2, main = "")
		dev.off()

	}else{
        	name_xlab = "RMSD (Å)"


		png (paste (file, colnames (d)[i+1], ".png", sep = ""), 600, 600)
		par(mar=c(5,5,5,2))
	        par (oma = c(1,1,1,1))
		hist (d[,i+1], xlim = c(0, max_x), breaks = brk, ylab = "Number of occurencies", xlab = name_xlab, col = "grey", cex.lab = 1.6, cex.axis = 1.2, main = "")
		dev.off()
	
		svg (paste (file, colnames (d)[i+1], ".svg", sep = ""), 6, 6)
		par(mar=c(4,4.2,1,1))
	   	par (oma = c(0.5,0.5,0.5,0.5))
		hist (d[,i+1], xlim = c(0, max_x), breaks = brk, ylab = "Number of occurencies", xlab = name_xlab, col = "grey", cex.lab = 1.6, cex.axis = 1.2, main = "")
		dev.off()

		png (paste (file, "identic_", colnames (d)[i+1], ".png", sep = ""), 600, 600)
		par(mar=c(5,5,5,2))
	        par (oma = c(1,1,1,1))
	        i_identic = which (d[,"identic"]==1)
		hist (d[i_identic,i+1], xlim = c(0, max_x), breaks = brk, ylab = "Number of occurencies", xlab = name_xlab, col = "grey", cex.lab = 1.6, cex.axis = 1.2, main = "")
		dev.off()

		svg (paste (file, "identic_", colnames (d)[i+1], ".svg", sep = ""), 6, 6)
		par(mar=c(4,4.2,1,1))
	        par (oma = c(0.5,0.5,0.5,0.5))
	        i_identic = which (d[,"identic"]==1)
		hist (d[i_identic,i+1], xlim = c(0, max_x), breaks = brk, ylab = "Number of occurencies", xlab = name_xlab, col = "grey", cex.lab = 1.6, cex.axis = 1.2, main = "")
		dev.off()

	}
}

