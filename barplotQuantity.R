#!/usr/bin/env Rscript


args <- commandArgs(TRUE)
file = args[1]


data = read.table(file)
data=data[order(data[,2],decreasing = T),]

large = dim(data)[1]*20

if (large < 300){
	large = 300
}

png(filename=paste(file,".png",sep = ""), width=as.integer(large), 300)
barplot(data[,2], names.arg=data[,1], main="Count Ligand", xlab="", ylab="Quantity", axes=TRUE, cex.axis=1, cex.lab=1, cex.main=1.5, cex.names = 0.6, col = "blue", las=2,space = 0.7, cex = 1)

dev.off()

