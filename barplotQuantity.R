#!/usr/bin/env Rscript


args <- commandArgs(TRUE)
file = args[1]


data = read.table(file)
data=data[order(data[,2],decreasing = T),]

large = dim(data)[1]*20

if (large < 600){
	large = 600
}

png(filename=paste(file,".png",sep = ""), width=as.integer(large), 600)
barplot(data[,2], names.arg=data[,1], main="", xlab="", ylab="Quantity", axes=TRUE, cex.axis=1, cex.lab=1, cex.main=1.5, cex.names = 1.2, col = "blue", las=2, space = 0.7, cex = 2)

dev.off()

