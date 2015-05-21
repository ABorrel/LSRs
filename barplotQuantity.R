#!/usr/bin/env Rscript


args <- commandArgs(TRUE)
file = args[1]


data = read.table(file, header = TRUE, sep = "\t")
data = data[order(data[,2],decreasing = T),]

large = dim(data)[1]*20

if (large < 600){
	large = 600
}

png(filename=paste(file,".png",sep = ""), width=as.integer(large), 600)
par(mar = c(5,5,5,2))
par(oma = c(1,1,1,1))
barplot(data[,2], names.arg=data[,1], main="", xlab="", ylab="Number of occurencies", axes=TRUE, cex.axis = 1.2, cex.lab=1.6, cex.main=1.5, cex.names = 1.6, col = "grey", las=2, space = 0.8, cex = 1.6)
dev.off()

svg(filename=paste(file,".svg",sep = ""), width=(dim(data)[1]), height = (dim(data)[1]))
par(mar = c(2,6,1,1))
par(oma = c(1,1,1,1))
barplot(data[,2], names.arg=data[,1], main="", xlab="", ylab="Number of occurencies", axes=TRUE, cex.axis = 1.2, cex.lab = 1.6, cex.main = 1.5, cex.names = 1.6, col = "grey", las = 2, space = 0.8, cex = 1.6)
dev.off()

