#!/usr/bin/env Rscript




pieType = function (d, path_out){
	
	# print (d)
	colors = seq (1,dim(d)[2])

	leg = NULL
	for (l in names (d)){
		leg = append (leg, paste (l, "\n", d[l], sep = ""))
	}

	par (lwd = 1000)
	png(filename=paste(path_out,".png",sep = ""),400, 400)
	try(pie(as.double(d), col = colors, label = names(d), lwd = 10))
	dev.off()

	svg(filename=paste(path_out,".svg",sep = ""))
	try(pie(as.double(d), col = colors, label = leg))
	dev.off()
}



# MAIN #
args <- commandArgs(TRUE)
p_filin = args[1]


d = read.table (p_filin, header = TRUE)

pieType  (d, p_filin)
