#!/usr/bin/env Rscript



vectCol = function(vin){
  
  vout = c()
  for (classLSR in vin){
    if (classLSR == "P" || classLSR == "cycle-P"){
      vout = append(vout, "#f78e04")
    }
    else if (classLSR == "B" || classLSR == "cycle-B"){
      vout = append(vout, "#c1fff5")
    }
    else if (classLSR == "F" || classLSR == "cycle-F"){
      vout = append(vout, "#faff7f")
    }
    else if (classLSR == "Cl" || classLSR == "cycle-Cl"){
      vout = append(vout, "#00680a")
    }
    else if (classLSR == "Br" || classLSR == "cycle-Br"){
      vout = append(vout, "#ef8e0e")
    }
    else if (classLSR == "Be" || classLSR == "cycle-Be"){
      vout = append(vout, "#ff19fb")
    }
    else if (classLSR == "NO2" || classLSR == "cycle-NO2"){
      vout = append(vout, "#6d0305")
    }
    else if (classLSR == "SO2" || classLSR == "cycle-SO2"){
      vout = append(vout, "#847e02")
    }
    else if (classLSR == "S" || classLSR == "cycle-S"){
      vout = append(vout, "#fff200")
    }
    else if (classLSR == "_CON" || classLSR == "cycle-CON"){
      vout = append(vout, "#00ff37")
    }
    else if (classLSR == "COO" || classLSR == "cycle-COO"){
      vout = append(vout, "#ff0008")
    }
    else if (classLSR == "onlyC" || classLSR == "cycle-onlyC"){
      vout = append(vout, "#383a38")
    }
    else if (classLSR == "C+O" || classLSR == "cycle-C+O"){
      vout = append(vout, "#965b4c")
    }
    else if (classLSR == "C+N" || classLSR == "cycle-C+N"){
      vout = append(vout, "#040d89")
    }
    else if (classLSR == "C+O+N" || classLSR == "cycle-C+O+N"){
      vout = append(vout, "#c489b9")
    }
    else if (classLSR == "other" || classLSR == "cycle-other"){
      vout = append(vout, "#afafa8")
    }
    else{
      vout = append(vout, "#ffffff")
    }
  }  
  
  names(vout) = vin
  
  vunique = unique(vout)
  names (vunique) = unique(names(vout))

  return(list(vout,vunique))
}


correlationPlot = function(din, pname){

  #print (din)
  vcol = vectCol(din[,1])
    
  svg(paste(pname, "_ESPVSshape.svg", sep = ""), 10, 10, bg = "transparent")
  par(mar=c(5,5,2,2))
  plot(din[,2], din[,3], ylab = "Shape", xlab = "ESP", cex.lab = 2, cex.axis = 1.5, col = vcol[[1]], pch = 20, cex = 2, ylim = c(0.45, 1), xlim = c(0.4, 1))
  legend("bottomright", legend = (names(vcol[[2]])), pch = 20, col = vcol[[2]], cex = 1.5)
  dev.off 
  
}



boxplotFormat = function(din, pname){

  # control order in boxplot
  din$ClassLSR<-factor(din$ClassLSR, levels=unique (din[,1]))
  
  #color
  vcol = vectCol(din[,1])[[2]]
  
  # Boxplot by type
  svg(paste(pname, "_ESPboxplot.svg", sep = ""), 10, 10, bg = "transparent")
  par(mar=c(9,5,2,2))
  boxplot(ESP~ClassLSR, data = din, las = 2, ylab = "ESP", cex.lab = 2, cex.axis = 1.5, col = vcol, ylim = c(0.4, 1))
  dev.off()
  
  # Boxplot Sheap
  svg(paste(pname, "_shapeboxplot.svg", sep = ""), 10, 10, bg = "transparent")
  par(mar=c(9,5,2,2))
  boxplot(shape~ClassLSR, data = din, las = 2, ylab = "Shape", cex.lab = 2, cex.axis = 1.5, col = vcol, ylim = c(0.5, 1))
  dev.off()
}



orderanddividedtable = function(din){
  lout = list()
  vorder = c("P", "B", "F","Cl", "Br", "Be", "NO2", "SO2", "S", "_CON", "COO", "onlyC", "C+O", "C+N", "C+O+N", "other")
  vordercycle = c("cycle-P", "cycle-B", "cycle-F", "cycle-Cl", "cycle-Br", "cycle-Be", "cycle-NO2", "cycle-SO2", "cycle-S", "cycle-CON", "cycle-COO", "cycle-onlyC", "cycle-C+O", "cycle-C+N", "cycle-C+O+N", "cycle-other")
  
  dout = c()
  for (order in vorder){
    for(i in seq(1, dim(din)[1])){
      if (order == din[i,1]){
        dout = rbind(dout, din[i,])
      }
    }
  }
  rownames(dout) = seq (1, dim(dout)[1])

  doutcycle = c()
  for (ordercycle in vordercycle){
    for(i in seq(1, dim(din)[1])){
      if (ordercycle == din[i,1]){
        doutcycle = rbind(doutcycle, din[i,])
      }
    }
  }
  rownames(doutcycle) = seq (1, dim(doutcycle)[1])
  
  lout = list (dout, doutcycle)
  
  return (lout)
}


###########
#  MAIN   #
###########

args = commandArgs(TRUE)
ptable = args[1]

ptable = "C://Users/Alexandre\ Borrel/Desktop/LSR/Pi_LSR_SheapClassif/PiLSRType"

dsheap = read.table(ptable, header = TRUE, sep = "\t")
print(colnames(dsheap))
#print (dsheap)

ldorded = orderanddividedtable (dsheap)

# write results
write.table(ldorded[[1]], file = paste(ptable, "_nocycle", sep = ""), sep = "\t")
write.table(ldorded[[2]], file = paste(ptable, "_cycle", sep = ""), sep = "\t")

dcycle = read.table (paste(ptable, "_cycle", sep = ""), sep = "\t", header = TRUE)
dnocycle = read.table (paste(ptable, "_nocycle", sep = ""), sep = "\t", header = TRUE)

# boxplot
boxplotFormat(dcycle, paste(ptable, "_cycle", sep = ""))
boxplotFormat(dnocycle, paste(ptable, "_nocycle", sep = ""))

# correlation
correlationPlot(dcycle, paste(ptable, "_cycle", sep = ""))
correlationPlot(dnocycle, paste(ptable, "_nocycle", sep = ""))