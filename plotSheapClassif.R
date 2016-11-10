

vectCol = function(vin){
  
  vunique = unique(vin)
  lcol = seq(1, length(vunique))
  names(lcol) = vunique
  #print (lcol)
  
  vcol = NULL
  for(i in vin){
      vcol = append(vcol, lcol[i])      
  }
  return(list(vcol,lcol))
}





###########
#  MAIN   #
###########

args = commandArgs(TRUE)
ptable = args[1]

ptable = "C://Users/Alexandre\ Borrel/Desktop/LSR/PiLSRType"

dsheap = read.table(ptable, header = TRUE, sep = "\t")
print(colnames(dsheap))
# plot cross
vcol = vectCol(dsheap[,1])

png(paste(ptable, "_crossSheapScorestext.png", sep = ""), 800, 800)
par(mar=c(5,5,5,5))
plot(dsheap[,2],dsheap[,3], xlab = colnames(dsheap)[2], ylab = colnames(dsheap)[3], cex.lab = 2, type = "n")
text(dsheap[,2],dsheap[,3], labels = dsheap[,1], col = vcol[[1]], cex = 0.8)
legend("bottomright", legend = (names(vcol[[2]])), pch = 20, col = vcol[[2]])
dev.off()

png(paste(ptable, "_crossSheapScores.png", sep = ""), 800, 800)
par(mar=c(5,5,5,5))
plot(dsheap[,2],dsheap[,3], xlab = colnames(dsheap)[2], ylab = colnames(dsheap)[3], cex.lab = 2, pch = 19, col = vcol[[1]])
legend("bottomright", legend = (names(vcol[[2]])), pch = 20, col = vcol[[2]])
dev.off()

# Boxplot by type
png(paste(ptable, "_ESPboxplots.png", sep = ""), 800, 800)
par(mar=c(10,7,5,5))
boxplot(ESP~ClassLSR, data = dsheap, las = 2, ylab = "ESP", cex.lab = 2)
dev.off()

# Boxplot Sheap
png(paste(ptable, "_Shapeboxplots.png", sep = ""), 800, 800)
par(mar=c(10,7,5,5))
boxplot(shape~ClassLSR, data = dsheap, las = 2, ylab = "Shape", cex.lab = 2)
dev.off()



