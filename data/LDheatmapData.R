LDheatmapData<-read.table("LDheatmapData.txt", header=TRUE)
HapMap.dat <- LDheatmapData[-1,]
for (LDindex in 1:length(HapMap.dat))
     HapMap.dat[,LDindex] <- genetics::as.genotype(HapMap.dat[,LDindex])
distance <- as.vector(as.matrix(LDheatmapData[1,]), mode="numeric")
rm(LDindex, LDheatmapData)
